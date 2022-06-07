module GFit

using Printf, PrettyTables
using Statistics, Distributions
using DataStructures
using LsqFit
using MacroTools
using Dates
using ProgressMeter

import Base.show
import Base.ndims
import Base.size
import Base.length
import Base.haskey
import Base.keys
import Base.getindex
import Base.setindex!
import Base.reshape
import Base.propertynames
import Base.getproperty
import Base.setproperty!
import Base.iterate
import Base.push!

export Domain, CartesianDomain, coords, axis, roi, Measures,
    Model, @λ, SumReducer, domain,
    MultiModel, patch!, evaluate, isfixed, thaw, freeze, fit!


include("HashVector.jl")
include("domain.jl")

# ====================================================================
struct λFunct
    funct::Function
    display::String
    args::Vector{Symbol}   # positional arguments
    optargs::Vector{Expr}  # optional arguments with default values
end
(f::λFunct)(args...; kws...) = f.funct(args...; kws...)

macro λ(_expr)
    @assert isexpr(longdef(_expr), :function)
    expr = prettify(_expr)
    def  = splitdef(expr)
    args    = convert(Vector{Symbol}, filter(x -> isa(x, Symbol), def[:args]))
    optargs = convert(Vector{Expr}  , filter(x -> isa(x, Expr)  , def[:args]))
    return esc(:(GFit.λFunct($expr, string($(QuoteNode(expr))), $args, $optargs)))
end


# ====================================================================
# Parameter
mutable struct Parameter
    val::Float64
    low::Float64              # lower limit value
    high::Float64             # upper limit value
    step::Float64
    fixed::Bool
    patch::Union{Nothing, Symbol, λFunct}
    superpatch::Union{Nothing, λFunct}
    pval::Float64
    unc::Float64
    Parameter(value::Number) = new(float(value), -Inf, +Inf, NaN, false, nothing, nothing, NaN, NaN)
end


# ====================================================================
# Components:
#
# A *component* is a generic implementation of a building block for a
# model. It must inherit `AbstractComponent` and implement the
# `evaluate!` method (optionally also `prepare!`).  The structure may
# contain zero or more field of type Parameter (see above)
abstract type AbstractComponent end

# Fall back methods
function getparams(comp::AbstractComponent)
    out = OrderedDict{Symbol, Parameter}()
    for name in fieldnames(typeof(comp))
        par = getfield(comp, name)
        if isa(par, Parameter)
            out[name] = par
        end
    end
    return out
end

deps(comp::AbstractComponent) = Symbol[]

prepare!(comp::AbstractComponent, domain::AbstractDomain) =
    fill(NaN, length(domain))

evaluate!(buffer::Vector{Float64}, comp::T, domain::AbstractDomain, pars...) where T <: AbstractComponent=
    error("No evaluate!() method implemented for $T")

# Built-in components
include("components/SimplePar.jl")
include("components/LComp.jl")
include("components/OffsetSlope.jl")
include("components/Gaussian.jl")
include("components/Lorentzian.jl")
include("components/SumReducer.jl")


# ====================================================================
# CompEval: a wrapper for a component evaluated on a specific domain
#
mutable struct CompEval{TComp <: AbstractComponent, TDomain <: AbstractDomain}
    comp::TComp
    domain::TDomain
    counter::Int
    deps::Vector{Vector{Float64}}
    lastvalues::Vector{Float64}
    buffer::Vector{Float64}
    cfixed::Bool
    done::Bool

    function CompEval(_comp::AbstractComponent, domain::AbstractDomain)
        # Components internal state may be affected by `prepare!`
        # call.  Avoid overwriting input state with a deep copy.
        comp = deepcopy(_comp)
        buffer = prepare!(comp, domain)
        return new{typeof(comp), typeof(domain)}(
            comp, domain, 0,
            Vector{Vector{Float64}}(),
            fill(NaN, length(getparams(comp))),
            buffer, false, false)
    end
end


function evaluate!(c::CompEval, pvalues::Vector{Float64})
    c.done  &&  (return c.buffer)

    # Do we actually need a new evaluation?
    if (any(c.lastvalues .!= pvalues)  ||  (c.counter == 0)  ||  (length(c.deps) > 0))
        if length(c.deps) > 0
            evaluate!(c.buffer, c.comp, c.domain, c.deps, pvalues...)
        else
            evaluate!(c.buffer, c.comp, c.domain, pvalues...)
        end
        c.lastvalues .= pvalues
        c.counter += 1
    end
    c.done = true
    return c.buffer
end
evaluate!(c::CompEval) = evaluate!(c, getfield.(values(getparams(c.comp)), :val))


# Facility to easily evaluate a component
function evaluate(domain::AbstractDomain, comp::AbstractComponent)
    ceval = CompEval(comp, domain)
    evaluate(ceval)
    return ceval.buffer
end


# ====================================================================
# Model
#
abstract type AbstractMultiModel end

# A model prediction suitable to be compared to a single empirical dataset
struct Model
    parent::Union{Nothing, AbstractMultiModel}
    domain::AbstractDomain
    cevals::OrderedDict{Symbol, CompEval}
    params::HashHashVector{Parameter}
    pvalues::HashHashVector{Float64}
    patched::HashHashVector{Float64}
    ifree::Vector{Int}
    buffers::OrderedDict{Symbol, Vector{Float64}}
    maincomp::Vector{Symbol}  # it is a vector to avoid making Model mutable

    function Model(domain::AbstractDomain, args...)
        function parse_args(args::AbstractDict)
            out = OrderedDict{Symbol, AbstractComponent}()
            for (name, item) in args
                isa(item, Number)  &&  (item = SimplePar(item))
                @assert isa(name, Symbol)
                @assert isa(item, AbstractComponent)
                out[name] = item
            end
            return out
        end

        function parse_args(args::Vararg{Pair})
            out = OrderedDict{Symbol, AbstractComponent}()
            for arg in args
                if isa(arg[2], AbstractComponent)
                    out[arg[1]] = arg[2]
                elseif isa(arg[2], λFunct)
                    out[arg[1]] = λComp(arg[2])
                elseif isa(arg[2], Number)
                    out[arg[1]] = SimplePar(arg[2])
                else
                    error("Unsupported data type: " * string(typeof(arg[2])) *
                          ".  Must be an AbstractComponent, a λFunct or a real number.")
                end
            end
            return parse_args(out)
        end

        parse_args(arg::AbstractComponent) = parse_args(:main => arg)
        parse_args(arg::λFunct) = parse_args(:main => λComp(arg))
        parse_args(arg::Real) = parse_args(:main => SimplePar(arg))

        model = new(nothing, domain, OrderedDict{Symbol, CompEval}(),
                    HashHashVector{Parameter}(),
                    HashHashVector{Float64}(),
                    HashHashVector{Float64}(),
                    Vector{Int}(),
                    OrderedDict{Symbol, Vector{Float64}}(),
                    Vector{Symbol}())

        for (name, item) in parse_args(args...)
            model[name] = item
        end
        return model
    end
end


function evaluate(model::Model)
    if !isnothing(model.parent)
        evaluate(model.parent)
    else
        eval_step1(model)
        eval_step2(model)
        eval_step3(model)
        eval_step4(model)
    end
    return model
end

function eval_step1(model::Model)
    empty!(model.params)
    empty!(model.pvalues)
    empty!(model.patched)
    empty!(model.ifree)

    for (cname, ceval) in model.cevals
        for (pname, par) in getparams(ceval.comp)
            if !(par.low <= par.val <= par.high)
                s = "Value outside limits for param [$(cname)].$(pname):\n" * string(par)
                error(s)
            end
            if isnan(par.low)  ||  isnan(par.high)  ||  isnan(par.val)
                s = "NaN value detected for param [$(cname)].$(pname):\n" * string(par)
                error(s)
            end
            model.params[ cname][pname] = par
            model.pvalues[cname][pname] = par.val
            model.patched[cname][pname] = par.val
            if !par.fixed  &&  (ceval.cfixed == 0)  &&  !isa(par.patch, Symbol)
                push!(model.ifree, length(internal_data(model.params)))
            end
        end

        empty!(ceval.deps)
        for d in deps(ceval.comp)
            push!(ceval.deps, model.buffers[d])
        end
    end
end

function eval_step2(model::Model; fromparent=false)
    if !isnothing(model.parent)  &&  !fromparent
        eval_step2(model.parent)
    else
        # Copy pvalues into patched
        internal_data(model.patched) .= internal_data(model.pvalues)
        # Reset `done` flag
        for (cname, ceval) in model.cevals
            ceval.done = false
        end
        # Patch parameter values
        for (cname, hv) in model.params
            for (pname, par) in hv
                if !isnothing(par.patch)
                    if isa(par.patch, Symbol)
                        # Use same parameter from a different component
                        model.patched[cname][pname] = model.pvalues[par.patch][pname]
                    else
                        # Evaluate a patch function
                        model.patched[cname][pname] = par.patch(model.pvalues[cname][pname], model.pvalues)
                    end
                end
            end
        end
    end
end

# Evaluation of model component, starting from model.maincomp[1]
eval_step3(model::Model) = eval_step3(model, model.maincomp[1])

function eval_step3(model::Model, cname::Symbol)
    # Recursive evaluation of dependencies
    for d in deps(model.cevals[cname].comp)
        eval_step3(model, d)
    end
    # Evaluate current component
    evaluate!(model.cevals[cname], values(model.patched[cname]))
end

function eval_step4(model::Model, unc::Union{Nothing, Vector{Float64}}=nothing)
    # Update values, uncertainties and patched params from ModelEval
    # to the actual Model structure.
    i = 1
    for (cname, hv) in model.params
        for (pname, par) in hv
            par.val  = model.pvalues[cname][pname]
            par.pval = model.patched[cname][pname]
            if !isnothing(unc)  &&  (i in model.ifree)
                par.unc = unc[i]
            end
            i += 1
        end
    end
end


setindex!(model::Model, v::Real, cname::Symbol) = setindex!(model, SimplePar(v), cname)
setindex!(model::Model, f::λFunct, cname::Symbol) = setindex!(model, λComp(f), cname)
function setindex!(model::Model, comp::AbstractComponent, cname::Symbol)
    ceval = CompEval(comp, model.domain)
    model.cevals[cname] = ceval
    model.buffers[cname] = ceval.buffer

    # Last added component is the one to be used as main model
    if length(model.maincomp) == 0
        push!(model.maincomp, cname)
    else
        model.maincomp[1] = cname
    end
    evaluate(model)
end


function isfixed(model::Model, cname::Symbol)
    @assert cname in keys(model.cevals) "Component $cname is not defined"
    return model.cevals[cname].cfixed
end

function freeze(model::Model, cname::Symbol)
    @assert cname in keys(model.cevals) "Component $cname is not defined"
    model.cevals[cname].cfixed = true
    evaluate(model)
    nothing
end

function thaw(model::Model, cname::Symbol)
    @assert cname in keys(model.cevals) "Component $cname is not defined"
    model.cevals[cname].cfixed = false
    evaluate(model)
    nothing
end


Base.keys(p::Model) = collect(keys(p.cevals))
Base.haskey(p::Model, name::Symbol) = haskey(p.cevals, name)
function Base.getindex(model::Model, name::Symbol)
    if name in keys(model.cevals)
        return model.cevals[name].comp
    elseif name in keys(model.revals)
        return model.revals[name].red
    end
    error("Name $name not defined")
end
domain(model::Model) = model.domain
(model::Model)() = model.cevals[model.maincomp[1]].buffer
(model::Model)(name::Symbol) = model.cevals[name].buffer


include("multimodel.jl")
include("minimizers.jl")
include("fit.jl")
include("show.jl")

end
