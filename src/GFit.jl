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

export Domain, CartesianDomain, coords, axis, roi, int_tabulated, Measures,
    Model, @λ, SumReducer, select_reducer!, domain,
    MultiModel, patch!, evaluate!, isfixed, thaw, freeze, fit!


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

prepare!(comp::AbstractComponent, domain::AbstractDomain) =
    fill(NaN, length(domain))

dependencies(comp::AbstractComponent) = Symbol[]

evaluate!(buffer::Vector{Float64}, comp::AbstractComponent, domain::AbstractDomain, pars...) =
    error("No evaluate!() method implemented for $T")

# Built-in components
include("components/SimplePar.jl")
include("components/FuncWrap.jl")
include("components/OffsetSlope.jl")
include("components/Gaussian.jl")
include("components/Lorentzian.jl")


# ====================================================================
# CompEval: a wrapper for a component evaluated on a specific domain
#
mutable struct CompEval{TComp <: AbstractComponent, TDomain <: AbstractDomain}
    comp::TComp
    domain::TDomain
    counter::Int
    dependencies::Vector{Vector{Float64}}
    lastvalues::Vector{Float64}
    buffer::Vector{Float64}
    cfixed::Bool

    function CompEval(_comp::AbstractComponent, domain::AbstractDomain)
        # Components internal state may be affected by `prepare!`
        # call.  Avoid overwriting input state with a deep copy.
        comp = deepcopy(_comp)
        buffer = prepare!(comp, domain)
        return new{typeof(comp), typeof(domain)}(
            comp, domain, 0,
            Vector{Vector{Float64}}(),
            fill(NaN, length(getparams(comp))),
            buffer, false)
    end
end


function evaluate!(c::CompEval, pvalues::Vector{Float64})
    @assert length(getparams(c.comp)) == length(pvalues)

    # Do we actually need a new evaluation?
    if (any(c.lastvalues .!= pvalues)  ||  (c.counter == 0)  ||  (length(c.dependencies) > 0))
        c.lastvalues .= pvalues
        c.counter += 1
        if !all(.!isnan.(pvalues))
            println("One or more parameter value(s) are NaN:")
            println(pvalues)
            @assert all(.!isnan.(pvalues))
        end
        evaluate!(c.buffer, c.comp, c.domain, c.dependencies..., pvalues...)
    end
    return c.buffer
end
evaluate!(c::CompEval) = evaluate!(c, getfield.(values(getparams(c.comp)), :val))


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
                out[arg[1]] = arg[2]
            end
            return parse_args(out)
        end

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


function evaluate!(model::Model)
    if !isnothing(model.parent)
        evaluate!(model.parent)
    else
        eval1!(model)
        eval2!(model)
        eval3!(model)
        eval4!(model)
    end
    return model
end

function eval1!(model::Model)
    empty!(model.params)
    empty!(model.pvalues)
    empty!(model.patched)
    empty!(model.ifree)
    empty!(model.buffers)

    for (cname, ceval) in model.cevals
        for (pname, par) in getparams(ceval.comp)
            if !(par.low <= par.val <= par.high)
                s = "Value outside limits for param [$(cname)].$(pname):\n" * string(par)
                error(s)
            end
            model.params[ cname][pname] = par
            model.pvalues[cname][pname] = par.val
            model.patched[cname][pname] = par.val
            if !par.fixed  &&  (ceval.cfixed == 0)  &&  !isa(par.patch, Symbol)
                push!(model.ifree, length(internal_data(model.params)))
            end
        end
        model.buffers[cname] = ceval.buffer
    end
end

function eval2!(model::Model; fromparent=false)
    # Evaluate patched parameters
    if !isnothing(model.parent)  &&  !fromparent
        eval2!(model.parent)
    else
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

function eval3!(model::Model)
    # Evaluate model
    for (cname, ceval) in model.cevals
        evaluate!(ceval, values(model.patched[cname]))
    end
end

function eval4!(model::Model, unc::Union{Nothing, Vector{Float64}}=nothing)
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


function setindex!(model::Model, comp::AbstractComponent, cname::Symbol)
    ceval = CompEval(comp, model.domain)
    model.cevals[cname] = ceval
    (length(model.maincomp) == 0)  &&  push!(model.maincomp, cname)
    evaluate!(model)
end


function isfixed(model::Model, cname::Symbol)
    @assert cname in keys(model.cevals) "Component $cname is not defined"
    return model.cevals[cname].cfixed
end

function freeze(model::Model, cname::Symbol)
    @assert cname in keys(model.cevals) "Component $cname is not defined"
    model.cevals[cname].cfixed = true
    evaluate!(model)
    nothing
end

function thaw(model::Model, cname::Symbol)
    @assert cname in keys(model.cevals) "Component $cname is not defined"
    model.cevals[cname].cfixed = false
    evaluate!(model)
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
