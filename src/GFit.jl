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
include("utils.jl")

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
# Parameter and component identifiers
struct CompParamID
    cname::Symbol
    param::Symbol
end

mutable struct Parameter
    val::Float64
    low::Float64              # lower limit value
    high::Float64             # upper limit value
    step::Float64
    fixed::Bool
    patched::Float64
    unc::Float64
    Parameter(value::Number) = new(float(value), -Inf, +Inf, NaN, false, NaN, NaN)
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

# dependencies(comp::AbstractComponent) = Symbol[]

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
    params::OrderedDict{Symbol, Parameter}
    counter::Int
    # dependencies::Vector{Vector{Float64}}
    lastvalues::Vector{Float64}
    buffer::Vector{Float64}
    cfixed::Bool

    function CompEval(_comp::AbstractComponent, domain::AbstractDomain)
        # Components internal state may be affected by `prepare!`
        # call.  Avoid overwriting input state with a deep copy.
        comp = deepcopy(_comp)

        params = getparams(comp)
        buffer = prepare!(comp, domain)
        return new{typeof(comp), typeof(domain)}(
            comp, domain, params, 0,
            fill(NaN, length(params)),
            buffer, false)
    end
end


function evaluate!(c::CompEval, pvalues::Vector{Float64})
    @assert length(c.params) == length(pvalues)

    # Do we actually need a new evaluation?
    if (any(c.lastvalues .!= pvalues)  ||  (c.counter == 0))
        c.lastvalues .= pvalues
        c.counter += 1
        if !all(.!isnan.(pvalues))
            println("One or more parameter value(s) are NaN:")
            println(pvalues)
            @assert all(.!isnan.(pvalues))
        end
        evaluate!(c.buffer, c.comp, c.domain, pvalues...)
    end
    return c.buffer
end
evaluate!(c::CompEval) = evaluate!(c, getfield.(values(c.params), :val))



# ====================================================================
# Reducer
#

abstract type AbstractReducer end

prepare!(comp::AbstractReducer, domain::AbstractDomain) = fill(NaN, length(domain))

struct λReducer <: AbstractReducer
    f::λFunct
    λReducer(f::λFunct) = new(f)
end

function evaluate!(buffer::Vector{Float64}, red::λReducer,
                   domain::AbstractDomain, m::OrderedDict{Symbol, Vector{Float64}}, pars...)
    buffer .= red.f.funct(domain, m, pars...)
    nothing
end

struct SumReducer <: AbstractReducer
    list::Vector{Symbol}
end


function evaluate!(buffer::Vector{Float64}, red::SumReducer,
                   domain::AbstractDomain, args::OrderedDict{Symbol, Vector{Float64}})
    buffer .= 0.
    for name in red.list
        buffer .+= args[name]
    end
    nothing
end

mutable struct ReducerEval{T <: AbstractReducer}
    red::T
    counter::Int
    buffer::Vector{Float64}
end


# ====================================================================
# Model
#
const PatchComp = HashVector{Float64}

struct ModelEval
    params::Vector{Parameter}
    ifree::Vector{Int}
    pvalues::Vector{Float64}
    patched::Vector{Float64}
    patchcomps::OrderedDict{Symbol, PatchComp}
    reducer_args::OrderedDict{Symbol, Vector{Float64}}
end

abstract type AbstractMultiModel end

# A model prediction suitable to be compared to a single empirical dataset
mutable struct Model
    parent::Union{Nothing, AbstractMultiModel}
    domain::AbstractDomain
    cevals::OrderedDict{Symbol, CompEval}
    revals::OrderedDict{Symbol, ReducerEval}
    rsel::Symbol
    patchfuncts::Vector{λFunct}
    meval::Union{Nothing, ModelEval}

    function Model(domain::AbstractDomain, args...)
        function parse_args(args::AbstractDict)
            out = OrderedDict{Symbol, Union{AbstractComponent, AbstractReducer}}()
            for (name, item) in args
                isa(item, Number)  &&  (item = SimplePar(item))
                @assert isa(name, Symbol)
                @assert isa(item, AbstractComponent) || isa(item, AbstractReducer)
                out[name] = item
            end
            return out
        end

        function parse_args(args::Vararg{Pair})
            out = OrderedDict{Symbol, Union{AbstractComponent, AbstractReducer}}()
            for arg in args
                out[arg[1]] = arg[2]
            end
            return parse_args(out)
        end

        model = new(nothing, domain,
                    OrderedDict{Symbol, CompEval}(),
                    OrderedDict{Symbol, ReducerEval}(),
                    Symbol(""),
                    Vector{λFunct}(),
                    nothing)
        evaluate!(model)  # populate meval
        for (name, item) in parse_args(args...)
            model[name] = item
        end
        if length(model.revals) == 0
            model[:default_sum] = SumReducer(collect(keys(model.cevals)))
        end
        return model
    end
end


function ModelEval(model::Model)
    params = Vector{Parameter}()
    ifree = Vector{Int}()
    pvalues = Vector{Float64}()
    patched = Vector{Float64}()
    patchcomps = OrderedDict{Symbol, PatchComp}()
    reducer_args = OrderedDict{Symbol, Vector{Float64}}()

    i = 1
    for (cname, ceval) in model.cevals
        patchcomps[cname] = HashVector{Float64}(patched)
        for (pname, par) in ceval.params
            if !(par.low <= par.val <= par.high)
                s = "Value outside limits for param [$(cname)].$(pname):\n" * string(par)
                error(s)
            end
            if (!par.fixed)  &&  (model.cevals[cname].cfixed == 0)
                push!(ifree, i)
            end
            push!(params, par)
            push!(pvalues, par.val)
            push!(patchcomps[cname], pname, NaN)
            i += 1
        end
    end
    @assert length(patched) == length(pvalues)

    for (cname, ceval) in model.cevals
        reducer_args[cname] = ceval.buffer
    end
    for (rname, reval) in model.revals
        reducer_args[rname] = reval.buffer
    end

    return ModelEval(params, ifree, pvalues, patched, patchcomps, reducer_args)
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
    # Update ModelEval structure
    model.meval = ModelEval(model)
end

function eval2!(model::Model; fromparent=false)
    # Evaluate patched parameters
    if !isnothing(model.parent)  &&  !fromparent
        eval2!(model.parent)
    else
        model.meval.patched .= model.meval.pvalues  # copy all values by default
        for pf in model.patchfuncts
            pf.funct(model.meval.patchcomps)
        end
    end
end

function eval3!(model::Model)
    # Evaluate model
    i1 = 1
    for (cname, ceval) in model.cevals
        if length(ceval.params) > 0
            i2 = i1 + length(ceval.params) - 1
            evaluate!(ceval, model.meval.patched[i1:i2])
            i1 += length(ceval.params)
        else
            evaluate!(ceval, Float64[])
        end
    end

    for (rname, reval) in model.revals
        (rname == model.rsel)  &&  continue  # this should be evaluated as last one
        reval.counter += 1
        evaluate!(reval.buffer, reval.red, model.domain, model.meval.reducer_args)
    end

    if length(model.revals) > 0
        reval = model.revals[model.rsel]
        reval.counter += 1
        evaluate!(reval.buffer, reval.red, model.domain, model.meval.reducer_args)
    end
end

function eval4!(model::Model, unc::Union{Nothing, Vector{Float64}}=nothing)
    # Update values, uncertainties and patched params from ModelEval
    # to the actual Model structure.
    i = 1
    j = 1
    for (cname, ceval) in model.cevals
        for (pname, par) in ceval.params
            if i in model.meval.ifree
                if !isnothing(unc)
                    par.unc = unc[j]
                elseif par.val != model.meval.pvalues[i]
                    # overwrite only if param. value has changed
                    # (i.e. avoid loosing uncertainty when invoking
                    # evaluate!(model)
                    par.unc = NaN
                end
                j += 1
            else
                par.unc = NaN
            end
            par.val = model.meval.pvalues[i]
            par.patched = model.meval.patched[i]
            i += 1
        end
    end
end


function setindex!(model::Model, comp::AbstractComponent, cname::Symbol)
    @assert !haskey(model.revals, cname)  "Name $cname already exists as a reducer name"
    ceval = CompEval(deepcopy(comp), model.domain)
    model.cevals[cname] = ceval
    evaluate!(model)
end

function setindex!(model::Model, reducer::T, rname::Symbol) where T <: AbstractReducer
    @assert !haskey(model.cevals, rname) "Name $rname already exists as a component name"
    model.revals[rname] = ReducerEval{T}(reducer, 1, prepare!(reducer, model.domain))
    (model.rsel == Symbol(""))  &&  (model.rsel = rname)
    evaluate!(model)
    return model
end


function select_reducer!(model::Model, rname::Symbol)
    @assert haskey(model.revals, rname) "$rname is not a reducer name"
    model.rsel = rname
end


function geteval(model::Model, name::Symbol)
    if haskey(model.cevals, name)
        return model.cevals[name].buffer
    else
        if haskey(model.revals, name)
            return model.revals[name].buffer
        else
            return Float64[]
        end
    end
end

geteval(model::Model) = geteval(model, model.rsel)


function patch!(func::λFunct, model::Model)
    push!(model.patchfuncts, func)
    evaluate!(model)
    return model
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
(model::Model)() = geteval(model)
(model::Model)(name::Symbol) = geteval(model, name)


include("multimodel.jl")
include("minimizers.jl")
include("fit.jl")
include("show.jl")

end
