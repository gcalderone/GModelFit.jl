module GFit

using Printf, PrettyTables
using Statistics, Distributions
using DataStructures
using LsqFit
using MacroTools
using Dates
using ProgressMeter
using Random

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

export AbstractDomain, Domain, CartesianDomain, coords, axis, Measures, uncerts,
    Model, @λ, select_maincomp!, SumReducer, domain,
    MultiModel, evaluate, isfreezed, thaw!, freeze!, fit!

include("HashVector.jl")
include("domain.jl")

# ====================================================================
"""
    λFunct

A representation for a λ function, containing a reference to the function itself, a string representation of its source code (for displaying purposes) and a list of its arguments.
"""
struct λFunct
    funct::Function
    display::String
    args::Vector{Symbol}   # positional arguments
    optargs::Vector{Expr}  # optional arguments with default values
end
(f::λFunct)(args...; kws...) = f.funct(args...; kws...)

"""
    @λ expr

Macro to generate `λFunct` objects using the same syntax as in a standard Julia anonymous function.
"""
macro λ(_expr)
    @assert isexpr(longdef(_expr), :function)
    expr = prettify(_expr)
    def  = splitdef(expr)
    args    = convert(Vector{Symbol}, filter(x -> isa(x, Symbol), def[:args]))
    optargs = convert(Vector{Expr}  , filter(x -> isa(x, Expr)  , def[:args]))
    return esc(:(GFit.λFunct($expr, string($(QuoteNode(expr))), $args, $optargs)))
end


# ====================================================================
"""
    Parameter

A structure representing a model parameter.

# Fields:
 - `val::Float64`: parameter value (either initial guess when set before fitting, or best fit one after fitting);
 - `low::Float64`: lower limit for the value (default: `-Inf`);
 - `high::Float64`: upper limit for the value (default: `+Inf`);
 - `fixed::Bool`: whether the parameter is fixed during fitting (default: `false`);
 - `patch::Union{Nothing, Symbol, λFunct}`: patch prescription within the same model;
 - `mpatch::Union{Nothing, λFunct}`: patch prescription in a multi-model analysis;
 - `actual::Float64`: actual value for the parameter (i.e. after applying the patch prescription)`;
 - `unc::Float64`: 1σ uncertainty associated to the parameter value.

Note: the `Parameter` fields are supposed to be accessed directly by the user, without invoking any get/set method.
"""
mutable struct Parameter
    val::Float64
    low::Float64              # lower limit value
    high::Float64             # upper limit value
    fixed::Bool
    patch::Union{Nothing, Symbol, λFunct}
    mpatch::Union{Nothing, λFunct}
    actual::Float64
    unc::Float64
    Parameter(value::Number) = new(float(value), -Inf, +Inf, false, nothing, nothing, NaN, NaN)
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

dependencies(comp::AbstractComponent) = Symbol[]

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
    updated::Bool

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
    c.updated  &&  return

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
    c.updated = true
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

"""
    Model

A structure containing a model description, whose evaluation is suitable to be compared to a single empirical dataset.

Constructor is: `Model(domain::AbstractDomain, components...)`
where the first argument is either a `Domain` or `CartesianDomain` object, and the remaining one(s) is (are) the model component(s), which may be given as:
- a single `Dict{Symbol, AbstractComponent}`, where the keys are the names and the values the component objects;
- a single component, which will have a default name is assigned (`:main`);
- a single `λFunct`, which will be wrapped into an `LComp` component and a default name will be assigned (`:main`);
- one or more `Pair{Symbol, AbstractComponent}`, where the first element is the name and the second is the component.

You may access the individual component in a `Model` using the indexing syntax, as if it was a `Dict{Symbol, AbstractComponent}`.  Also, you may add new components to a `Model` after it has been created using the same synatx.  Finally, you may use the `keys()` and `haskey()` functions with their usual meanings.

Individual components may be *freezed* (i.e. have all its parameters fixed during fitting, despite the individual `Parameter` settings) or *thawed* using the `freeze!()` and `thaw!()` functions.  Use the `isfreezed()` function to check if a component is freezed.

The main component, i.e. the one whose evaluation corresponds to the overall model evaluation, is typically automatically identified by analyzing the component dependencies.  However a specific component may be forced to be the main one by invoking `select_maincomp!`.

The model is automatically evaluated whenever needed, however there are a few cases where it is not possible to trigger an automatic evaluation, e.g. immediately after the user modifies a `Parameter` value. In this case an evaluation can be forced by invoking `evaluate()`.

The most important function for a `Model` object is `fit!()`, which allows to fit the model against an empirical dataset. The `!` in the name reminds us that after fitting the new parameter values will be set to the best fit ones (rather than retaining their original values).

The model and all component evaluation can be obtained by using the `Model` object has if it was a function: with no arguments it will return the main component evaluation, while if a `Symbol` is given as argument it will return the evaluation of the component with the same name.
"""
mutable struct Model   # mutable because of parent and maincomp
    parent::Union{Nothing, AbstractMultiModel}
    domain::AbstractDomain
    cevals::OrderedDict{Symbol, CompEval}
    params::HashHashVector{Parameter}
    pvalues::HashHashVector{Float64}
    actual::HashHashVector{Float64}
    ifree::Vector{Int}
    buffers::OrderedDict{Symbol, Vector{Float64}}
    maincomp::Symbol

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
                    Symbol(""))

        for (name, item) in parse_args(args...)
            model[name] = item
        end
        return model
    end
end


function find_maincomp(model::Model)
    if model.maincomp != Symbol("")
        return model.maincomp
    end

    if length(model.cevals) == 1
        return collect(keys(model.cevals))[1]
    end

    maincomps = collect(keys(model.cevals))
    for (cname, ceval) in model.cevals
        for d in dependencies(ceval.comp)
            @assert d in maincomps "$cname depends on $d, but the latter is not a component in the model."
            i = findfirst(maincomps .== d)
            deleteat!(maincomps, i)
        end
    end

    if length(maincomps) > 1
        # Ignoring components with no dependencies
        for (cname, ceval) in model.cevals
            if length(dependencies(ceval.comp)) == 0
                i = findfirst(maincomps .== cname)
                if !isnothing(i)
                    deleteat!(maincomps, i)

                    # ...but keep the last
                    if length(maincomps) == 1
                        return maincomps[1]
                    end
                end
            end
        end
    end

    return maincomps[end]
end


"""
    evaluate(model::Model)

Evaluate a `Model` and update internal structures.
"""
function evaluate(model::Model)
    eval_step0(model)
    # eval_step1()
    eval_step2(model)
    eval_step3(model)
    eval_step4(model)
    return model
end


# Evaluation step 0: update internal structures before fitting
function eval_step0(model::Model)
    empty!(model.params)
    empty!(model.pvalues)
    empty!(model.actual)
    empty!(model.ifree)

    ipar = 1
    for (cname, ceval) in model.cevals
        for (pname, _par) in getparams(ceval.comp)
            # Parameter may be changed here, hence we take a copy of the original one
            par = deepcopy(_par)
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
            model.actual[cname][pname] = par.val

            if !isnothing(par.patch)
                @assert isnothing(par.mpatch) "Parameter [$cname].$pname has both patch and mpatch fields set, while only one is allowed"
                if isa(par.patch, Symbol)  # use same param. value from a different component
                    par.fixed = true
                else                       # invoke a patch function
                    if length(par.patch.args) == 1
                        par.fixed = true
                    else
                        par.fixed = false
                    end
                end
            elseif !isnothing(par.mpatch)
                @assert !isnothing(model.parent) "Parameter [$cname].$pname has the mpatch field set but no MultiModel has been created"
                if length(par.mpatch.args) == 1
                    par.fixed = true
                else
                    par.fixed = false
                end
            end
            if ceval.cfixed != 0
                par.fixed = true
            end
            if !par.fixed
                push!(model.ifree, ipar)
            end
            ipar += 1
        end

        empty!(ceval.deps)
        for d in dependencies(ceval.comp)
            push!(ceval.deps, model.buffers[d])
        end
    end
end


# Evaluation step 1: set new model parameters
function eval_step1(model::Model, pvalues::Vector{Float64})
    internal_data(model.pvalues)[model.ifree] .= pvalues
end


# Evaluation step 2: copy all parameter values into actual, then
# update the latter by invoking the user patch functions.
function eval_step2(model::Model)
    # Reset `updated` flag
    for (cname, ceval) in model.cevals
        ceval.updated = false
    end
    # Copy pvalues into actual
    internal_data(model.actual) .= internal_data(model.pvalues)
    # Patch parameter values
    for (cname, hv) in model.params
        for (pname, par) in hv
            if !isnothing(par.patch)
                @assert isnothing(par.mpatch) "Parameter [$cname].$pname has both patch and mpatch fields set, while only one is allowed"
                if isa(par.patch, Symbol)  # use same param. value from a different component
                    model.actual[cname][pname] = model.pvalues[par.patch][pname]
                else                       # invoke a patch function
                    if length(par.patch.args) == 1
                        model.actual[cname][pname] = par.patch(model.pvalues)
                    else
                        model.actual[cname][pname] = par.patch(model.pvalues[cname][pname], model.pvalues)
                    end
                end
            elseif !isnothing(par.mpatch)
                @assert !isnothing(model.parent) "Parameter [$cname].$pname has the mpatch field set but no MultiModel has been created"
                if length(par.mpatch.args) == 1
                    model.actual[cname][pname] = par.mpatch(model.parent.pvalues)
                else
                    model.actual[cname][pname] = par.mpatch(model.pvalues[cname][pname], model.parent.pvalues)
                end
            end
        end
    end
end


# Evaluation step 3: actual evaluation of model components, starting
# from the main one and following dependencies
eval_step3(model::Model) = eval_step3(model, find_maincomp(model))
function eval_step3(model::Model, cname::Symbol)
    # Recursive evaluation of dependencies
    for d in dependencies(model.cevals[cname].comp)
        eval_step3(model, d)
    end
    evaluate!(model.cevals[cname], values(model.actual[cname]))
end


# Evaluation step 4: copy back bestfit and actual values, as well as
# uncertainties, into their original Parameter structures.
function eval_step4(model::Model, uncerts=Vector{Float64}[])
    ipar = 1
    i = 1
    for (cname, hv) in model.params
        for (pname, par) in hv
            par.val    = model.pvalues[cname][pname]
            par.actual = model.actual[ cname][pname]
            if (length(uncerts) > 0)  &&  (ipar in model.ifree)
                par.unc = uncerts[i]
                i += 1
            else
                par.unc = NaN
            end
            ipar += 1
        end
    end

    # Also update Model's parameters
    for (cname, ceval) in model.cevals
        for (pname, par) in getparams(ceval.comp)
            par.val    = model.params[cname][pname].val
            par.actual = model.params[cname][pname].actual
            par.unc = NaN
        end
    end
end


# User interface
setindex!(model::Model, v::Real, cname::Symbol) = setindex!(model, SimplePar(v), cname)
setindex!(model::Model, f::λFunct, cname::Symbol) = setindex!(model, λComp(f), cname)
function setindex!(model::Model, comp::AbstractComponent, cname::Symbol)
    ceval = CompEval(comp, model.domain)
    model.cevals[cname] = ceval
    model.buffers[cname] = ceval.buffer
    evaluate(model)
end

free_params(model::Model) = internal_data(model.params)[model.ifree]

"""
    isfreezed(model::Model, cname::Symbol)

Check whether a component is *freezed* in the model.
"""
function isfreezed(model::Model, cname::Symbol)
    @assert cname in keys(model.cevals) "Component $cname is not defined"
    return model.cevals[cname].cfixed
end

"""
    freeze!(model::Model, cname::Symbol)

Freeze a component in the model (i.e. treat all component parameters as fixed for fitting).
"""
function freeze!(model::Model, cname::Symbol)
    @assert cname in keys(model.cevals) "Component $cname is not defined"
    model.cevals[cname].cfixed = true
    evaluate(model)
    nothing
end

"""
    thaw!(model::Model, cname::Symbol)

Thaw a freezed component in the model (i.e. treat component parameters as fixed only if explicitly set in the corresponding `Parameter` structure).
"""
function thaw!(model::Model, cname::Symbol)
    @assert cname in keys(model.cevals) "Component $cname is not defined"
    model.cevals[cname].cfixed = false
    evaluate(model)
    nothing
end


Base.keys(p::Model) = collect(keys(p.cevals))

"""
    haskey(m::Model, name::Symbol)

Check whether a component exists in model.
"""
Base.haskey(m::Model, name::Symbol) = haskey(m.cevals, name)
function Base.getindex(model::Model, name::Symbol)
    if name in keys(model.cevals)
        return model.cevals[name].comp
    end
    error("Name $name not defined")
end

"""
    domain(model::Model)

Return the domain where the model is evaluated.
"""
domain(model::Model) = model.domain
(model::Model)() = model.cevals[find_maincomp(model)].buffer
(model::Model)(name::Symbol) = model.cevals[name].buffer

"""
    select_maincomp!(model::Model, cname::Symbol)

Force a component to be the final one for model evaluation.
"""
function select_maincomp!(model::Model, cname::Symbol)
    @assert haskey(model, cname)
    model.maincomp = cname
end

include("multimodel.jl")

abstract type AbstractFitProblem end
include("minimizers.jl")
include("fit.jl")

include("show.jl")
include("utils.jl")

end
