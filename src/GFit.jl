module GFit

using Printf, PrettyTables
using Statistics, Distributions
using DataStructures
using LsqFit
using MacroTools
using Dates
using ProgressMeter
using Random
using JSON, GZip

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
    Model, @λ, select_maincomp!, SumReducer, domain, comptype,
    MultiModel, update!, isfreezed, thaw!, freeze!, fit

include("HashVector.jl")
include("domain.jl")

# ====================================================================
"""
    FunctDesc

A "Julia function" descriptor containing the reference to the function itself, a string representation of its source code definition (for displaying purposes) and the lists of its arguments.  It can be invoked using the standard syntax for functions

### Example:
```
julia> f = GFit.FunctDesc((x, p=0) -> x + p,    # actual function definition
                          "(x, p=0) -> x + p",  # string representation
                          [:x],                 # vector of argument namess
                          [:(p = 0)])           # vector of `Expr` with arguments default values
julia> f(1, 2)
3
```

Note that it is unpractical to directly create a `FunctDescr` using its constructor, and the above results can be obtained using the @λ macro:
```
f = @λ (x, p=0) -> x + p
```

"""
struct FunctDesc
    funct::Function
    display::String
    args::Vector{Symbol}   # positional arguments
    optargs::Vector{Expr}  # optional arguments with default values
end
(f::FunctDesc)(args...; kws...) = f.funct(args...; kws...)

"""
    @λ expr

Macro to generate a `FunctDesc` object using the same syntax as in a standard Julia anonymous function.

### Example
```
julia> f = @λ (x, p=0) -> x + p

julia> f(1, 2)
3
```
"""
macro λ(_expr)
    @assert isexpr(longdef(_expr), :function)
    expr = prettify(_expr)
    def  = splitdef(expr)
    args    = convert(Vector{Symbol}, filter(x -> isa(x, Symbol), def[:args]))
    optargs = convert(Vector{Expr}  , filter(x -> isa(x, Expr)  , def[:args]))
    return esc(:(GFit.FunctDesc($expr, string($(QuoteNode(expr))), $args, $optargs)))
end


# ====================================================================
"""
    Parameter

A structure representing a model parameter.

# Fields:
 - `val::Float64`: parameter value (initial guess before fitting, or best fit one after fitting);
 - `low::Float64`: lower limit for the value (default: `-Inf`);
 - `high::Float64`: upper limit for the value (default: `+Inf`);
 - `fixed::Bool`: whether the parameter is fixed during fitting (default: `false`);
 - `patch::Union{Nothing, Symbol, FunctDesc}`: patch prescription within the same model;
 - `mpatch::Union{Nothing, FunctDesc}`: patch prescription in a multi-model analysis;
 - `actual::Float64`: actual value for the parameter (i.e. after applying the patch prescription)`;
 - `unc::Float64`: 1σ uncertainty associated to the parameter value.

Note: the `Parameter` fields are supposed to be accessed directly by the user, without invoking any get/set method.
"""
mutable struct Parameter
    val::Float64
    low::Float64              # lower limit value
    high::Float64             # upper limit value
    fixed::Bool
    patch::Union{Nothing, Symbol, FunctDesc}
    mpatch::Union{Nothing, FunctDesc}
    actual::Float64
    unc::Float64
end
Parameter(value::Number) = Parameter(float(value), -Inf, +Inf, false, nothing, nothing, NaN, NaN)


# ====================================================================
# Components:
#
# A *component* is a generic implementation of a building block for a
# model. It must inherit `AbstractComponent` and implement the
# `evaluate!` method.  The structure should contain zero or more field
# of type Parameter, or have all parameters collected in a single
# field of type OrderedDict{Symbol, Parameter}()
abstract type AbstractComponent end

# Note: this function must mirror setparams!()
function getparams(comp::AbstractComponent)
    out = OrderedDict{Symbol, Parameter}()
    for name in fieldnames(typeof(comp))
        field = getfield(comp, name)
        if isa(field, Parameter)
            out[name] = field
        elseif isa(field, OrderedDict{Symbol, Parameter})
            @assert length(out) == 0  # avoid parameter name clash
            return field
        end
    end
    return out
end

# Note: this function must mirror getparams()
function setparams!(comp::AbstractComponent, params::HashVector{Parameter})
    for name in fieldnames(typeof(comp))
        field = getfield(comp, name)
        if isa(field, Parameter)
            field.val    = params[name].val
            field.unc    = params[name].unc
            field.actual = params[name].actual
        elseif isa(field, OrderedDict{Symbol, Parameter})
            for (name, par) in field
                par.val    = params[name].val
                par.unc    = params[name].unc
                par.actual = params[name].actual
            end
        end
    end
    nothing
end


# Fall back methods
dependencies(comp::AbstractComponent) = Symbol[]

prepare!(comp::AbstractComponent, domain::AbstractDomain) =
    fill(NaN, length(domain))

evaluate!(buffer::Vector{Float64}, comp::T, domain::AbstractDomain, pars...) where T <: AbstractComponent=
    error("No evaluate!() method implemented for $T")

# Built-in components
include("components/FComp.jl")
include("components/OffsetSlope.jl")
include("components/Polynomial.jl")
include("components/Gaussian.jl")
include("components/Lorentzian.jl")
include("components/SumReducer.jl")


# ====================================================================
# CompEval: a wrapper for a component evaluated on a specific domain
#
mutable struct CompEval{TComp <: AbstractComponent, TDomain <: AbstractDomain}
    comp::TComp
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
            comp, 0,
            Vector{Vector{Float64}}(),
            fill(NaN, length(getparams(comp))),
            buffer, false, false)
    end
end


function update!!(c::CompEval, domain::AbstractDomain, pvalues::Vector{Float64})
    c.updated  &&  return

    # Do we actually need a new evaluation?
    if (any(c.lastvalues .!= pvalues)  ||  (c.counter == 0)  ||  (length(c.deps) > 0))
        if length(c.deps) > 0
            evaluate!(c.buffer, c.comp, domain, c.deps, pvalues...)
        else
            evaluate!(c.buffer, c.comp, domain, pvalues...)
        end
        c.lastvalues .= pvalues
        c.counter += 1
    end
    c.updated = true
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
- a single `FunctDesc`, which will be wrapped into an `LComp` component and a default name will be assigned (`:main`);
- one or more `Pair{Symbol, AbstractComponent}`, where the first element is the name and the second is the component.

You may access the individual component in a `Model` using the indexing syntax, as if it was a `Dict{Symbol, AbstractComponent}`.  Also, you may add new components to a `Model` after it has been created using the same synatx.  Finally, you may use the `keys()` and `haskey()` functions with their usual meanings.

Individual components may be *freezed* (i.e. have all its parameters fixed during fitting, despite the individual `Parameter` settings) or *thawed* using the `freeze!()` and `thaw!()` functions.  Use the `isfreezed()` function to check if a component is freezed.

The main component, i.e. the one whose evaluation corresponds to the overall model evaluation, is typically automatically identified by analyzing the component dependencies.  However a specific component may be forced to be the main one by invoking `select_maincomp!`.

The model is automatically evaluated whenever needed, however there are a few cases where it is not possible to trigger an automatic evaluation, e.g. immediately after the user modifies a `Parameter` value. In this case an evaluation can be forced by invoking `update!()`.

The most important function for a `Model` object is `fit()`, which allows to fit the model against an empirical dataset. The `!` in the name reminds us that, after fitting, the parameter values will be set to the best fit ones (rather than retaining their original values).

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
                # isa(item, Number)  &&  (item = SimplePar(item))
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
                elseif isa(arg[2], FunctDesc)
                    out[arg[1]] = FComp(arg[2])
                # elseif isa(arg[2], Number)
                #     out[arg[1]] = SimplePar(arg[2])
                else
                    error("Unsupported data type: " * string(typeof(arg[2])) *
                          ".  Must be an AbstractComponent, a FunctDesc or a real number.")
                end
            end
            return parse_args(out)
        end

        parse_args(arg::AbstractComponent) = parse_args(:main => arg)
        parse_args(arg::FunctDesc) = parse_args(:main => FComp(arg))
        # parse_args(arg::Real) = parse_args(:main => SimplePar(arg))

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
        for d in dependencies(model, cname)
            i = findfirst(maincomps .== d)
            deleteat!(maincomps, i)
        end
    end

    if length(maincomps) > 1
        # Ignoring components with no dependencies
        for (cname, ceval) in model.cevals
            if length(dependencies(model, cname)) == 0
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


function dependencies(model::Model, cname::Symbol; select_domain=false)
    domdeps = Vector{Symbol}()
    compdeps = Vector{Symbol}()
    nd = ndims(domain(model))
    for d in dependencies(model.cevals[cname].comp)
        if haskey(model.cevals, d)
            # Dependency with known name
            push!(compdeps, d)
        else
            # Dependency with unknown name is intended as a domain dimension
            @assert length(compdeps) == 0 "Domain dependencies must be listed first"
            @assert length(domdeps) < nd "Component $cname depends on $d, but the latter is not a component in the model."
            push!(domdeps, d)
        end
    end
    @assert (length(domdeps) == 0)  ||  (length(domdeps) == nd) "Domain has $nd dimensions but only $(length(domdeps)) are listed as dependencies"
    return (select_domain ? domdeps : compdeps)
end



"""
    update!(model::Model)

Evaluate a `Model` and update internal structures.
"""
function update!(model::Model)
    update_step0(model)
    # update_step1()
    update_step2(model)
    update_step3(model)
    update_step4(model)
    return model
end


# Evaluation step 0: update internal structures before fitting
function update_step0(model::Model)
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
            model.actual[ cname][pname] = par.val

            if !isnothing(par.patch)
                @assert isnothing(par.mpatch) "Parameter [$cname].$pname has both patch and mpatch fields set, while only one is allowed"
                if isa(par.patch, Symbol)  # use same param. value from a different component
                    par.fixed = true
                else                       # invoke a patch function
                    @assert length(par.patch.args) in [1,2]
                    if length(par.patch.args) == 1
                        par.fixed = true
                    else
                        par.fixed = false
                    end
                end
            elseif !isnothing(par.mpatch)
                @assert !isnothing(model.parent) "Parameter [$cname].$pname has the mpatch field set but no MultiModel has been created"
                @assert length(par.mpatch.args) in [1,2]
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
        i = 1
        for d in dependencies(model, cname, select_domain=true)
            push!(ceval.deps, coords(domain(model), i))
        end
        for d in dependencies(model, cname, select_domain=false)
            push!(ceval.deps, model.buffers[d])
        end
    end
end


# Evaluation step 1: set new model parameters
function update_step1(model::Model, pvalues::Vector{Float64})
    internal_data(model.pvalues)[model.ifree] .= pvalues
end


# Evaluation step 2: copy all parameter values into actual, then
# update the latter by invoking the user patch functions.
function update_step2(model::Model)
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
                        model.actual[cname][pname] = par.patch(model.pvalues, model.pvalues[cname][pname])
                    end
                end
            elseif !isnothing(par.mpatch)
                @assert !isnothing(model.parent) "Parameter [$cname].$pname has the mpatch field set but no MultiModel has been created"
                if length(par.mpatch.args) == 1
                    model.actual[cname][pname] = par.mpatch(model.parent.pvalues)
                else
                    model.actual[cname][pname] = par.mpatch(model.parent.pvalues, model.pvalues[cname][pname])
                end
            end
        end
    end
end


# Evaluation step 3: actual evaluation of model components, starting
# from the main one and following dependencies
update_step3(model::Model) = update_step3(model, find_maincomp(model))
function update_step3(model::Model, cname::Symbol)
    # Recursive evaluation of dependencies
    for d in dependencies(model, cname)
        update_step3(model, d)
    end
    update!!(model.cevals[cname], model.domain, values(model.actual[cname]))
end


# Evaluation step 4: copy back bestfit and actual values, as well as
# uncertainties, into their original Parameter structures.
function update_step4(model::Model, uncerts=Vector{Float64}[])
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
        setparams!(ceval.comp, model.params[cname])
    end
end


# User interface
# setindex!(model::Model, v::Real, cname::Symbol) = setindex!(model, SimplePar(v), cname)
setindex!(model::Model, f::FunctDesc, cname::Symbol) = setindex!(model, FComp(f), cname)
function setindex!(model::Model, comp::AbstractComponent, cname::Symbol)
    ceval = CompEval(comp, model.domain)
    model.cevals[cname] = ceval
    model.buffers[cname] = ceval.buffer
    update!(model)
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
    update!(model)
    nothing
end

"""
    thaw!(model::Model, cname::Symbol)

Thaw a freezed component in the model (i.e. treat component parameters as fixed only if explicitly set in the corresponding `Parameter` structure).
"""
function thaw!(model::Model, cname::Symbol)
    @assert cname in keys(model.cevals) "Component $cname is not defined"
    model.cevals[cname].cfixed = false
    update!(model)
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


"""
    comptype(model::Model, cname::Symbol)

Return a component type as a string
"""
comptype(model::Model, cname::Symbol) = string(typeof(model[cname]))


"""
    comptypes(model::Model)

Return a `OrderedDict{Symbol, String}` with the model component types.
"""
comptypes(model::Model) = OrderedDict([cname => comptype(model, cname) for cname in keys(model)])


# Return model evaluations
(model::Model)() = reshape(domain(model), model.cevals[find_maincomp(model)].buffer)
(model::Model)(name::Symbol) = reshape(domain(model), model.cevals[name].buffer)

"""
    select_maincomp!(model::Model, cname::Symbol)

Force a component to be the final one for model evaluation.
"""
function select_maincomp!(model::Model, cname::Symbol)
    @assert haskey(model, cname)
    model.maincomp = cname
end


struct ModelSnapshot
    domain::AbstractDomain
    params::HashHashVector{Parameter}
    buffers::OrderedDict{Symbol, Vector{Float64}}
    maincomp::Symbol
    comptypes::OrderedDict{Symbol, String}
    show::String
end
function ModelSnapshot(model::Model)
    io = IOBuffer()
    if showsettings.plain
        show(io , model)
    else
        ctx = IOContext(io, :color => true)
        show(ctx, model)
    end
    s = String(take!(io))
    ModelSnapshot(deepcopy(domain(model)), deepcopy(model.params),
                  deepcopy(model.buffers), find_maincomp(model),
                  comptypes(model), s)
end

domain(model::ModelSnapshot) = model.domain
Base.keys(model::ModelSnapshot) = collect(keys(model.buffers))
(model::ModelSnapshot)() = reshape(domain(model), model.buffers[model.maincomp])
(model::ModelSnapshot)(name::Symbol) = reshape(domain(model), model.buffers[name])
comptype(model::ModelSnapshot, cname::Symbol) = model.comptypes[cname]
Base.haskey(m::ModelSnapshot, name::Symbol) = haskey(m.params, name)
function Base.getindex(model::ModelSnapshot, name::Symbol)
    if name in keys(model.params)
        return model.params[name]
    end
    error("Name $name not defined")
end


include("multimodel.jl")

abstract type AbstractFitProblem end
include("minimizers.jl")
include("fit.jl")
include("serialize.jl")

include("show.jl")
include("utils.jl")

end
