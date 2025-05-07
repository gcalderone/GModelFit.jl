module GModelFit

using Printf, PrettyTables
using Statistics
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
import Base.iterate
import Base.values
import Base.push!
import Base.empty!

export AbstractDomain, Domain, CartesianDomain, coords, axis, Measures, uncerts,
    CompEval, Model, @fd, select_maincomp!, SumReducer, domain, comptype,
    isfreezed, thaw!, freeze!, fit, fit!, compare

include("PV.jl")
using .PV

include("domain.jl")

# ====================================================================
"""
    FunctDesc

A "Julia function" descriptor containing the reference to the function itself, a string representation of its source code definition (for displaying purposes) and the lists of its arguments.  It can be invoked using the standard syntax for functions

### Example:
```
julia> f = GModelFit.FunctDesc( (x, p=0) -> x + p,   # actual function definition
                               "(x, p=0) -> x + p",  # string representation
                               [:x],                 # vector of argument namess
                               [:(p = 0)])           # vector of `Expr` with arguments default values
julia> f(1, 2)
3
```

Note that it is inconvenient to directly create a `FunctDescr` using its constructor, and the above results can be obtained by using the @fd macro:
```
f = @fd (x, p=0) -> x + p
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
    @fd expr

Macro to generate a `FunctDesc` object using the same syntax as for a standard Julia anonymous function.

### Example
```
julia> f = @fd (x, p=0) -> x + p

julia> f(1, 2)
3
```
"""
macro fd(_expr)
    @assert isexpr(longdef(_expr), :function)
    expr = prettify(_expr)
    def  = splitdef(expr)
    args    = convert(Vector{Symbol}, filter(x -> isa(x, Symbol), def[:args]))
    optargs = convert(Vector{Expr}  , filter(x -> isa(x, Expr)  , def[:args]))
    return esc(:(GModelFit.FunctDesc($expr, string($(QuoteNode(expr))), $args, $optargs)))
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
abstract type AbstractCompWDeps end

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
function setparams!(comp::AbstractComponent, params::PVComp{Parameter})
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


"""
    prepare!(comp::AbstractComponent, domain::AbstractDomain)

Allocate the buffer for a component evaluation on a specific domain.  Return value must be a `Vector{Float64}`.

This function is invoked only once when the `ModelEval` structure is created (typically within a `fit` of `fit!` call), hence it is the perfect place to pre-compute quantities associated to a component evaluation on a specific domain.  Default implementation returns a vector filled with `NaN`s with the same length as the domain.
"""
prepare!(comp::AbstractComponent, domain::AbstractDomain) = fill(NaN, length(domain))

"""
    dependencies(comp::AbstractComponent)

Return the name of dependecies for a component. Return value must be a `Vector{Symbol}`.

Default implementation returns `Symbol[]` (i.e. no dependencies.
"""
dependencies(comp::AbstractComponent) = Symbol[]


# ====================================================================
# Model
#

"""
    Model

A structure containing a model description.

Constructor is: `Model(components...)`.  Components may be specified as:
- a single `Dict{Symbol, AbstractComponent}`, where the keys are the names and the values the component objects;
- a single component (the default `:main` is automatically assigned);
- a single `FunctDesc` which will be wrapped into an `FComp` component and a default name will be assigned (`:main`);
- one or more `Pair{Symbol, AbstractComponent}`, where the first element is the name and the second is the component.

You may access the individual component in a `Model` using the indexing syntax, as if it was a `Dict{Symbol, AbstractComponent}`.  Also, you may add new components to a `Model` after it has been created using the same syntax.  Finally, you may use the `keys()` and `haskey()` functions with their usual meanings.

Individual components may be *freezed* (i.e. have all its parameters fixed during fitting, despite the individual `Parameter` settings) or *thawed* using the `freeze!()` and `thaw!()` functions.  Use the `isfreezed()` function to check if a component is freezed.

The main component, i.e. the one whose evaluation corresponds to the overall model evaluation, is automatically identified by analyzing the component dependencies.  However a specific component may be forced to be the main one by invoking `select_maincomp!`.

The most important function for a `Model` object is `fit()` which allows to fit the model against an empirical dataset.  The `fit!()` function has the same purpose, with the only difference that it stores the best fit parameter values into the original `Model` object.

The model and all component evaluation can be evaluated has if they were a function by simply passing a `Domain` object.
"""
mutable struct Model
    comps::OrderedDict{Symbol, AbstractComponent}
    fixed::OrderedDict{Symbol, Bool}
    maincomp::Union{Nothing, Symbol}

    Model() = new(OrderedDict{Symbol, AbstractComponent}(), OrderedDict{Symbol, Bool}(), nothing)

    function Model(dict::AbstractDict)
        model = Model()
        for (name, item) in dict
            # isa(item, Number)  &&  (item = SimplePar(item))
            @assert isa(name, Symbol)
            @assert isa(item, AbstractComponent)
            model[name] = item
        end
        return model
    end

    function Model(args::Vararg{Pair})
        model = Model()
        for arg in args
            @assert isa(arg[1], Symbol)
            if isa(arg[2], AbstractComponent)
                model[arg[1]] = arg[2]
            elseif isa(arg[2], FunctDesc)
                model[arg[1]] = FComp(arg[2])
                # elseif isa(arg[2], Number)
                #     out[arg[1]] = SimplePar(arg[2])
            else
                error("Unsupported data type: " * string(typeof(arg[2])) *
                    ".  (accepted types ar T <: AbstractComponent or FunctDesc.")
            end
        end
        return model
    end
    Model(arg::AbstractComponent) = Model(:main => arg)
    Model(arg::FunctDesc) = Model(:main => FComp(arg))
end


function find_maincomp(model::Model)
    # Identify parent for all comps
    parent = OrderedDict{Symbol, Symbol}()
    for (cname, comp) in model.comps
        for d in GModelFit.dependencies(model, cname)
            @assert !haskey(parent, d) "Component $d has two parent nodes: $(parent[d]) and $cname"
            parent[d] = cname
        end
    end

    # Ensure no circular dependency is present by checking all parent
    # nodes of a given component to be different from the component
    # itself.  Also collect components with no parent.
    comps_with_no_parent = Vector{Symbol}()
    for cname in keys(model.comps)
        if haskey(parent, cname)
            p = parent[cname]
            @assert p != cname "Component $cname depends on itself"
            while haskey(parent, p)
                p = parent[p]
                if cname == p
                    display(parent)
                    error("Circular dependency detected for component $cname")
                end
            end
        else
            push!(comps_with_no_parent, cname)
        end
    end

    # If multiple possibilities are stll available neglect components
    # with no dependencies
    while length(comps_with_no_parent) > 1
        if length(dependencies(model, comps_with_no_parent[1])) == 0
            deleteat!(comps_with_no_parent, 1)
        end
    end

    # The above check is always performed even if an explicit maincomp
    # has been set
    if isnothing(model.maincomp)
        return comps_with_no_parent[end]
    else
        return model.maincomp
    end
end


function dependencies(model::Model, cname::Symbol; select_domain=false)
    domdeps = Vector{Symbol}()
    compdeps = Vector{Symbol}()
    # nd = ndims(domain(model))
    for d in dependencies(model.comps[cname])
        if haskey(model.comps, d)  # dependency with known name
            push!(compdeps, d)
        else # dependency with unknown name is intended as a domain dimension
            @assert length(compdeps) == 0 "Domain dependencies must be listed first"
            # @assert length(domdeps) < nd "Component $cname depends on $d, but the latter is not a component in the model."
            push!(domdeps, d)
        end
    end
    # @assert (length(domdeps) == 0)  ||  (length(domdeps) == nd) "Domain has $nd dimensions but only $(length(domdeps)) are listed as dependencies"
    return (select_domain  ?  domdeps  :  compdeps)
end


# User interface
setindex!(model::Model, f::FunctDesc, cname::Symbol) = model[cname] = FComp(f)
function setindex!(model::Model, comp::AbstractComponent, cname::Symbol)
    model.comps[cname] = deepcopy(comp)
    model.fixed[cname] = false
end

function iterate(model::Model, i=1)
    k = collect(keys(model))
    (i > length(k))  &&  return nothing
    return (k[i] => model[k[i]], i+1)
end


"""
    isfreezed(model::Model, cname::Symbol)

Check whether a component is *freezed* in the model.
"""
function isfreezed(model::Model, cname::Symbol)
    @assert cname in keys(model.fixed) "Component $cname is not defined"
    return model.fixed[cname]
end

"""
    freeze!(model::Model, cname::Symbol)

Freeze a component in the model (i.e. treat all component parameters as fixed for fitting).
"""
function freeze!(model::Model, cname::Symbol)
    @assert cname in keys(model.fixed) "Component $cname is not defined"
    model.fixed[cname] = true
    nothing
end

"""
    thaw!(model::Model, cname::Symbol)

Thaw a freezed component in the model (i.e. treat component parameters as fixed only if explicitly set in the corresponding `Parameter` structure).
"""
function thaw!(model::Model, cname::Symbol)
    @assert cname in keys(model.fixed) "Component $cname is not defined"
    model.fixed[cname] = false
    nothing
end


Base.keys(model::Model) = collect(keys(model.comps))


"""
    haskey(m::Model, name::Symbol)

Check whether a component exists in model.
"""
Base.haskey(model::Model, cname::Symbol) = haskey(model.comps, cname)


"""
    getindex(model::Model, cname::Symbol)

Return the model component with name `cname`.
"""
function Base.getindex(model::Model, cname::Symbol)
    @assert cname in keys(model.comps) "Component $cname is not defined"
    return model.comps[cname]
end

"""
    length(model::Model)

Return number of components in a model.
"""
Base.length(model::Model) = length(model.comps)


"""
    comptype(model::Model, cname::Symbol)

Return a component type as a string.
"""
comptype(model::Model, cname::Symbol) = replace(string(typeof(model[cname])), "GModelFit." => "")


"""
    comptypes(model::Model)

Return a `OrderedDict{Symbol, String}` with the model component types.
"""
comptypes(model::Model) = OrderedDict([cname => comptype(model, cname) for cname in keys(model)])


"""
    select_maincomp!(model::Model, cname::Symbol)

Force a component to be the final one for model evaluation.
"""
function select_maincomp!(model::Model, cname::Symbol)
    @assert haskey(model, cname) "Component $cname is not defined"
    model.maincomp = cname
end


include("evaluation.jl")
include("snapshot.jl")

abstract type AbstractFitProblem end
include("minimizers.jl")
include("fit.jl")
# TODO include("multimodel.jl")
include("serialize.jl")
include("show.jl")
include("utils.jl")
include("gnuplot_recipe.jl")
# TODO include("precompile.jl")

end
