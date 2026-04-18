module GModelFit

using Printf, PrettyTables
using Statistics
using DataStructures
using LsqFit
using MacroTools
using Dates
using Random

import Base.show
import Base.ndims
import Base.size
import Base.axes
import Base.length
import Base.iterate
import Base.keys
import Base.haskey
import Base.getindex
import Base.setindex!
import Base.values

import ForwardDiff: Dual

export AbstractDomain, Domain, CartesianDomain, coords, Measures, uncerts,
    Model, ModelSet, @fd, SumReducer, domain, compnames,
    isfrozen, thaw!, freeze!, set_IR!, fit, fit!, fitstat, select_maincomp!


include("PV.jl")
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
 - `patch::Union{Nothing, Symbol, FunctDesc}`: patch prescription;
 - `reparam::Union{Nothing, FunctDesc}`: re-interpret prescription;
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
    reparam::Union{Nothing, FunctDesc}
    actual::Float64
    unc::Float64
    actually_fixed::Bool
end
Parameter(value::Number) = Parameter(float(value), -Inf, +Inf, false, nothing, nothing, NaN, NaN, false)


# ====================================================================
# Components:
#
# A *component* is a generic implementation of a building block for a
# model. It must inherit `AbstractComponent` and implement the
# `evaluate!` method.  The structure is expected to have `Parameter`
# among its fields, or to implement a new method for the `getparams()`
# function.
abstract type AbstractComponent end

function getparams(comp::T) where T <: AbstractComponent
    out = OrderedDict{Symbol, Parameter}()
    for name in fieldnames(T)
        if fieldtype(T, name) == Parameter
            out[name] = getfield(comp, name)
        end
    end

    dict_fields = Vector{Symbol}()
    for name in fieldnames(T)
        if fieldtype(T, name) <: AbstractDict{Symbol, Parameter}
            push!(dict_fields, name)
        end
    end

    (length(dict_fields) == 0)  &&  (return out)
    @assert (length(out) != 0)  ||  (length(dict_fields) != 1) "Ambiguous set of Parameters in structure $T"
    return getfield(comp, dict_fields[1])
end


"""
    prepare!(comp::AbstractComponent, domain::AbstractDomain)

Invoked to precompute component-specific quantities

This method is invoked only once when the component is first evaluated hence it is the perfect place to pre-compute quantities associated to a component evaluation on a specific domain.
"""
prepare!(comp::AbstractComponent, domain::AbstractDomain) = nothing


"""
    dependencies(comp::AbstractComponent)

Return the name of dependecies for a component. Return value must be a `Vector{Symbol}`.

Default implementation returns `Symbol[]` (i.e. no dependencies).
"""
dependencies(comp::AbstractComponent) = Symbol[]


include("instrument_response.jl")


# ====================================================================
"""
    Model

A structure containing a model description.

Constructor is: `Model(components...)`.  Components may be specified as:
- a single `Dict{Symbol, AbstractComponent}`, where the keys are the names and the values the component objects;
- a single component (the default `:main` is automatically assigned);
- a single `FunctDesc` which will be wrapped into an `FComp` component and a default name will be assigned (`:main`);
- one or more `Pair{Symbol, AbstractComponent}`, where the first element is the name and the second is the component.

You may access the individual component in a `Model` using the indexing syntax, as if it was a `Dict{Symbol, AbstractComponent}`.  Also, you may add new components to a `Model` after it has been created using the same syntax.  Finally, you may use the `compnames()` function to retrieve the list of component names.

Individual components may be *frozen* (i.e. have all its parameters fixed during fitting, despite the individual `Parameter` settings) or *thawed* using the `freeze!()` and `thaw!()` functions.  Use the `isfrozen()` function to check if a component is frozen.

The main component, i.e. the one whose evaluation corresponds to the overall model evaluation, is automatically identified by analyzing the component dependencies.  However a specific component may be forced to be the main one by invoking `select_maincomp!`.

The most important function for a `Model` object is `fit()` which allows to fit the model against an empirical dataset.  The `fit!()` function has the same purpose, with the only difference that it stores the best fit parameter values into the original `Model` object.

The model and all component evaluation can be evaluated has if they were a function by simply passing a `Domain` object.
"""
mutable struct Model
    comps::OrderedDict{Symbol, AbstractComponent}
    frozen::OrderedDict{Symbol, Bool}
    maincomp::Union{Nothing, Symbol}
    IR::AbstractInstrumentResponse

    Model() = new(OrderedDict{Symbol, AbstractComponent}(),
                  OrderedDict{Symbol, Bool}(), nothing, IdealInstrument())
end

function Model(dict::AbstractDict)
    model = Model()
    for (name, item) in dict
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
        else
            error("Unsupported data type: " * string(typeof(arg[2])) *
                ".  (accepted types ar T <: AbstractComponent or FunctDesc.")
        end
    end
    return model
end
Model(arg::AbstractComponent) = Model(:main => arg)
Model(arg::FunctDesc) = Model(:main => FComp(arg))

# User interface
setindex!(model::Model, f::FunctDesc, cname::Symbol) = model[cname] = FComp(f)
function setindex!(model::Model, comp::AbstractComponent, cname::Symbol)
    comp = deepcopy(comp)
    model.comps[cname] = comp
    model.frozen[cname] = false
end


"""
    comptype(model::Model, cname::Symbol)

Return the type of a component as a String.
"""
comptype(model::Model, cname::Symbol) = replace(string(typeof(model[cname])), "GModelFit." => "")

"""
    isfrozen(model::Model, cname::Symbol)

Check whether a component is *frozen* in the model.
"""
function isfrozen(model::Model, cname::Symbol)
    @assert cname in keys(model.frozen) "Component $cname is not defined"
    return model.frozen[cname]
end

"""
    freeze!(model::Model, cname::Symbol)

Freeze a component in the model (i.e. treat all component parameters as fixed for fitting).
"""
function freeze!(model::Model, cname::Symbol)
    @assert cname in keys(model.frozen) "Component $cname is not defined"
    model.frozen[cname] = true
    nothing
end

"""
    thaw!(model::Model, cname::Symbol)

Thaw a frozen component in the model (i.e. treat component parameters as fixed only if explicitly set in the corresponding `Parameter` structure).
"""
function thaw!(model::Model, cname::Symbol)
    @assert cname in keys(model.frozen) "Component $cname is not defined"
    model.frozen[cname] = false
    nothing
end

"""
    select_maincomp!(model::Model, cname::Symbol)

Force a component to be the final one for model evaluation.
"""
function select_maincomp!(model::Model, cname::Symbol)
    @assert haskey(model, cname) "Component $cname is not defined"
    model.maincomp = cname
end

"""
    set_IR!(model::Model, IR::AbstractInstrumentResponse)

Set a model to use the instrument response passed as argument.

If this method is not invoked the GModelFit.IdealInstrument response will be used.
"""
function set_IR!(model::Model, IR::AbstractInstrumentResponse)
    model.IR = IR
end

# ====================================================================
struct ModelSet
    dict::OrderedDict{Symbol, Model}
    ModelSet() = new(OrderedDict{Symbol, Model}())
end

function ModelSet(dict::AbstractDict{Symbol, Model})
    ms = ModelSet()
    for (name, item) in dict
        ms[name] = item
    end
    return ms
end

function ModelSet(args::Vararg{Pair{Symbol, Model}})
    ms = ModelSet()
    for arg in args
        ms[arg[1]] = arg[2]
    end
    return ms
end

function setindex!(ms::ModelSet, model::Model, key::Symbol)
    ms.dict[key] = model
end

select_free(d::AbstractDict{K, Parameter}) where K = filter(p -> !p.second.actually_fixed, d)

# ====================================================================
include("dependencies.jl")
include("evaluation.jl")
include("snapshot.jl")
include("interface.jl")
include("fit.jl")
include("serialize.jl")
include("show.jl")
include("utils.jl")
include("ppl.jl")
include("gnuplot_recipe.jl")
include("precompile.jl")

end
