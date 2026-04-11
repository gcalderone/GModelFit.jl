module GModelFit

using Printf, PrettyTables
using Statistics
using DataStructures
using LsqFit
using MacroTools
using Dates
using Random
using StaticDictTrees
using DataFrames

import Base: show, ndims, size, length, iterate, keys, haskey, getindex, setindex!, values, empty!, view, propertynames, setproperty!

import ForwardDiff: Dual
import StaticDictTrees: get_layer

export AbstractDomain, Domain, CartesianDomain, Measures, GaussianData, PoissonCounts,
    coords, gridcoords, uncerts,
    @fd, SumReducer, domain, stype, dependencies, getcomp, getresp,
    getgroups, getcomps, getparams, getfreepars, flattenkeys,
    Model, ModelEval, Likelihood, snapshot,
    isfrozen, thaw!, freeze!, fit, fit!, loglh, gofstat, dof

include("domain.jl")

# --------------------------------------------------------------------
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

# --------------------------------------------------------------------
# Abstract types

abstract type AbstractComponent end
abstract type AbstractModel{     GT <: Tuple, CT <: Tuple, PT <: Tuple} <: AbstractDict{Any, Any} end
abstract type AbstractModelEval{ GT <: Tuple, CT <: Tuple, PT <: Tuple} <: AbstractModel{GT, CT, PT} end
abstract type AbstractLikelihood{GT <: Tuple, CT <: Tuple, PT <: Tuple} <: AbstractModelEval{GT, CT, PT} end

# --------------------------------------------------------------------
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
abstract type AbstractParameter end

macro _Parameter_fields()
    esc(quote
            val::Float64
            low::Float64              # lower limit value
            high::Float64             # upper limit value
            fixed::Bool
            patch::Union{Nothing, Symbol, FunctDesc}
            reparam::Union{Nothing, FunctDesc}
        end)
end

mutable struct Parameter <: AbstractParameter
    @_Parameter_fields()
end

mutable struct ParameterEval <: AbstractParameter
    @_Parameter_fields()
    _dirty::Bool
    actually_fixed::Bool
    unc::Float64
    mval::Float64
end

# --------------------------------------------------------------------
abstract type AbstractCompSlot     <: AbstractDict{Symbol, Parameter} end
abstract type AbstractCompSlotEval <: AbstractCompSlot end

macro _CompSlot_fields()
    esc(quote
            comp::TComp
            frozen::Bool
        end)
end

mutable struct CompSlot{TComp <: AbstractComponent} <: AbstractCompSlot
    rootmodel::AbstractModel
    key::Tuple
    @_CompSlot_fields()
    CompSlot(rootmodel::AbstractModel, key::Tuple, comp::TComp) where {TComp <: AbstractComponent} =
        new{TComp}(rootmodel, key, deepcopy(comp), false)
end

mutable struct CompSlotEval{T <: Real, TDomain <: AbstractDomain, TComp <: AbstractComponent} <: AbstractCompSlotEval
    rootmodel::AbstractModelEval
    key::Tuple
    @_CompSlot_fields()
    pactual::AbstractVector{T}
    _dirty::Bool
    domain::TDomain
    counter::Int
    lastparvalues::Vector{T}
    lastdepscounter::Vector{Int}
    deps::Vector{CompSlotEval}
    evalbuf::Array{T}

    function CompSlotEval(rootmodel::AbstractModelEval, key::Tuple, comp::CompSlot{TComp}, domain::TDomain) where {TComp <: AbstractComponent, TDomain <: AbstractDomain}
        T = NumType(rootmodel)
        prepare!(comp.comp, domain)
        return new{T, TDomain, TComp}(rootmodel, key, deepcopy(comp.comp), comp.frozen, view(Float64[], Int64[]),
                                      true, domain, 0,
                                      Vector{T}(undef, length(getparams(comp.comp))),
                                      Vector{Int}(),
                                      Vector{CompSlotEval}(),
                                      domain2buffer(T, domain))
    end
    function CompSlotEval(rootmodel::AbstractModelEval, key::Tuple, comp::TComp, domain::TDomain) where {TComp <: AbstractComponent, TDomain <: AbstractDomain}
        T = NumType(rootmodel)
        prepare!(comp, domain)
        return new{T, TDomain, TComp}(rootmodel, key, deepcopy(comp), false, true, domain, 0,
                                      Vector{T}(undef, length(getparams(comp))),
                                      Vector{Int}(),
                                      Vector{CompSlotEval}(),
                                      domain2buffer(T, domain))
    end
end

# --------------------------------------------------------------------
abstract type AbstractGroup <: AbstractDict{Any, Any} end
abstract type AbstractGroupEval <: AbstractGroup end
abstract type AbstractGroupLH <: AbstractGroupEval end

struct Group <: AbstractGroup
    rootmodel::AbstractModel
    key::Tuple
    Group(rootmodel::AbstractModel, key::Tuple) = new(rootmodel, key)
end

macro _GroupEval_fields()
    esc(quote
            resp::AbstractResponse
            patched::OrderedDict{NTuple{2, Symbol}, ParameterEval}
            seq::Vector{Symbol}
            evalbuf::Array{T}
        end)
end
struct GroupEval{T} <: AbstractGroupEval
    rootmodel::AbstractModelEval
    key::Tuple
    @_GroupEval_fields()

    function GroupEval(rootmodel::AbstractModelEval, key::Tuple, resp::AbstractResponse)
        T = NumType(rootmodel)
        return new{T}(rootmodel, key, resp,
                      OrderedDict{NTuple{2, Symbol}, ParameterEval}(),
                      Vector{Symbol}(),
                      domain2buffer(T, folded_domain(resp)))
    end
end

struct GroupLH{T <: Real, DT <: AbstractData} <: AbstractGroupLH
    rootmodel::AbstractModelEval
    key::Tuple
    @_GroupEval_fields()
    data::DT
    residuals::Array{T}

    function GroupLH(rootmodel::AbstractModelEval, key::Tuple, data::DT) where {DT <: AbstractData}
        T = NumType(rootmodel)
        resp = getresp(data)
        return new{T, DT}(rootmodel, key, resp,
                          OrderedDict{NTuple{2, Symbol}, Parameter}(),
                          Vector{Symbol}(),
                          domain2buffer(T, folded_domain(resp)),
                          data, Vector{T}(undef, length(folded_domain(resp))))
    end
end

# --------------------------------------------------------------------
"""
    Model

A structure containing a model description.

Constructor is: `Model(components...)`.  Components may be specified as:
- a single `Dict{Symbol, AbstractComponent}`, where the keys are the names and the values the component objects;
- a single component (the default `:main` is automatically assigned);
- a single `FunctDesc` which will be wrapped into an `FComp` component and a default name will be assigned (`:main`);
- one or more `Pair{Symbol, AbstractComponent}`, where the first element is the name and the second is the component.

You may access the individual component in a `Model` using the indexing syntax, as if it was a `Dict{Symbol, AbstractComponent}`.  Also, you may add new components to a `Model` after it has been created using the same syntax.  Finally, you may use the `keys()` function to retrieve the list of component names.

Individual components may be *frozen* (i.e. have all its parameters fixed during fitting, despite the individual `Parameter` settings) or *thawed* using the `freeze!()` and `thaw!()` functions.  Use the `isfrozen()` function to check if a component is frozen.

The main component, i.e. the one whose evaluation corresponds to the overall model evaluation, is automatically identified by analyzing the component dependencies.

The most important function for a `Model` object is `fit()` which allows to fit the model against an empirical dataset.  The `fit!()` function has the same purpose, with the only difference that it stores the best fit parameter values into the original `Model` object.

The model and all component evaluation can be evaluated has if they were a function by simply passing a `Domain` object.
"""

macro _Model_fields()
    esc(quote;
            dt::DictTree
        end)
end

struct Model{GT, CT, PT} <: AbstractModel{GT, CT, PT}
    @_Model_fields()

    Model() = Model{Tuple{}}()
    function Model{GT}() where {GT <: Tuple}
        CT = Tuple{fieldtypes(GT)..., Symbol}
        PT = Tuple{fieldtypes(GT)..., Symbol, Symbol}

        dt = DictTree()
        add_layer!(dt, SDTree{GT, Group}(    on_delete=forbid_delete) ,label=:groups)
        add_layer!(dt, SDTree{CT, CompSlot}( on_delete=forbid_delete), label=:comps)
        add_layer!(dt, SDTree{PT, Parameter}(on_delete=forbid_delete), label=:params)
        return new{GT, CT, PT}(dt)
    end
end

macro _ModelEval_fields()
    esc(quote;
            @_Model_fields()
            pvalues::SDTree{PT, T}
            pactual::SDTree{PT, T}
            ifree::Vector{Int}
        end)
end

struct ModelEval{T <: Real, GT, CT, PT} <: AbstractModelEval{GT, CT, PT}
    @_ModelEval_fields()

    function ModelEval(resp::AbstractDict{GT, <: AbstractResponse}; autodiff=false) where {GT <: Tuple}
        CT = Tuple{fieldtypes(GT)..., Symbol}
        PT = Tuple{fieldtypes(GT)..., Symbol, Symbol}

        dt = DictTree()
        add_layer!(dt, SDTree{GT, GroupEval}(    on_delete=forbid_delete) ,label=:groups)
        add_layer!(dt, SDTree{CT, CompSlotEval}( on_delete=forbid_delete), label=:comps)
        add_layer!(dt, SDTree{PT, ParameterEval}(on_delete=forbid_delete), label=:params)
        T = NumType(autodiff)
        model = new{T, GT, CT, PT}(dt, SDTree{PT, Float64}(), SDTree{PT, Float64}(), Vector{Int}())
        for (gkey, resp) in resp
            model.dt[gkey] = GroupEval(model, gkey, resp)
        end
        return model
    end

    function ModelEval(resp::AbstractDict{GT, <: AbstractResponse}, model::Model{GT}; autodiff=false) where {GT <: Tuple}
        (keys(resp) == keys(getgroups(model)))  ||  error("Keys in input dict and model are different")
        out = ModelEval(resp; autodiff=autodiff)
        for (ckey, comp) in getcomps(model)
            gkey = group_key(model, ckey)
            out.dt[ckey] = CompSlotEval(out, ckey, comp, unfolded_domain(resp[gkey]))
        end
        for (pkey, par) in getparams(model)
            out.dt[pkey] = ParameterEval(par)
        end
        evaluate!(out)
        return out
    end
end

struct Likelihood{T <: Real, GT, CT, PT} <: AbstractLikelihood{GT, CT, PT}
    @_ModelEval_fields()
    evalbuf::Vector{T}  # local buffer containing residuals from all groups

    function Likelihood(data::AbstractDict{GT, <: AbstractData}; autodiff=false) where {GT <: Tuple}
        CT = Tuple{fieldtypes(GT)..., Symbol}
        PT = Tuple{fieldtypes(GT)..., Symbol, Symbol}

        dt = DictTree()
        add_layer!(dt, SDTree{GT, GroupLH}(      on_delete=forbid_delete) ,label=:groups)
        add_layer!(dt, SDTree{CT, CompSlotEval}( on_delete=forbid_delete), label=:comps)
        add_layer!(dt, SDTree{PT, ParameterEval}(on_delete=forbid_delete), label=:params)
        T = NumType(autodiff)
        model = new{T, GT, CT, PT}(dt, SDTree{PT, T}(), SDTree{PT, T}(), Vector{Int}(), Vector{T}(undef, sum(length.(collect(values(data))))))
        for (gkey, data) in data
            model.dt[gkey] = GroupLH(model, gkey, data)
        end
        evaluate!(model)
        return model
    end

    function Likelihood(data::AbstractDict{GT, <: AbstractData}, model::Model{GT}; autodiff=false) where {GT <: Tuple}
        (keys(data) == keys(getgroups(model)))  ||  error("Keys in input dict and model are different")
        out = Likelihood(data; autodiff=autodiff)
        T = NumType(autodiff)
        for (ckey, comp) in getcomps(model)
            gkey = group_key(model, ckey)
            out.dt[ckey] = CompSlotEval(out, ckey, comp, unfolded_domain(getresp(data[gkey])))
        end
        for (pkey, par) in getparams(model)
            out.dt[pkey] = ParameterEval(par)
        end
        evaluate!(out)
        return out
    end
end

# --------------------------------------------------------------------
include("components/FComp.jl")
include("components/OffsetSlope.jl")
include("components/Polynomial.jl")
include("components/Gaussian.jl")
include("components/Lorentzian.jl")
include("components/SumReducer.jl")

include("snapshot.jl")
include("dependencies.jl")
include("utils.jl")
include("methods.jl")
include("show.jl")
include("serialize.jl")
include("gnuplot_recipe.jl")
include("export.jl")
include("precompile.jl")

end
