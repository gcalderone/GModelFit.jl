module GFit

using Printf, PrettyTables
using Statistics, Distributions
using DataStructures
using LsqFit
using ExprTools

import Base.show
import Base.ndims
import Base.size
import Base.length
import Base.haskey
import Base.keys
import Base.getindex
import Base.reshape
import Base.propertynames
import Base.getproperty
import Base.setproperty!
import Base.iterate

⋄ = getfield


export Domain, CartesianDomain, Measures,
    Prediction, Reducer, @reducer, add!, domain,
    Model, patch!, evaluate, isfixed, thaw, freeze, fit!,
    savelog

const MDict = OrderedDict{Symbol, Any}

include("domain.jl")

# ====================================================================
# Parameter
#
mutable struct Parameter
    val::Float64
    low::Float64              # lower limit value
    high::Float64             # upper limit value
    step::Float64
    fixed::Bool
    Parameter(value::Number) = new(float(value), -Inf, +Inf, NaN, false)
end

# ====================================================================
# A *component* is a generic implementation of a constituent part of a
# model.
#
# A component must inherit `AbstractComponent`, and implement the
# `ceval_data` and `evaluate` methods
abstract type AbstractComponent end


# ====================================================================
struct ParamID
    name::Symbol
    index::Int
    ParamID(pname::Symbol) = new(pname, 0)                  # scalar param
    ParamID(pname::Symbol, index::Int) = new(pname, index)  # vector of params
end

struct CompID
    id::Int
    name::Symbol
    CompID(cname::Symbol) = new(1, cname)  # default on first prediction
    CompID(id::Int, cname::Symbol) = new(id, cname)
end

struct CompParamID
    comp::CompID
    param::ParamID
end



# ====================================================================
function getparams(comp::AbstractComponent)
    params = OrderedDict{ParamID, Parameter}()
    for pname in fieldnames(typeof(comp))
        par = getfield(comp, pname)
        if isa(par, Parameter)
            params[ParamID(pname)] = par
        elseif isa(par, Vector{Parameter})
            for i in 1:length(par)
                params[ParamID(pname, i)] = par[i]
            end
        end
    end
    params
end

# ====================================================================
# CompEval: a wrapper for a component evaluated on a specific domain
#
mutable struct CompEval{TComp <: AbstractComponent, TDomain <: AbstractDomain}
    comp::TComp
    domain::TDomain
    params::OrderedDict{ParamID, Parameter}
    cdata
    counter::Int
    lastvalues::Vector{Float64}
    buffer::Vector{Float64}
    cfixed::Int8

    function CompEval(comp::AbstractComponent, domain::AbstractDomain)
        params = getparams(comp)
        cdata  = compeval_cdata(comp, domain)
        buffer = compeval_array(comp, domain)
        return new{typeof(comp), typeof(domain)}(
            comp, domain, params, cdata, 0,
            fill(NaN, length(params)),
            buffer, false)
    end
end


evaluate_cached(c::CompEval) = evaluate_cached(c, [par.val for par in values(c.params)])
function evaluate_cached(c::CompEval, pvalues::Vector{Float64})
    @assert length(c.params) == length(pvalues)

    # Do we actually need a new evaluation?
    if (any(c.lastvalues .!= pvalues)  ||  (c.counter == 0))
        c.lastvalues .= pvalues
        c.counter += 1
        @assert all(.!isnan.(pvalues))
        evaluate(c, pvalues...)
    end
    return c.buffer
end


# ====================================================================
# Component fall back methods
function compeval_cdata(comp::AbstractComponent, domain::AbstractDomain) end
#    error("Component " * string(typeof(comp)) * " must implement its own method for `compeval_cdata`.")

function compeval_array(comp::AbstractComponent, domain::AbstractDomain) end
#    error("Component " * string(typeof(comp)) * " must implement its own method for `compeval_array`.")

function evaluate(c::CompEval{TComp, TDomain}, args...) where {TComp, TDomain} end
#    error("Component " * string(TComp) * " must implement its own method for `evaluate`.")

function evaluate(c::CompEval{TComp, TDomain}) where {TComp, TDomain} end
#    error("Component " * string(TComp) * " must implement its own method for `evaluate`.")


# ====================================================================
# Built-in components
#
include("components/CDomain.jl")
include("components/SimplePar.jl")
include("components/FuncWrap.jl")
include("components/OffsetSlope.jl")
include("components/Gaussian.jl")


# ====================================================================
# Parse a dictionary or a collection of `Pair`s to extract all
# components
function extract_components(comp_iterable::AbstractDict)
    out = OrderedDict{Symbol, AbstractComponent}()
    for (name, comp) in comp_iterable
        isa(comp, Number)  &&  (comp = SimplePar(comp))
        @assert isa(name, Symbol)
        @assert isa(comp, AbstractComponent)
        out[name] = comp
    end
    return out
end

function extract_components(comp_iterable::Vararg{Pair})
    out = OrderedDict{Symbol, AbstractComponent}()
    for thing in comp_iterable
        name = thing[1]
        comp = thing[2]
        isa(comp, Number)  &&  (comp = SimplePar(comp))
        @assert isa(name, Symbol)
        @assert isa(comp, AbstractComponent)
        out[name] = comp
    end
    return out
end



# ====================================================================
mutable struct Reducer
    funct::Function
    allargs::Bool
    args::Vector{Symbol}
end

Reducer(f::Function) = Reducer(f, true, Vector{Symbol}())
Reducer(f::Function, args::AbstractVector{Symbol}) =  Reducer(f, false, collect(args))

macro reducer(ex)
    @assert ex.head == Symbol("->") "Not an anonmymous function"
    f = splitdef(ex; throw=false)
    @assert !isnothing(f) "Not an function definition"

    allargs = false
    args = Vector{Symbol}()
    if haskey(f, :args)
        if (length(f[:args]) == 1)  &&
            isa(f[:args][1], Expr)  &&  (f[:args][1].head == :...)
            allargs= true
        end
        if !allargs
            for arg in f[:args]
                if isa(arg, Symbol)
                    push!(args, arg)
                elseif isa(arg, Expr)  &&  (arg.head = Symbol("::"))
                    push!(args, arg.args[1])
                else
                    error("Unexpected input: $arg")
                end
            end
        end
    end
    # Here `esc` is necessary to evaluate the function expression in
    # the caller scope
    return esc(:(Reducer($ex, $allargs, $args)))
end


# ====================================================================
mutable struct ReducerEval
    args::Vector{Vector{Float64}}
    funct::Function
    counter::Int
    buffer::Vector{Float64}
end


# ====================================================================
# A model prediction suitable to be compared to experimental data
mutable struct Prediction
    orig_domain::AbstractDomain
    domain::AbstractLinearDomain
    cevals::OrderedDict{Symbol, CompEval}
    revals::OrderedDict{Symbol, ReducerEval}
    rsel::Symbol
    counter::Int

    function Prediction(domain::AbstractDomain, comp_iterable...)
        @assert length(comp_iterable) > 0
        pred = new(domain, flatten(domain),
                   OrderedDict{Symbol, CompEval}(),
                   OrderedDict{Symbol, ReducerEval}(),
                   Symbol(""), 0)
        add_comps!(  pred, comp_iterable...)
        add_reducer!(pred, :sum1 => Reducer(sum_of_array))
        return pred
    end

    function Prediction(domain::AbstractDomain,
                        redpair::Pair{Symbol, Reducer}, comp_iterable...)
        @assert length(comp_iterable) > 0
        pred = new(domain, flatten(domain),
                   OrderedDict{Symbol, CompEval}(),
                   OrderedDict{Symbol, ReducerEval}(),
                   Symbol(""), 0)
        add_comps!(  pred, comp_iterable...)
        add_reducer!(pred, redpair)
        return pred
    end

    Prediction(domain::AbstractDomain, reducer::Reducer, comp_iterable...) =
        Prediction(domain, :reducer1 => reducer, comp_iterable...)
end

function add_comps!(pred::Prediction, comp_iterable...)
    for (cname, comp) in extract_components(comp_iterable...)
        @assert !haskey(pred.cevals, cname)  "Name $cname already exists"
        @assert !haskey(pred.revals, cname)  "Name $cname already exists"
        pred.cevals[cname] = CompEval(deepcopy(comp), pred.domain)
    end
end

sum_of_array( arg::Array) = arg
sum_of_array( args...) = .+(args...)
prod_of_array( arg::Array) = arg
prod_of_array(args...) = .*(args...)

add_reducer!(pred::Prediction, reducer::Reducer) =
    add_reducer!(pred, Symbol(:reducer, length(pred.revals)+1) => reducer)

function add_reducer!(pred::Prediction, redpair::Pair{Symbol, Reducer})
    rname = redpair[1]
    reducer = redpair[2]
    @assert !haskey(pred.cevals, rname)  "Name $rname already exists"
    haskey(pred.revals, rname)  &&  delete!(pred.revals, rname)
    if reducer.allargs
        append!(reducer.args, keys(pred.cevals))
        append!(reducer.args, keys(pred.revals))
    end

    args = Vector{Vector{Float64}}()
    for arg in reducer.args
        if haskey(pred.cevals, arg)
            push!(args, pred.cevals[arg].buffer)
        else
            push!(args, pred.revals[arg].buffer)
        end
    end

    (reducer.funct == sum)   &&  (reducer.funct = sum_of_array)
    (reducer.funct == prod)  &&  (reducer.funct = prod_of_array)
    eval = reducer.funct(args...)
    pred.revals[rname] = ReducerEval(args, reducer.funct, 1, eval)
    pred.rsel = rname
    evaluate(pred)
    return pred
end


function evaluate(pred::Prediction)
    for (cname, ceval) in pred.cevals
        evaluate_cached(ceval)
    end
    reduce(pred)
    return pred
end


function reduce(pred::Prediction)
    for (rname, reval) in pred.revals
        reval.counter += 1
        reval.buffer .= reval.funct(reval.args...)
    end
    pred.counter += 1
end


geteval(pred::Prediction) = geteval(pred, pred.rsel)
function geteval(pred::Prediction, name::Symbol)
    if haskey(pred.cevals, name)
        return pred.cevals[name].buffer
    else
        return pred.revals[name].buffer
    end
end


# ====================================================================
include("HashVector.jl")

# ====================================================================
# Model and ModelInternals structures
#
const PatchComp = HashVector{Float64}

struct ModelInternals
    cevals::OrderedDict{CompID, CompEval}
    params::OrderedDict{CompParamID, Parameter}
    par_indices::OrderedDict{CompID, Vector{Int}}
    pvalues::Vector{Float64}
    patched::Vector{Float64}
    patchcomps::OrderedDict{CompID, PatchComp}
    buffer::Vector{Float64}
end
ModelInternals() = ModelInternals(OrderedDict{CompID, CompEval}(),
                                  OrderedDict{CompParamID, Parameter}(),
                                  OrderedDict{CompID, Vector{Int}}(),
                                  Vector{Float64}(), Vector{Float64}(),
                                  OrderedDict{CompID, PatchComp}(),
                                  Vector{Float64}())

mutable struct Model
    meta::Vector{MDict}
    preds::Vector{Prediction}
    priv::ModelInternals
    patchfuncts::Vector{Function}

    function Model(v::Vector{Prediction})
        model = new(Vector{MDict}(), v, ModelInternals(), Vector{Function}())
        evaluate(model)
        return model
    end
    Model(p::Prediction) = Model([p])
    Model(args...) = Model([Prediction(args...)])
end


struct PredRef
    model::Model
    id::Int
end
deref(p::PredRef) = (p ⋄ :model).preds[p ⋄ :id]
Base.getproperty(p::PredRef, name::Symbol) = getproperty(deref(p), name)



function ModelInternals(model::Model)
    cevals = OrderedDict{CompID, CompEval}()
    params = OrderedDict{CompParamID, Parameter}()
    par_indices = OrderedDict{CompID, Vector{Int}}()

    ndata = 0
    i = 1
    for id in 1:length(model.preds)
        pred = model.preds[id]
        for (cname, ceval) in pred.cevals
            cid = CompID(id, cname)
            cevals[cid] = ceval
            par_indices[cid] = Vector{Int}()
            for (pid, par) in ceval.params
                params[CompParamID(cid, pid)] = par
                push!(par_indices[cid], i)
                i += 1
            end
            evaluate_cached(ceval)
        end
        reduce(pred)
        ndata += length(geteval(pred))
    end

    pvalues = getfield.(values(params), :val)
    patched = deepcopy(pvalues)

    patchcomps = OrderedDict{CompID, PatchComp}()
    i = 1
    for (cid, ceval) in cevals
        ipar  = OrderedDict{Symbol, Union{Int, Vector{Int}}}()
        for (pid, par) in ceval.params
            if pid.index == 0
                ipar[pid.name] = i
            elseif pid.index == 1
                ipar[pid.name] = [i]
            else
                push!(ipar[pid.name], i)
            end
            i += 1
        end
        patchcomps[cid] = HashVector(ipar, patched)
    end

    return ModelInternals(cevals, params, par_indices, pvalues, patched, patchcomps,
                          fill(NaN, ndata))
end


function evaluate(model::Model)
    @assert length(model.preds) >= 1
    model.priv = ModelInternals(model)
    quick_evaluate(model)

    # Prepare meta data
    while length(model.meta) < length(model.preds)
        push!(model.meta, MDict())
    end
    for id in 1:length(model.preds)
        pred = model.preds[id]
        meta = model.meta[id]

        haskey(meta, :components)  ||  (meta[:components] = MDict())
        for (cname, ceval) in pred.cevals
            haskey(meta[:components], cname)  ||  (meta[:components][cname] = MDict(:params => MDict()))
            for (pid, param) in ceval.params
                if pid.index >= 1
                    pname = Symbol(pid.name, "[", pid.index, "]")
                else
                    pname = pid.name
                end
                haskey(meta[:components][cname][:params], pname)  ||  (meta[:components][cname][:params][pname] = MDict())
            end
        end

        haskey(meta, :reducers)  ||  (meta[:reducers] = MDict())
        for (rname, reval) in pred.revals
            haskey(meta[:reducers], rname)  ||  (meta[:reducers][rname] = MDict())
        end
    end

    return model
end


# This is supposed to be called from `fit!`, not by user
function quick_evaluate(model::Model)
    model.priv.patched .= model.priv.pvalues  # copy all values by default
    for func in model.patchfuncts
        func(model.priv.patchcomps)
    end

    for (cid, ceval) in model.priv.cevals
        evaluate_cached(ceval, model.priv.patched[model.priv.par_indices[cid]])
    end
    for pred in model.preds
        reduce(pred)
    end
    nothing
end


function add!(model::Model, p::Prediction)
    push!(model.preds, p)
    evaluate(model)
end

add!(model::Model, args...) = add!(model[1], args...)

function add!(p::PredRef, comp_iterable...)
    @assert length(comp_iterable) > 0
    add_comps!(deref(p), comp_iterable...)
    evaluate(p ⋄ :model)
end

function add!(p::PredRef, reducer::Reducer, comp_iterable...)
    if length(comp_iterable) > 0
        add_comps!(deref(p), comp_iterable...)
    end
    add_reducer!(deref(p), reducer)
    evaluate(p ⋄ :model)
end

function add!(p::PredRef, redpair::Pair{Symbol, Reducer}, comp_iterable...)
    if length(comp_iterable) > 0
        add_comps!(deref(p), comp_iterable...)
    end
    add_reducer!(deref(p), redpair)
    evaluate(p ⋄ :model)
end


# ====================================================================
function patch!(func::Function, model::Model)
    push!(model.patchfuncts, func)
    evaluate(model)
    return model
end


isfixed(pref::PredRef, cname::Symbol) = (pref.cevals[cname].cfixed >= 1)
isfixed(model::Model, cname::Symbol) = isfixed(model[1], cname)

freeze(model::Model, cname::Symbol) = freeze(model[1], cname)
function freeze(pref::PredRef, cname::Symbol)
    @assert cname in keys(pref.cevals) "Component $cname is not defined on prediction $(pref ⋄ :id)"
    pref.cevals[cname].cfixed = 1
    evaluate(pref ⋄ :model)
    nothing
end


thaw(model::Model, cname::Symbol) = thaw(model[1], cname)
function thaw(pref::PredRef, cname::Symbol)
    @assert cname in keys(pref.cevals) "Component $cname is not defined on prediction $(pref ⋄ :id)"
    pref.cevals[cname].cfixed = 0
    evaluate(pref ⋄ :model)
    nothing
end


# ====================================================================
# Fit results
#
struct BestFitPar
    val::Float64
    unc::Float64
    fixed::Bool
    patched::Float64  # value after transformation
end

const BestFitComp = HashVector{BestFitPar}


struct BestFitResult
    preds::Vector{OrderedDict{Symbol, BestFitComp}}
    ndata::Int
    dof::Int
    cost::Float64
    status::Symbol      #:OK, :Warn, :Error
    log10testprob::Float64
    elapsed::Float64
end


struct BestFitPredRef
    result::BestFitResult
    id::Int
end
deref(p::BestFitPredRef) = (p ⋄ :result).preds[p ⋄ :id]
Base.getindex(p::BestFitPredRef, name::Symbol) = getindex(deref(p), name)
Base.iterate(p::BestFitPredRef, args...) = iterate(deref(p), args...)


# ====================================================================
function data1D(model::Model, data::Vector{T}) where T<:AbstractMeasures
    out = Vector{Measures_1D}()
    for i in 1:length(model.preds)
        pred = model.preds[i]
        @assert(length(data[i]) == length(geteval(pred)),
                "Length of dataset $i do not match corresponding model prediction.")
        push!(out, flatten(data[i], pred.domain))
    end
    return out
end


function residuals1d(model::Model, data1d::Vector{Measures_1D})
    c1 = 1
    for i in 1:length(model.preds)
        pred = model.preds[i]
        eval = geteval(pred)
        c2 = c1 + length(eval) - 1
        model.priv.buffer[c1:c2] .= ((eval .- data1d[i].val) ./ data1d[i].unc)
        c1 = c2 + 1
    end
    return model.priv.buffer
end


# ====================================================================
abstract type AbstractMinimizer end

using LsqFit
mutable struct lsqfit <: AbstractMinimizer
end

function minimize(minimizer::lsqfit, func::Function, params::Vector{Parameter})
    ndata = length(func(getfield.(params, :val)))
    bestfit = LsqFit.curve_fit((dummy, pvalues) -> func(pvalues),
                               1.:ndata, fill(0., ndata),
                               getfield.(params, :val),
                               lower=getfield.(params, :low),
                               upper=getfield.(params, :high))
    status = :Error
    (bestfit.converged)  &&  (status = :OK)
    error = LsqFit.margin_error(bestfit, 0.6827)
    return (status, getfield.(Ref(bestfit), :param), error)
end


using CMPFit;

mutable struct cmpfit <: AbstractMinimizer;
    config::CMPFit.Config;
    cmpfit() = new(CMPFit.Config());
end;

function minimize(minimizer::cmpfit, func::Function, params::Vector{Parameter});
    guess = getfield.(params, :val);
    low   = getfield.(params, :low);
    high  = getfield.(params, :high);
    parinfo = CMPFit.Parinfo(length(guess));
    for i in 1:length(guess);
        llow  = isfinite(low[i])   ?  1  :  0;
        lhigh = isfinite(high[i])  ?  1  :  0;
        parinfo[i].limited = (llow, lhigh);
        parinfo[i].limits  = (low[i], high[i]);
    end;
    bestfit = CMPFit.cmpfit((pvalues) -> func(pvalues),
                            guess, parinfo=parinfo, config=minimizer.config);
    return (:OK, getfield.(Ref(bestfit), :param), getfield.(Ref(bestfit), :perror));
end;


# ====================================================================
fit!(model::Model, data::T; kw...) where T<:AbstractMeasures =
    fit!(model, [data]; kw...)

function fit!(model::Model, data::Vector{T};
              only_id::Int=0,
              minimizer=lsqfit()) where T<:AbstractMeasures
    elapsedTime = Base.time_ns()
    evaluate(model)

    if only_id != 0
        for id in 1:length(model.preds)
            (id == only_id)  &&  continue
            pred = model.preds[id]
            for (cname, ceval) in pred.cevals
                ceval.cfixed += 1
            end
        end
    end

    free = Vector{Bool}()
    for (cpid, par) in model.priv.params
        push!(free, (!par.fixed)  &&  (model.priv.cevals[cpid.comp].cfixed == 0))
    end
    ifree = findall(free)
    @assert length(ifree) > 0 "No free parameter in the model"

    # Flatten empirical data
    data1d = data1D(model, data)

    # Evaluate normalized residuals starting from free parameter values
    function pval2resid(pvalues_free::Vector{Float64})
        model.priv.pvalues[ifree] .= pvalues_free  # update parameter values
        quick_evaluate(model)
        return residuals1d(model, data1d)
    end

    (status, best_val, best_unc) = minimize(minimizer, pval2resid,
                                            collect(values(model.priv.params))[ifree])
    model.priv.pvalues[ifree] .= best_val

    # Copy best fit values back into components.  This is needed since
    # the evaluated components are stored in the Model, not in
    # BestFitResult, hence I do this to maintain a coherent status.
    setfield!.(values(model.priv.params), :val, model.priv.pvalues)
    uncerts = fill(NaN, length(model.priv.pvalues))
    uncerts[ifree] .= best_unc

    # Prepare output
    quick_evaluate(model)  # ensure best fit values are used

    bfpars = Vector{BestFitPar}()
    i = 1
    for (cpid, par) in model.priv.params
        push!(bfpars, BestFitPar(model.priv.pvalues[i], uncerts[i],
                                 !(i in ifree), model.priv.patched[i]))
        i += 1
    end

    preds = fill(OrderedDict{Symbol, BestFitComp}(), length(model.preds))
    for (cid, ceval) in model.priv.cevals
        preds[cid.id][cid.name] = BestFitComp(getfield(model.priv.patchcomps[cid], :dict), bfpars)
    end

    cost = sum(abs2, model.priv.buffer)
    dof = length(model.priv.buffer) - length(ifree)

    result = BestFitResult(preds, length(model.priv.buffer), dof, cost, status,
                           logccdf(Chisq(dof), cost) * log10(exp(1)),
                           float(Base.time_ns() - elapsedTime) / 1.e9)

    if only_id != 0
        for id in 1:length(model.preds)
            (id == only_id)  &&  continue
            pred = model.preds[id]
            for (cname, ceval) in pred.cevals
                (ceval.cfixed > 1)  &&  (ceval.cfixed = 1)
            end
        end
        evaluate(model)
    end
    return result
end


# ====================================================================
# User interface

##
(p::PredRef)() = geteval(deref(p))
(p::PredRef)(name::Symbol) = geteval(deref(p), name)
function (m::Model)()
    evaluate(m)
    m[1]()
end
function (m::Model)(name::Symbol)
    evaluate(m)
    m[1](name)
end

##
Base.keys(m::Model) = keys(m[1])
Base.keys(p::PredRef) = keys(p.cevals)
Base.haskey(p::PredRef, name::Symbol) = haskey(p.cevals, name)
Base.haskey(m::Model, name::Symbol) = haskey(m[1], name)

##
Base.getindex(pref::PredRef, cname::Symbol) = pref.cevals[cname].comp
Base.getindex(m::Model, id::Int) = PredRef(m, id)
Base.getindex(m::Model, cname::Symbol) = m[1][cname]
Base.getindex(res::BestFitResult, id::Int) = BestFitPredRef(res, id)
Base.getindex(res::BestFitResult, cname::Symbol) = res[1][cname]

##
domain(pref::PredRef; dim::Int=1) = pref.orig_domain[dim]
domain(m::Model; dim::Int=1) = domain(m[1], dim=dim)


# ====================================================================
include("todict.jl")
include("show.jl")

end
