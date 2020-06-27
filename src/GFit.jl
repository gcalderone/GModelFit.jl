module GFit

using Printf, PrettyTables
using Statistics, Distributions
using DataStructures
using LsqFit
using JSON
using ExprTools

import Base.push!
import Base.show
import Base.ndims
import Base.size
import Base.length
import Base.getindex
import Base.reshape
import Base.propertynames
import Base.getproperty
import Base.setindex!
import Base.iterate
import Base.dump


export Domain, CartesianDomain, Measures,
    Prediction, Reducer, @reducer, add!, domain,
    Model, constraint!, evaluate, parindex, thaw, freeze, fit!,
    savelog

include("domain.jl")

# ====================================================================
# Parameter
#
mutable struct Parameter
    val::Float64
    low::Float64              # lower limit value
    high::Float64             # upper limit value
    step::Float64
    free::Bool
    Parameter(value::Number) = new(float(value), -Inf, +Inf, NaN, true)
end

# ====================================================================
# A *component* is a generic implementation of a constituent part of a
# model.
#
# A component must inherit `AbstractComponent`, and implement the
# `ceval_data` and `evaluate` methods
abstract type AbstractComponent end


function getparams(comp::AbstractComponent)
    params = OrderedDict{Tuple{Symbol,Int}, Parameter}()
    for pname in fieldnames(typeof(comp))
        par = getfield(comp, pname)
        if isa(par, Parameter)
            params[(pname, 0)] = par
        elseif isa(par, Vector{Parameter})
            for i in 1:length(par)
                params[(pname, i)] = par[i]
            end
        end
    end
    params
end

# ====================================================================
# CompEval: a wrapper for a component evaluated on a specific domain
#
mutable struct CompEval{TDomain <: AbstractDomain, TComp <: AbstractComponent}
    domain::TDomain
    comp::TComp
    params::OrderedDict{Tuple{Symbol,Int}, Parameter}
    cdata
    counter::Int
    lastvalues::Vector{Float64}
    eval::Vector{Float64}
    ipar::Vector{Int}  # handled by Model

    function CompEval(domain::AbstractDomain, comp::AbstractComponent)
        params = getparams(comp)
        (cdata, len) = ceval_data(domain, comp)
        return new{typeof(domain), typeof(comp)}(
            domain, comp, params, cdata, 0,
            fill(NaN, length(params)),
            fill(NaN, len), Vector{Int}())
    end
end


# This is called to `update` to distinguish it from component's `evaluate`.
update(c::CompEval) = update(c, [par.val for par in values(c.params)])
function update(c::CompEval, pvalues::Vector{Float64})
    @assert length(c.params) == length(pvalues)

    # Do we actually need a new evaluation?
    if (any(c.lastvalues .!= pvalues)  ||  (c.counter == 0))
        c.lastvalues .= pvalues
        c.counter += 1
        @assert all(.!isnan.(pvalues))
        evaluate(c, pvalues...)
    end
    return c.eval
end


# ====================================================================
# Component fall back methods
ceval_data(domain::AbstractDomain, comp::AbstractComponent) =
    error("Component " * string(typeof(comp)) * " must implement its own method for `ceval_data`.")

evaluate(c::CompEval{TDomain, TComp}, args...) where {TDomain, TComp} =
    error("Component " * string(TComp) * " must implement its own method for `evaluate`.")


# ====================================================================
# Built-in components
#
include("components/SimplePar.jl")
include("components/FuncWrap.jl")
include("components/OffsetSlope.jl")
include("components/Gaussian.jl")


# ====================================================================
# Parse a user defined structure or dictionary to extract all
# components
function extract_components(things::AbstractDict)
    out = OrderedDict{Symbol, AbstractComponent}()
    for (name, comp) in things
        isa(comp, Number)  &&  (comp = SimplePar(comp))
        @assert isa(name, Symbol)
        @assert isa(comp, AbstractComponent)
        out[name] = comp
    end
    return out
end

function extract_components(things::Vararg{Pair})
    out = OrderedDict{Symbol, AbstractComponent}()
    for thing in things
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
Reducer(f::Function, args::Vector{Symbol}) =  Reducer(f, false, args)

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
    eval::Vector{Float64}
end


# ====================================================================
# A model prediction suitable to be compared to experimental data
mutable struct Prediction
    domain::AbstractDomain
    cevals::OrderedDict{Symbol, CompEval}
    revals::OrderedDict{Symbol, ReducerEval}
    rsel::Symbol
    counter::Int

    function Prediction(domain::AbstractDomain, things...)
        pred = new(domain,
                  OrderedDict{Symbol, CompEval}(),
                  OrderedDict{Symbol, ReducerEval}(),
                  Symbol(""), 0)
        add!(pred, things...)
        return pred
    end

    Prediction(domain::AbstractDomain, reducer::Reducer, things...) =
        Prediction(domain, :_1 => reducer, things...)

    function Prediction(domain::AbstractDomain, redpair::Pair{Symbol, Reducer}, things...)
        pred = Prediction(domain, things...)
        add!(pred, redpair)
        return pred
    end
end

function add!(pred::Prediction, things...)
    for (cname, comp) in extract_components(things...)
        @assert !haskey(pred.cevals, cname)  "Name $cname already exists"
        @assert !haskey(pred.revals, cname)  "Name $cname already exists"
        pred.cevals[cname] = CompEval(pred.domain, comp)
    end
end

sum_of_array(args...) = .+(args...)

add!(pred::Prediction, reducer::Reducer) =
    add!(pred, Symbol(:_, length(pred.revals)+1) => reducer)

function add!(pred::Prediction, redpair::Pair{Symbol, Reducer})
    rname = redpair[1]
    reducer = redpair[2]
    if rname == Symbol("")
        rname = Symbol(:_, length(pred.revals)+1)
    end
    @assert !haskey(pred.cevals, rname)  "Name $rname already exists"
    haskey(pred.revals, rname)  &&  delete!(pred.revals, rname)
    if reducer.allargs
        append!(reducer.args, keys(pred.cevals))
        append!(reducer.args, keys(pred.revals))
    end

    args = Vector{Vector{Float64}}()
    for arg in reducer.args
        if haskey(pred.cevals, arg)
            push!(args, pred.cevals[arg].eval)
        else
            push!(args, pred.revals[arg].eval)
        end
    end

    (reducer.funct == sum)  &&  (reducer.funct = sum_of_array)
    eval = reducer.funct(args...)
    pred.revals[rname] = ReducerEval(args, reducer.funct, eval)
    pred.rsel = rname
    evaluate(pred)
    return pred
end

add!(pred::Prediction, reducer::Reducer, things...) =
    add!(pred, Symbol(:_, length(pred.revals)+1) => reducer, things...)

function add!(pred::Prediction, redpair::Pair{Symbol, Reducer}, things...)
    add!(pred, things...)
    add!(pred, redpair)
end

function evaluate(pred::Prediction)
    for (name, ceval) in pred.cevals
        update(ceval)
    end
    reduce(pred)
    return pred
end

function reduce(pred::Prediction)
    for (rname, reval) in pred.revals
        reval.eval .= reval.funct(reval.args...)
    end
    pred.counter += 1
end

(pred::Prediction)() = pred(pred.rsel)
function (pred::Prediction)(name::Symbol)
    if haskey(pred.cevals, name)
        return pred.cevals[name].eval
    else
        return pred.revals[name].eval
    end
end
Base.getindex(pred::Prediction, cname::Symbol) = pred.cevals[cname].comp
domain(pred::Prediction, dim::Int=1) = pred.domain[dim]


# ====================================================================
# Global model, actually a collection of `Prediction`s.
mutable struct Model
    preds::Vector{Prediction}
    comps::OrderedDict{Symbol, AbstractComponent}
    cfree::OrderedDict{Symbol, Bool}
    params::OrderedDict{Tuple{Symbol, Symbol, Int}, Parameter}
    pvalues::Vector{Float64}
    actual::Vector{Float64}
    buffer::Vector{Float64}
    constraints::Vector{Function}

    function Model(v::Vector{Prediction})
        model = new(v, OrderedDict{Symbol, AbstractComponent}(),
                    OrderedDict{Symbol, Bool}(),
                    OrderedDict{Tuple{Symbol, Symbol, Int}, Parameter}(),
                    Vector{Float64}(), Vector{Float64}(), Vector{Float64}(), Vector{Function}())
        evaluate(model)
        return model
    end
    Model(p::Prediction) = Model([p])
    Model(args...) = Model([Prediction(args...)])
end

function evaluate(model::Model)
    @assert length(model.preds) >= 1

    # Save list of previously free components
    cfree = deepcopy(model.cfree)

    # Collect components and parameters
    empty!(model.comps)
    empty!(model.cfree)
    empty!(model.params)
    for pred in model.preds
        for (cname, ceval) in pred.cevals
            model.comps[cname] = ceval.comp
            model.cfree[cname] = get(cfree, cname, true)
            for (pname, par) in ceval.params
                cpname = (cname, pname[1], pname[2])
                model.params[cpname] = par
            end
        end
    end

    # Populate CompEval.ipar and evaluate all predictions
    ndata = 0
    cpnames = keys(model.params)
    for pred in model.preds
        for (cname, ceval) in pred.cevals
            empty!(ceval.ipar)
            for (pname, par) in ceval.params
                cpname = (cname, pname[1], pname[2])
                push!(ceval.ipar, findfirst(cpnames .== Ref(cpname)))
            end
            update(ceval)
        end
        reduce(pred)
        ndata += length(pred())
    end

    model.pvalues = [par.val for par in values(model.params)]
    model.actual = deepcopy(model.pvalues)
    model.buffer = Vector{Float64}(undef, ndata)
    quick_evaluate(model)
    return model
end

# This is supposed to be called from `fit!`, not by user
function quick_evaluate(model::Model)
    model.actual .= model.pvalues  # copy all values by default
    for func in model.constraints
        func(model, model.pvalues, model.actual)
    end

    for pred in model.preds
        for (cname, ceval) in pred.cevals
            update(ceval, model.actual[ceval.ipar])
        end
    end
    for pred in model.preds
        reduce(pred)
    end
    nothing
end


function add!(model::Model, p::Prediction)
    push!(model.preds, p)
    evaluate(model)
    return model
end


function constraint!(model::Model, func::Function)
    push!(model.constraints, func)
    evaluate(model)
    return model
end


(m::Model)() = m(1)
(m::Model)(i::Int) = m.preds[i]()
(m::Model)(i::Int, name::Symbol) = m.preds[i](name)
Base.getindex(m::Model, cname::Symbol) = m.comps[cname]
domain(m::Model, i::Int=1, dim::Int=1) = domain(m.preds[i], dim)

parindex(model::Model, cname::Symbol, pname::Symbol, i::Int=0) =
    findfirst(keys(model.params) .== Ref((cname, pname, i)))

function freeze(model::Model, cname::Symbol)
    evaluate(model)
    @assert cname in keys(model.cfree) "Component $cname is not defined"
    model.cfree[cname] = false
    model
end

function thaw(model::Model, cname::Symbol)
    evaluate(model)
    @assert cname in keys(model.cfree) "Component $cname is not defined"
    model.cfree[cname] = true
    model
end


# ====================================================================
# Fit results
#
struct BestFitPar
    val::Float64
    unc::Float64
    free::Bool
    actual::Float64  # value after transformation
end

struct BestFitComp
    params::OrderedDict{Symbol, Union{BestFitPar, Vector{BestFitPar}}}
    BestFitComp() = new(OrderedDict{Symbol, Union{BestFitPar, Vector{BestFitPar}}}())
end

Base.propertynames(comp::BestFitComp) = keys(getfield(comp, :params))
Base.getproperty(comp::BestFitComp, p::Symbol) = getfield(comp, :params)[p]
Base.getindex(comp::BestFitComp, p::Symbol) = getfield(comp, :params)[p]
Base.length(comp::BestFitComp) = length(getfield(comp, :params))
Base.iterate(comp::BestFitComp, args...) = iterate(getfield(comp, :params), args...)
Base.setindex!(comp::BestFitComp, x, p::Symbol) = getfield(comp, :params)[p] = x


struct BestFitResult
    comps::OrderedDict{Symbol, BestFitComp}
    ndata::Int
    dof::Int
    cost::Float64
    status::Symbol      #:OK, :Warn, :Error
    log10testprob::Float64
    elapsed::Float64
end

Base.getindex(res::BestFitResult, cname::Symbol) = res.comps[cname]

# ====================================================================
function data1D(model::Model, data::Vector{T}) where T<:AbstractMeasures
    out = Vector{Measures_1D}()
    for i in 1:length(model.preds)
        pred = model.preds[i]
        @assert(length(data[i]) == length(pred()),
                "Length of dataset $i do not match corresponding model prediction.")
        push!(out, flatten(data[i], pred.domain))
    end
    return out
end


function residuals1d(model::Model, data1d::Vector{Measures_1D})
    c1 = 1
    for i in 1:length(model.preds)
        pred = model.preds[i]
        eval = pred()
        c2 = c1 + length(eval) - 1
        model.buffer[c1:c2] .= ((eval .- data1d[i].val) ./ data1d[i].unc)
        c1 = c2 + 1
    end
    return model.buffer
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


macro with_CMPFit()
    return esc(:(
        using CMPFit;
        import GFit.minimize;

        mutable struct cmpfit <: GFit.AbstractMinimizer;
        config::CMPFit.Config;
        cmpfit() = new(CMPFit.Config());
        end;

        function minimize(minimizer::cmpfit, func::Function, params::Vector{GFit.Parameter});
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
    ))
end


fit!(model::Model, data::T; kw...) where T<:AbstractMeasures =
    fit!(model, [data]; kw...)

function fit!(model::Model, data::Vector{T};
              minimizer=lsqfit()) where T<:AbstractMeasures
    elapsedTime = Base.time_ns()
    evaluate(model)

    free = Vector{Bool}()
    for (cpname, par) in model.params
        push!(free, par.free  &&  model.cfree[cpname[1]])
    end
    ifree = findall(free)
    @assert length(ifree) > 0 "No free parameter in the model"

    # Flatten empirical data
    data1d = data1D(model, data)

    # Evaluate normalized residuals starting from free parameter values
    function pval2resid(pvalues_free::Vector{Float64})
        model.pvalues[ifree] .= pvalues_free  # update parameter values
        quick_evaluate(model)
        return residuals1d(model, data1d)
    end

    (status, best_val, best_unc) = minimize(minimizer, pval2resid,
                                            collect(values(model.params))[ifree])

    model.pvalues[ifree] .= best_val
    setfield!.(values(model.params), :val, model.pvalues)
    uncerts = fill(NaN, length(model.pvalues))
    uncerts[ifree] .= best_unc

    # Prepare output
    quick_evaluate(model)  # ensure best fit values are used
    comps = OrderedDict{Symbol, BestFitComp}()
    for cname in keys(model.comps)
        comps[cname] = BestFitComp()
    end
    i = 1
    for (cpname, par) in model.params
        cname = cpname[1]
        pname = cpname[2]
        parid = cpname[3]
        bfpar = BestFitPar(model.pvalues[i], uncerts[i],
                           (i in ifree), model.actual[i])
        if parid == 0
            comps[cname][pname] = bfpar
        else
            if parid == 1
                comps[cname][pname] = [bfpar]
            else
                push!(comps[cname][pname], bfpar)
            end
        end
        i += 1
    end
    cost = sum(abs2, model.buffer)
    dof = length(model.buffer) - length(ifree)

    result = BestFitResult(comps, length(model.buffer), dof, cost, status,
                           logccdf(Chisq(dof), cost) * log10(exp(1)),
                           float(Base.time_ns() - elapsedTime) / 1.e9)
    return result
end

include("show.jl")
include("viewer.jl")

end
