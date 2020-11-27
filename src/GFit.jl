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
import Base.getindex
import Base.reshape
import Base.propertynames
import Base.getproperty
import Base.setproperty!
import Base.iterate


export Domain, CartesianDomain, Measures,
    Prediction, Reducer, @reducer, add!, domain,
    Model, patch!, evaluate, thaw, freeze, fit!,
    metadict, savelog

const MDict = OrderedDict{Symbol, Any}


include("domain.jl")


# ====================================================================
# Parameter
#
mutable struct Parameter
    meta::MDict
    val::Float64
    low::Float64              # lower limit value
    high::Float64             # upper limit value
    step::Float64
    fixed::Bool
    Parameter(value::Number) = new(MDict(), float(value), -Inf, +Inf, NaN, false)
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
mutable struct CompEval{TComp <: AbstractComponent, TDomain <: AbstractDomain}
    comp::TComp
    domain::TDomain
    params::OrderedDict{Tuple{Symbol,Int}, Parameter}
    cdata
    counter::Int
    lastvalues::Vector{Float64}
    buffer::Vector{Float64}
    ipar::Vector{Int}  # handled by Model

    function CompEval(comp::AbstractComponent, domain::AbstractDomain)
        params = getparams(comp)
        cdata  = compeval_cdata(comp, domain)
        buffer = compeval_array(comp, domain)
        return new{typeof(comp), typeof(domain)}(
            comp, domain, params, cdata, 0,
            fill(NaN, length(params)), 
            buffer, Vector{Int}())
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
compeval_cdata(comp::AbstractComponent, domain::AbstractDomain) =
    error("Component " * string(typeof(comp)) * " must implement its own method for `compeval_cdata`.")

compeval_array(comp::AbstractComponent, domain::AbstractDomain) =
    error("Component " * string(typeof(comp)) * " must implement its own method for `compeval_array`.")

evaluate(c::CompEval{TComp, TDomain}, args...) where {TComp, TDomain} =
    error("Component " * string(TComp) * " must implement its own method for `evaluate`.")

evaluate(c::CompEval{TComp, TDomain}) where {TComp, TDomain} = 
    error("Component " * string(TComp) * " must implement its own method for `evaluate`.")


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
    meta::MDict
    orig_domain::AbstractDomain
    domain::AbstractLinearDomain
    cevals::OrderedDict{Symbol, CompEval}
    revals::OrderedDict{Symbol, ReducerEval}
    rsel::Symbol
    counter::Int

    function Prediction(domain::AbstractDomain, comp_iterable...)
        @assert length(comp_iterable) > 0
        pred = new(MDict(), domain, flatten(domain),
                   OrderedDict{Symbol, CompEval}(),
                   OrderedDict{Symbol, ReducerEval}(),
                   Symbol(""), 0)
        add_comps!(  pred, comp_iterable...)
        add_reducer!(pred, :autogen1 => Reducer(sum_of_array))
        return pred
    end

    function Prediction(domain::AbstractDomain,
                        redpair::Pair{Symbol, Reducer}, comp_iterable...)
        @assert length(comp_iterable) > 0
        pred = new(MDict(), domain, flatten(domain),
                   OrderedDict{Symbol, CompEval}(),
                   OrderedDict{Symbol, ReducerEval}(),
                   Symbol(""), 0)
        add_comps!(  pred, comp_iterable...)
        add_reducer!(pred, redpair)
        return pred
    end

    Prediction(domain::AbstractDomain, reducer::Reducer, comp_iterable...) =
        Prediction(domain, :autogen1 => reducer, comp_iterable...)
end

function add_comps!(pred::Prediction, comp_iterable...)
    for (cname, comp) in extract_components(comp_iterable...)
        @assert !haskey(pred.cevals, cname)  "Name $cname already exists"
        @assert !haskey(pred.revals, cname)  "Name $cname already exists"
        pred.cevals[cname] = CompEval(comp, pred.domain)
    end
end

sum_of_array( arg::Array) = arg
sum_of_array( args...) = .+(args...)
prod_of_array( arg::Array) = arg
prod_of_array(args...) = .*(args...)

add_reducer!(pred::Prediction, reducer::Reducer) =
    add_reducer!(pred, Symbol(:autogen, length(pred.revals)+1) => reducer)

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
    for (name, ceval) in pred.cevals
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
struct PatchComp
    pvalues::Vector{Float64}
    ipar::OrderedDict{Symbol, Vector{Int}}
end

function Base.getproperty(comp::PatchComp, pname::Symbol)
    v = getfield(comp, :pvalues)
    i = getfield(comp, :ipar)[pname]
    (length(i) == 1)  &&  (i = i[1])
    return v[i]
end

function Base.setproperty!(comp::PatchComp, pname::Symbol, value::Real)
    v = getfield(comp, :pvalues)
    i = getfield(comp, :ipar)[pname]
    @assert length(i) == 1
    v[i[1]] = value
end



# ====================================================================
# Global model, actually a collection of `Prediction`s.
mutable struct Model
    meta::MDict
    preds::Vector{Prediction}
    comps::OrderedDict{Symbol, AbstractComponent}
    cfixed::OrderedDict{Symbol, Bool}
    params::OrderedDict{Tuple{Symbol, Symbol, Int}, Parameter}
    patchbay::OrderedDict{Symbol, PatchComp}
    pvalues::Vector{Float64}
    patched::Vector{Float64}
    buffer::Vector{Float64}
    patchfuncts::Vector{Function}

    function Model(v::Vector{Prediction})
        model = new(MDict(), v, OrderedDict{Symbol, AbstractComponent}(),
                    OrderedDict{Symbol, Bool}(),
                    OrderedDict{Tuple{Symbol, Symbol, Int}, Parameter}(),
                    OrderedDict{Symbol, PatchComp}(),
                    Vector{Float64}(), Vector{Float64}(), Vector{Float64}(), Vector{Function}())
        evaluate(model)
        return model
    end
    Model(p::Prediction) = Model([p])
    Model(args...) = Model([Prediction(args...)])
end

function evaluate(model::Model)
    @assert length(model.preds) >= 1

    # Save list of previously fixed components
    cfixed = deepcopy(model.cfixed)

    # Collect components and parameters
    empty!(model.comps)
    empty!(model.cfixed)
    empty!(model.params)
    for pred in model.preds
        for (cname, ceval) in pred.cevals
            model.comps[cname] = ceval.comp
            model.cfixed[cname] = get(cfixed, cname, false)
            for (pname, par) in ceval.params
                cpname = (cname, pname[1], pname[2])
                model.params[cpname] = par
            end
        end
    end

    # Prepare vectors of parameter values (pvalues) and "patched"
    # values (patched)
    model.pvalues = [par.val for par in values(model.params)]
    model.patched = deepcopy(model.pvalues)

    # Populate CompEval.ipar and evaluate all predictions
    empty!(model.patchbay)
    ndata = 0
    cpnames = keys(model.params)
    for pred in model.preds
        for (cname, ceval) in pred.cevals
            empty!(ceval.ipar)
            liveipar = OrderedDict{Symbol, Vector{Int}}()

            for (pname, par) in ceval.params
                cpname = (cname, pname[1], pname[2])
                i = findall(cpnames .== Ref(cpname))
                @assert length(i) == 1
                push!(ceval.ipar, i[1])

                haskey(liveipar, pname[1])  ||  (liveipar[pname[1]] = Vector{Int}())
                push!( liveipar[pname[1]], i[1])
            end
            evaluate_cached(ceval)

            model.patchbay[cname] = PatchComp(model.patched, liveipar)
        end
        reduce(pred)
        ndata += length(geteval(pred))
    end

    model.buffer = Vector{Float64}(undef, ndata)
    quick_evaluate(model)
    return model
end

# This is supposed to be called from `fit!`, not by user
function quick_evaluate(model::Model)
    model.patched .= model.pvalues  # copy all values by default
    for func in model.patchfuncts
        func(model.patchbay)
    end

    for pred in model.preds
        for (cname, ceval) in pred.cevals
            evaluate_cached(ceval, model.patched[ceval.ipar])
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
end


function add!(model::Model, reducer::Reducer, comp_iterable...; id::Int=1)
    @assert length(comp_iterable) > 0
    add_comps!(  model.preds[id], comp_iterable...)
    add_reducer!(model.preds[id], reducer)
    evaluate(model)
end


function add!(model::Model, redpair::Pair{Symbol, Reducer}, comp_iterable...; id::Int=1)
    @assert length(comp_iterable) > 0
    add_comps!(  model.preds[id], comp_iterable...)
    add_reducer!(model.preds[id], redpair)
    evaluate(model)
end


# ====================================================================
function patch!(func::Function, model::Model)
    push!(model.patchfuncts, func)
    evaluate(model)
    return model
end


function freeze(model::Model, cname::Symbol)
    evaluate(model)
    @assert cname in keys(model.cfixed) "Component $cname is not defined"
    model.cfixed[cname] = true
    model
end


function thaw(model::Model, cname::Symbol)
    evaluate(model)
    @assert cname in keys(model.cfixed) "Component $cname is not defined"
    model.cfixed[cname] = false
    model
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

struct BestFitComp
    params::OrderedDict{Symbol, Union{BestFitPar, Vector{BestFitPar}}}
    BestFitComp() = new(OrderedDict{Symbol, Union{BestFitPar, Vector{BestFitPar}}}())
end

Base.length(comp::BestFitComp) = length(getfield(comp, :params))
Base.iterate(comp::BestFitComp, args...) = iterate(getfield(comp, :params), args...)


struct BestFitResult
    comps::OrderedDict{Symbol, BestFitComp}
    ndata::Int
    dof::Int
    cost::Float64
    status::Symbol      #:OK, :Warn, :Error
    log10testprob::Float64
    elapsed::Float64
end


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
              id::Int=0,
              minimizer=lsqfit()) where T<:AbstractMeasures
    elapsedTime = Base.time_ns()
    evaluate(model)

    if id != 0
        origcfixed = deepcopy(model.cfixed)
        for (cname, comp) in model.comps
            if !haskey(model.preds[id].cevals, cname)
                model.cfixed[cname] = true
            end
        end
    end

    free = Vector{Bool}()
    for (cpname, par) in model.params
        push!(free, (!par.fixed)  &&  (!model.cfixed[cpname[1]]))
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
    # Copy best fit values back into components.  This is needed since
    # the evaluated components are stored in the Model, not in
    # BestFitResult, hence I do this to maintain a coherent status.
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
                           !(i in ifree), model.patched[i])
        if parid == 0
            #Base.setindex!(comp::BestFitComp, x, p::Symbol) = getfield(comp, :params)[p] = x
            getfield(comps[cname], :params)[pname] = bfpar
        else
            if parid == 1
                getfield(comps[cname], :params)[pname] = [bfpar]
            else
                push!(getfield(comps[cname], :params)[pname], bfpar)
            end
        end
        i += 1
    end
    cost = sum(abs2, model.buffer)
    dof = length(model.buffer) - length(ifree)

    result = BestFitResult(comps, length(model.buffer), dof, cost, status,
                           logccdf(Chisq(dof), cost) * log10(exp(1)),
                           float(Base.time_ns() - elapsedTime) / 1.e9)

    if id != 0
        for cname in keys(origcfixed)
            model.cfixed[cname] = origcfixed[cname]
        end
    end

    return result
end


# ====================================================================
# User interface
#
Base.propertynames(comp::BestFitComp) = keys(getfield(comp, :params))
Base.getproperty(comp::BestFitComp, p::Symbol) = getfield(comp, :params)[p]

##
(m::Model)(; id::Int=1) = geteval(m.preds[id])
(m::Model)(name::Symbol; id::Int=1) = geteval(m.preds[id], name)

##
Base.getindex(m::Model, cname::Symbol) = m.comps[cname]
Base.getindex(res::BestFitResult, cname::Symbol) = res.comps[cname]

##
domain(pred::Prediction; dim::Int=1) = pred.orig_domain[dim]
domain(m::Model; id::Int=1, dim::Int=1) = m.preds[id].orig_domain[dim]


##
metadict(d::AbstractData) = d.meta
metadict(param::Parameter) = param.meta
metadict(m::Model; id::Int=1) = m.preds[id].meta

function metadict(m::Model, name::Symbol)
    if haskey(m.comps, name)
        haskey(m.meta, name)  ||  (m.meta[name] = MDict())
        comp = m.comps[name]
        if name in fieldnames(typeof(comp))
            field = getfield(comp, name)
            if isa(field, MDict)
                m.meta[name] = field
            end
        end

        return m.meta[name]
    else
        for pred in m.preds
            if haskey(pred.revals, name)
                haskey(m.meta, name)  ||  (m.meta[name] = MDict())
                return m.meta[name]
            end
        end
    end
    error("Name $name is not defined")
end

include("todict.jl")
include("show.jl")

end
