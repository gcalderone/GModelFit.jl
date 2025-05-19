import ForwardDiff: Dual

"""
    CompEval(comp::AbstractComponent, domain::AbstractDomain)

A container for a component to be evaluated on a specific domain.

# Fields:
 - `comp::AbstractComponent`: the wrapped component;
 - `domain::AbstractDomain`: the domain where the component is supposed to be evaluated;
 - `counter::Int`: the number of times the component has been evaluated since creatio of the `CompEval` object;
 - `lastparvalues::Vector{Float64}`: the parameter values used in the last evaluation.  A call to `evaluate()` with the same values stored in `lastparvalues` will not result in a new evaluation;
 - `deps::Vector{Vector{Float64}}`: the evaluation buffers of all dependencies;
 - `buffer::Vector{Float64}`: the buffer to store the outcome of the component.
"""
mutable struct CompEvalT{T <: Real}
    counter::Int
    lastparvalues::Vector{Union{T, Float64}}
    deps::Vector{Vector{Union{T, Float64}}}
    buffer::Vector{T}

    CompEvalT{T}(npar::Int, nres::Int) where T <: Real =
        new{T}(0,
            Vector{       Union{T, Float64}}( undef, npar),
            Vector{Vector{Union{T, Float64}}}(),
            Vector{             T}(           undef, nres))
end

mutable struct CompEval{TComp <: AbstractComponent, TDomain <: AbstractDomain}
    comp::TComp
    domain::TDomain
    tpar::CompEvalT{Float64}
    tparad::CompEvalT{Dual}

    function CompEval(comp::AbstractComponent, domain::AbstractDomain)
        prepare!(comp, domain)
        npar = length(getparams(comp))
        nres = result_length(comp, domain)
        return new{typeof(comp), typeof(domain)}(
            comp, domain,
            CompEvalT{Float64}(npar, nres),
            CompEvalT{Dual}(   npar, nres))
    end
end


"""
    evaluate(ceval::CompEval, pvalues::Vector{Float64})

Evaluate a component using the provided parameter values.  Outcomes shall be stored in the `CompEval.buffer` vector.
"""
evaluate(::AbstractComponent, ::AbstractDomain, args...) =
    error("No evaluate method implemented for $(TComp), $(TDomain)")


"""
    evaluate_comp!(ceval::CompEval, pvalues::Vector{Float64})

Update a `CompEval` structure using the provided parameter values.

The component is actually evaluated if one of the following applies:
- the component has never been evaluated;
- the component has at least one dependency (whose evaluation may have changed since its last evaluation);
- at least one parameter value has changed since last evaluation.

If none of the above applies, no evaluation occurs.
"""
function evaluate(ceval::CompEval, pvalues::AbstractVector{Float64})
    if length(ceval.tpar.deps) > 0
        evaluate!(ceval.comp, ceval.domain, ceval.tpar.buffer, ceval.tpar.deps, pvalues...)
        ceval.tpar.counter += 1
    elseif (ceval.tpar.counter == 0)  ||  any(ceval.tpar.lastparvalues .!= pvalues)
        evaluate!(ceval.comp, ceval.domain, ceval.tpar.buffer, pvalues...)
        ceval.tpar.lastparvalues .= pvalues
        ceval.tpar.counter += 1
    end
    return ceval.tpar.buffer
end

function evaluate(ceval::CompEval, pvalues::AbstractVector)
    @info typeof(pvalues)
    if length(ceval.tparad.deps) > 0
        evaluate!(ceval.comp, ceval.domain, ceval.tparad.buffer, ceval.tparad.deps, pvalues...)
        ceval.tparad.counter += 1
    elseif (ceval.tparad.counter == 0)  ||  any(ceval.tparad.lastparvalues .!= pvalues)
        evaluate!(ceval.comp, ceval.domain, ceval.tparad.buffer, pvalues...)
        ceval.tparad.lastparvalues .= pvalues
        ceval.tparad.counter += 1
    end
    return ceval.tparad.buffer
end


# Evaluate component on the given domain.  Parameter values are the
# ones stored in the component unless a custom value is provided via a
# keyword.
function (comp::AbstractComponent)(domain::AbstractDomain; kws...)
    @assert length(dependencies(comp)) == 0 "Can't evaluate a component with dependencies as a stand-alone one."
    ceval = CompEval(comp, domain)
    pvalues = OrderedDict([(pname, par.val) for (pname, par) in getparams(comp)])
    for (pname, pval) in kws
        if pname in keys(pvalues)
            pvalues[pname] = pval
        else
            @warn "$pname is not a parameter name for $(typeof(comp)). Valid names are: " * join(string.(keys(pvalues)), ", ")
        end
    end
    return evaluate(ceval, collect(values(pvalues)))
end


# ====================================================================
# Built-in components
include("components/FComp.jl")
include("components/OffsetSlope.jl")
include("components/Polynomial.jl")
include("components/Gaussian.jl")
include("components/Lorentzian.jl")
include("components/SumReducer.jl")


# ====================================================================
"""
    ModelEval(model::Model, domain::AbstractDomain)

A structure containing the required informations to evaluate a model on a specific domain, and to compare the outcomes to a single empirical dataset.

The model and all component evaluation can be obtained by using the `Model` object has if it was a function: with no arguments it will return the main component evaluation, while if a `Symbol` is given as argument it will return the evaluation of the component with the same name.
"""
struct ModelEvalT{T <: Real}
    pvalues::PVModel{Union{T, Float64}}
    pactual::PVModel{Union{T, Float64}}
    pvmulti::Vector{PVModel{Union{T, Float64}}}

    ModelEvalT{T}() where T  <: Real =
        new(PVModel{Union{T, Float64}}(),
            PVModel{Union{T, Float64}}(),
            Vector{PVModel{Union{T, Float64}}}())
end

function empty!(v::ModelEvalT)
    empty!(v.pvalues)
    empty!(v.pactual)
    empty!(v.pvmulti)
end


struct ModelEval
    model::Model
    domain::AbstractDomain
    cevals::OrderedDict{Symbol, CompEval}
    ifree::Vector{Int}
    patched::Vector{NTuple{2, Symbol}}
    tpar::ModelEvalT{Float64}
    tparad::ModelEvalT{Dual}
    seq::Vector{Symbol}

    function ModelEval(model::Model, domain::AbstractDomain)
        meval = new(model, domain, OrderedDict{Symbol, CompEval}(),
                    Vector{Int}(), Vector{NTuple{2, Symbol}}(),
                    ModelEvalT{Float64}(), ModelEvalT{Dual}(),
                    Vector{CompEval}())
        scan_model!(meval, evaluate=false)
        return meval
    end
end


function scan_model!(meval::ModelEval; evaluate=true)
    function isParamFixed(par::Parameter)
        if !isnothing(par.patch)
            @assert isnothing(par.mpatch) "Parameter [$cname].$pname has both patch and mpatch fields set, while only one is allowed"
            if isa(par.patch, Symbol)  # use same param. value from a different component
                return true
            else                       # invoke a patch function
                @assert length(par.patch.args) in [1,2]
                return  length(par.patch.args) == 1
            end
        elseif !isnothing(par.mpatch)
            @assert length(par.mpatch.args) in [1,2]
            return  length(par.mpatch.args) == 1
        end
        return par.fixed
    end

    empty!(meval.ifree)
    empty!(meval.patched)
    empty!(meval.tpar)
    empty!(meval.tparad)

    isfixed = Vector{Bool}()
    for (cname, comp) in meval.model.comps
        for (pname, par) in getparams(comp)
            if !(par.low <= par.val <= par.high)
                s = "Value outside limits for param [$(cname)].$(pname):\n" * string(par)
                error(s)
            end
            if isnan(par.low)  ||  isnan(par.high)  ||  isnan(par.val)
                s = "NaN value detected for param [$(cname)].$(pname):\n" * string(par)
                error(s)
            end

            push!(isfixed, isParamFixed(par)  ||  meval.model.fixed[cname])
            if !isnothing(par.patch)  ||  !isnothing(par.mpatch)
                push!(meval.patched, (cname, pname))
            end

            push!(meval.tpar.pvalues  , cname, pname, par.val)
            push!(meval.tpar.pactual  , cname, pname, par.val)
            push!(meval.tparad.pvalues, cname, pname, Dual{:tag}(0., 0.))
            push!(meval.tparad.pactual, cname, pname, Dual{:tag}(0., 0.))
        end
        if !(cname in keys(meval.cevals))
            ceval = CompEval(comp, meval.domain)
            meval.cevals[cname] = ceval
        end
    end
    append!(meval.ifree, findall(.! isfixed))

    for (cname, ceval) in meval.cevals
        empty!(ceval.tpar.deps)
        empty!(ceval.tparad.deps)
        i = 1
        for d in dependencies(meval.model, cname, select_domain=true)
            push!(ceval.tpar.deps  , coords(meval.domain, i))
            push!(ceval.tparad.deps, coords(meval.domain, i))
            i += 1
        end
        for d in dependencies(meval.model, cname, select_domain=false)
            push!(ceval.tpar.deps  , meval.cevals[d].tpar.buffer)
            push!(ceval.tparad.deps, meval.cevals[d].tparad.buffer)
        end
    end

    function compeval_sequence!(meval::ModelEval, cname::Symbol)
        for d in dependencies(meval.model, cname)
            compeval_sequence!(meval, d)
        end
        push!(meval.seq, cname)
    end
    function compeval_sequence!(meval::ModelEval)
        empty!(meval.seq)
        compeval_sequence!(meval, find_maincomp(meval.model))
    end
    compeval_sequence!(meval)

    evaluate  &&  GModelFit.evaluate(meval)
    nothing
end


function free_params(meval::ModelEval)
    out = Vector{Parameter}()
    for (cname, comp) in meval.model.comps
        for (pname, par) in getparams(comp)
            push!(out, par)
        end
    end
    return out[meval.ifree]
end
free_params_val(meval::ModelEval) = getfield.(free_params(meval), :val)
nfree(meval::ModelEval) = length(meval.ifree)


# Set new model parameters
function set_pvalues!(meval::ModelEval, pvalues::AbstractVector{Float64})
    items(meval.tpar.pvalues)[meval.ifree] .= pvalues
    items(meval.tpar.pactual)[meval.ifree] .= pvalues
end

function set_pvalues!(meval::ModelEval, pvalues::AbstractVector)
    items(meval.tparad.pvalues)[meval.ifree] .= pvalues
    items(meval.tparad.pactual)[meval.ifree] .= pvalues
end


function run_patch_functs!(meval::ModelEval, tpar::ModelEvalT)
    for (cname, pname) in meval.patched
        par = getfield(meval.model[cname], pname)
        if !isnothing(par.patch)
            @assert isnothing(par.mpatch) "Parameter [:$(cname)].$pname has both patch and mpatch fields set, while only one is allowed"
            if isa(par.patch, Symbol)  # use same param. value from a different component
                tpar.pactual[cname][pname] = tpar.pvalues[par.patch][pname]
            else                       # invoke a patch function
                if length(par.patch.args) == 1
                    tpar.pactual[cname][pname] = par.patch(tpar.pvalues)
                else
                    tpar.pactual[cname][pname] = par.patch(tpar.pvalues, tpar.pvalues[cname][pname])
                end
            end
        elseif !isnothing(par.mpatch)
            @assert length(tpar.pvmulti) > 0 "Parameter [:$(cname)].$pname has the mpatch field set but no other Model is being considered"
            if length(par.mpatch.args) == 1
                tpar.pactual[cname][pname] = par.mpatch(tpar.pvmulti)
            else
                tpar.pactual[cname][pname] = par.mpatch(tpar.pvmulti, tpar.pvalues[cname][pname])
            end
        end
    end
end


"""
    evaluate(meval::ModelEval)

Update a `ModelEval` structure by evaluating all components in the model.
"""
evaluate(meval::ModelEval) = evaluate(meval, free_params_val(meval))

function evaluate(meval::ModelEval, cname::Symbol)
    evaluate(meval)
    return meval.cevals[cname].tpar.buffer
end

function evaluate(meval::ModelEval, pvalues::AbstractVector{Float64})
    set_pvalues!(meval, pvalues)
    run_patch_functs!(meval, meval.tpar)
    for cname in meval.seq
        evaluate(meval.cevals[cname], items(meval.tpar.pactual[cname]))
    end
    return meval.cevals[meval.seq[end]].tpar.buffer
end

function evaluate(meval::ModelEval, pvalues::AbstractVector)
    set_pvalues!(meval, pvalues)
    run_patch_functs!(meval, meval.tparad)
    for cname in meval.seq
        evaluate(meval.cevals[cname], items(meval.tparad.pactual[cname]))
    end
    return meval.cevals[meval.seq[end]].tparad.buffer
end


"""
    evalcounter(meval::ModelEval, cname::Symbol)

Return the number of times a component has been evaluated.
"""
evalcounter(meval::ModelEval, cname::Symbol) = meval.cevals[cname].tpar.counter + meval.cevals[cname].tparad.counter
evalcounter(model::Model, cname::Symbol) = "???"


"""
    evalcounters(meval::ModelEval)

Return a `OrderedDict{Symbol, Int}` with the number of times each model component has been evaluated.
"""
evalcounters(meval::ModelEval) = OrderedDict([cname => evalcounter(meval, cname) for cname in keys(meval.cevals)])


# ====================================================================
# Evaluate Model on the given domain
function (model::Model)(domain::AbstractDomain, cname::Union{Nothing, Symbol}=nothing)
    meval = ModelEval(model, domain)
    isnothing(cname)  &&  (cname = meval.seq[end])
    return meval.cevals[cname].tpar.buffer
end
