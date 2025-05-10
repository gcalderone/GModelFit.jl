
"""
    CompEval(comp::AbstractComponent, domain::AbstractDomain)

A container for a component to be evaluated on a specific domain.

# Fields:
 - `comp::AbstractComponent`: the wrapped component;
 - `domain::AbstractDomain`: the domain where the component is supposed to be evaluated;
 - `counter::Int`: the number of times the component has been evaluated since creatio of the `CompEval` object;
 - `lastparvalues::Vector{Float64}`: the parameter values used in the last evaluation.  A call to `evaluate!()` with the same values stored in `lastparvalues` will not result in a new evaluation;
 - `deps::Vector{Vector{Float64}}`: the evaluation buffers of all dependencies;
 - `buffer::Vector{Float64}`: the buffer to store the outcome of the component.
"""
mutable struct CompEval{TComp <: AbstractComponent, TDomain <: AbstractDomain, T}
    comp::TComp
    domain::TDomain
    counter::Int
    lastparvalues::Vector{Float64}
    deps::Vector{Vector{Float64}}
    buffer::Vector{T}

    function CompEval(comp::AbstractComponent, domain::AbstractDomain, T)
        prepare!(comp, domain)
        return new{typeof(comp), typeof(domain), T}(
            comp, domain, 0,
            fill(NaN, length(getparams(comp))),
            Vector{Vector{Float64}}(),
            Vector{T}(undef, result_length(comp, domain)))
    end
end


"""
    evaluate!(ceval::CompEval, pvalues::Vector{Float64})

Evaluate a component using the provided parameter values.  Outcomes shall be stored in the `CompEval.buffer` vector.
"""
evaluate!(::TComp, ::TDomain, args...) where {TComp <: AbstractComponent, TDomain <: AbstractDomain} =
    error("No evaluate! method implemented for $(TComp), $(TDomain)")


"""
    evaluate_comp!(ceval::CompEval, pvalues::Vector{Float64})

Update a `CompEval` structure using the provided parameter values.

The component is actually evaluated if one of the following applies:
- the component has never been evaluated;
- the component has at least one dependency (whose evaluation may have changed since its last evaluation);
- at least one parameter value has changed since last evaluation.

If none of the above applies, no evaluation occurs.
"""
function evaluate_comp!(ceval::CompEval, pvalues::AbstractVector{T}) where T
    # This is named evaluate_comp! rather than evaluate! to distinguish it from evaluate!(ceval::GModelFit.CompEval{GModelFit.FComp}, params...)
    if length(ceval.deps) > 0
        evaluate!(ceval.comp, ceval.domain, ceval.buffer, ceval.deps, pvalues...)
        ceval.counter += 1
    elseif any(ceval.lastparvalues .!= pvalues)  ||  (ceval.counter == 0)
        evaluate!(ceval.comp, ceval.domain, ceval.buffer, pvalues...)
        (T == Float64)  &&  ceval.lastparvalues .= pvalues
        ceval.counter += 1
    end
    return ceval.buffer
end


# Evaluate component on the given domain.  Parameter values are the
# ones stored in the component unless a custom value is provided via a
# keyword.
function (comp::AbstractComponent)(domain::AbstractDomain; kws...)
    @assert length(dependencies(comp)) == 0 "Can't evaluate a component with dependencies as a stand-alone one."
    ceval = CompEval(comp, domain, Float64)
    pvalues = OrderedDict([(pname, par.val) for (pname, par) in getparams(comp)])
    for (pname, pval) in kws
        if pname in keys(pvalues)
            pvalues[pname] = pval
        else
            @warn "$pname is not a parameter name for $(typeof(comp)). Valid names are: " * join(string.(keys(pvalues)), ", ")
        end
    end
    evaluate_comp!(ceval, collect(values(pvalues)))
    return ceval.buffer
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
struct ModelEval
    model::Model
    domain::AbstractDomain
    cevals::OrderedDict{Symbol, CompEval}
    maincomp::Symbol
    pvalues::PVModel{Float64}
    pactual::PVModel{Float64}
    patched::Vector{NTuple{2, Symbol}}
    ifree::Vector{Int}
    pvmulti::Vector{PVModel{Float64}}
    seq::Vector{Symbol}

    function ModelEval(model::Model, domain::AbstractDomain)
        meval = new(model, domain, OrderedDict{Symbol, CompEval}(), find_maincomp(model),
                    PVModel{Float64}(), PVModel{Float64}(), Vector{NTuple{2, Symbol}}(),
                    Vector{Int}(), Vector{PVModel{Float64}}(), Vector{Symbol}())
        update!(meval, evaluate=false)
        return meval
    end
end


function update!(meval::ModelEval; evaluate=true)
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

    @assert meval.maincomp == find_maincomp(meval.model) "Main component is not allowed to change after ModelEval creation"
    empty!(meval.pvalues)
    empty!(meval.pactual)
    empty!(meval.patched)
    empty!(meval.ifree)
    empty!(meval.pvmulti)

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

            push!(meval.pvalues, cname, pname, par.val)
            push!(meval.pactual, cname, pname, par.val)
            push!(isfixed, isParamFixed(par)  ||  meval.model.fixed[cname])
            if !isnothing(par.patch)  ||  !isnothing(par.mpatch)
                push!(meval.patched, (cname, pname))
            end
        end
        if !(cname in keys(meval.cevals))
            ceval = CompEval(comp, meval.domain, Float64)
            meval.cevals[cname] = ceval
        end
    end
    append!(meval.ifree, findall(.! isfixed))

    for (cname, ceval) in meval.cevals
        empty!(ceval.deps)
        i = 1
        for d in dependencies(meval.model, cname, select_domain=true)
            push!(ceval.deps, coords(meval.domain, i))
            i += 1
        end
        for d in dependencies(meval.model, cname, select_domain=false)
            push!(ceval.deps, meval.cevals[d].buffer)
        end
    end

    function compeval_sequence!(seq::Vector{Symbol}, meval::ModelEval, cname::Symbol)
        for d in dependencies(meval.model, cname)
            compeval_sequence!(seq, meval, d)
        end
        push!(seq, cname)
        return seq
    end
    function compeval_sequence!(meval::ModelEval)
        empty!(meval.seq)
        compeval_sequence!(meval.seq, meval, meval.maincomp)
    end
    compeval_sequence!(meval)

    evaluate  &&  evaluate!(meval)
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
nfree(meval::ModelEval) = length(meval.ifree)


# Set new model parameters
function set_pvalues!(meval::ModelEval, pvalues::Vector{Float64})
    items(meval.pvalues)[meval.ifree] .= pvalues
    items(meval.pactual)[meval.ifree] .= pvalues
end


"""
    evaluate!(meval::ModelEval)

Update a `ModelEval` structure by evaluating all components in the model.
"""
function evaluate!(meval::ModelEval)
    # Patch parameter values
    for (cname, pname) in meval.patched
        par = getfield(meval.model[cname], pname)
        if !isnothing(par.patch)
            @assert isnothing(par.mpatch) "Parameter [:$(cname)].$pname has both patch and mpatch fields set, while only one is allowed"
            if isa(par.patch, Symbol)  # use same param. value from a different component
                meval.pactual[cname][pname] = meval.pvalues[par.patch][pname]
            else                       # invoke a patch function
                if length(par.patch.args) == 1
                    meval.pactual[cname][pname] = par.patch(meval.pvalues)
                else
                    meval.pactual[cname][pname] = par.patch(meval.pvalues, meval.pvalues[cname][pname])
                end
            end
        elseif !isnothing(par.mpatch)
            @assert length(meval.pvmulti) > 0 "Parameter [:$(cname)].$pname has the mpatch field set but no other Model is being considered"
            if length(par.mpatch.args) == 1
                meval.pactual[cname][pname] = par.mpatch(meval.pvmulti)
            else
                meval.pactual[cname][pname] = par.mpatch(meval.pvmulti, meval.pvalues[cname][pname])
            end
        end
    end

    # Evaluation of all components
    for cname in meval.seq
        evaluate_comp!(meval.cevals[cname], items(meval.pactual[cname]))
    end
    return meval
end


"""
    evalcounter(meval::ModelEval, cname::Symbol)

Return the number of times a component has been evaluated.
"""
evalcounter(meval::ModelEval, cname::Symbol) = meval.cevals[cname].counter
evalcounter(model::Model, cname::Symbol) = "???"


"""
    evalcounters(meval::ModelEval)

Return a `OrderedDict{Symbol, Int}` with the number of times each model component has been evaluated.
"""
evalcounters(meval::ModelEval) = OrderedDict([cname => evalcounter(meval, cname) for cname in keys(meval.cevals)])


"""
    last_evaluation(meval::ModelEval)
    last_evaluation(meval::ModelEval, name::Symbol)

Return last evaluation of a component whose name is `cname` in a `ModelEval` object.  If `cname` is not provided the evaluation of the main component is returned.
"""
last_evaluation(meval::ModelEval) = last_evaluation(meval, meval.maincomp)
last_evaluation(meval::ModelEval, name::Symbol) = reshape(meval.domain, meval.cevals[name].buffer)


# ====================================================================
# Evaluate Model on the given domain
function (model::Model)(domain::AbstractDomain, cname::Union{Nothing, Symbol}=nothing)
    meval = ModelEval(model, domain)
    evaluate!(meval)
    if isnothing(cname)
        return last_evaluation(meval)
    else
        return last_evaluation(meval, cname)
    end
end
