
"""
    CompEval(comp::AbstractComponent, domain::AbstractDomain)

A container for a component to be evaluated on a specific domain.

# Fields:
 - `comp::AbstractComponent`: the wrapped component;
 - `domain::AbstractDomain`: the domain where the component is supposed to be evaluated;
 - `counter::Int`: the number of times the component has been evaluated since creatio of the `CompEval` object;
 - `lastparvalues::Vector{Float64}`: the parameter values used in the last evaluation.  A call to `update!()` with the same values stored in `lastparvalues` will not result in a new evaluation;
 - `deps::Vector{Vector{Float64}}`: the evaluation buffers of all dependencies;
 - `buffer::Vector{Float64}`: the buffer to store the outcome of the component.
"""
mutable struct CompEval{TComp <: AbstractComponent, TDomain <: AbstractDomain}
    comp::TComp
    domain::TDomain
    counter::Int
    lastparvalues::Vector{Float64}
    deps::Vector{Vector{Float64}}
    buffer::Vector{Float64}

    function CompEval(comp::AbstractComponent, domain::AbstractDomain)
        buffer = prepare!(comp, domain)
        return new{typeof(comp), typeof(domain)}(
            comp, domain, 0,
            fill(NaN, length(getparams(comp))),
            Vector{Vector{Float64}}(),
            buffer)
    end
end


"""
    evaluate!(ceval::CompEval, pvalues::Vector{Float64})

Evaluate a component using the provided parameter values.  Outcomes shall be stored in the `CompEval.buffer` vector.
"""
evaluate!(ceval::CompEval{T, D}, par_values...) where {T <: AbstractComponent, D <: AbstractDomain} =
    error("No evaluate! method implemented for CompEval{$(T), $(D)}")


"""
    update!(ceval::CompEval, pvalues::Vector{Float64})

Update a `CompEval` structure using the provided parameter values.

The component is actually evaluated if one of the following applies:
- the component has never been evaluated;
- the component has at least one dependency (whose evaluation may have changed since its last evaluation);
- at least one parameter value has changed since last evaluation.

If none of the above applies, no evaluation occurs.
"""
function update!(ceval::CompEval{<: AbstractComponent, <: AbstractDomain},
                 pvalues::AbstractVector{Float64})
    if any(ceval.lastparvalues .!= pvalues)  ||  (ceval.counter == 0)  ||  (length(ceval.deps) > 0)
        evaluate!(ceval, pvalues...)
        ceval.lastparvalues .= pvalues
        ceval.counter += 1
    end
end


# Evaluate component on the given domain.  Parameter values are the
# ones stored in the component unless a custom value is provided via a
# keyword.
function (comp::AbstractComponent)(domain::AbstractDomain; kws...)
    @assert length(dependencies(comp)) == 0 "Can't evaluate a stand-alone component with dependencies."
    ceval = CompEval(comp, domain)
    pvalues = OrderedDict([(pname, par.val) for (pname, par) in getparams(comp)])
    for (pname, pval) in kws
        if pname in keys(pvalues)
            pvalues[pname] = pval
        else
            @warn "$pname is not a parameter name for $(typeof(comp)). Valid names are: " * join(string.(keys(pvalues)), ", ")
        end
    end
    update!(ceval, collect(values(pvalues)))
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
    ifree::Vector{Int}
    pvmulti::Vector{PVModel{Float64}}
    bestfit::PVModel{Parameter}

    function ModelEval(model::Model, domain::AbstractDomain)
        function isParamFixed(par::Parameter)
            if !isnothing(par.patch)
                @assert isnothing(par.mpatch) "Parameter [$cname].$pname has both patch and mpatch fields set, while only one is allowed"
                if isa(par.patch, Symbol)  # use same param. value from a different component
                    return true
                else                       # invoke a patch function
                    @assert length(par.patch.args) in [1,2]
                    if length(par.patch.args) == 1
                        return true
                    else
                        return false
                    end
                end
            elseif !isnothing(par.mpatch)
                @assert length(par.mpatch.args) in [1,2]
                if length(par.mpatch.args) == 1
                    return true
                else
                    return false
                end
            end
            return par.fixed
        end

        cevals = OrderedDict{Symbol, CompEval}()
        pvalues = PVModel{Float64}()
        pactual = PVModel{Float64}()
        isfixed = Vector{Bool}()

        for (cname, comp) in model.comps
            ceval = CompEval(comp, domain)
            cevals[cname] = ceval

            for (pname, par) in getparams(comp)
                if !(par.low <= par.val <= par.high)
                    s = "Value outside limits for param [$(cname)].$(pname):\n" * string(par)
                    error(s)
                end
                if isnan(par.low)  ||  isnan(par.high)  ||  isnan(par.val)
                    s = "NaN value detected for param [$(cname)].$(pname):\n" * string(par)
                    error(s)
                end

                push!(pvalues, cname, pname, par.val)
                push!(pactual, cname, pname, par.val)
                push!(isfixed, isParamFixed(par)  ||  model.fixed[cname])
            end

            i = 1
            for d in dependencies(model, cname, select_domain=true)
                push!(ceval.deps, coords(domain, i))
                i += 1
            end
            for d in dependencies(model, cname, select_domain=false)
                push!(ceval.deps, cevals[d].buffer)
            end
        end

        return new(model, domain, cevals, find_maincomp(model),
                   pvalues, pactual, findall(.! isfixed),
                   Vector{PVModel{Float64}}(), PVModel{Parameter}())
    end
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
    if !isa(meval.pvalues, PVModel{Float64})
        meval.pvalues = PVModel{Float64}()
        meval.pactual = PVModel{Float64}()
        for (cname, comp) in meval.model.comps
            for (pname, par) in getparams(comp)
                push!(meval.pvalues, cname, pname, par.val)
                push!(meval.pactual, cname, pname, par.val)
            end
        end
    end
    items(meval.pvalues)[meval.ifree] .= pvalues
    items(meval.pactual)[meval.ifree] .= pvalues
end


"""
    update!(meval::ModelEval)

Update a `ModelEval` structure by evaluating all components in the model.
"""
function update!(meval::ModelEval)
    # Patch parameter values
    for (cname, comp) in meval.model.comps
        for (pname, par) in getparams(comp)
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
    end

    # Evaluation of all components, starting from the main one and
    # following dependencies
    function update_compeval_recursive(meval::ModelEval, cname::Symbol)
        for d in dependencies(meval.model, cname)
            update_compeval_recursive(meval, d)
        end
        update!(meval.cevals[cname], items(meval.pactual[cname]))
    end
    update_compeval_recursive(meval, meval.maincomp)
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


function set_bestfit!(meval::ModelEval, pvalues::Vector{Float64}, uncerts::Vector{Float64})
    set_pvalues!(meval, pvalues)
    update!(meval)

    empty!(meval.bestfit)
    i = 1
    for (cname, comp) in meval.model.comps
        for (pname, _par) in getparams(comp)
            par = deepcopy(_par)
            par.val    = meval.pvalues[cname][pname]
            par.actual = meval.pactual[cname][pname]
            push!(meval.bestfit, cname, pname, par)
            if length(meval.bestfit.data) in meval.ifree
                par.unc = uncerts[i]
                par.fixed = false
                i += 1
            else
                par.fixed = true
            end
        end
    end
    nothing
end


# ====================================================================
# Evaluate Model on the given domain
function (model::Model)(domain::AbstractDomain, cname::Union{Nothing, Symbol}=nothing)
    meval = ModelEval(model, domain)
    update!(meval)
    if isnothing(cname)
        return last_evaluation(meval)
    else
        return last_evaluation(meval, cname)
    end
end
