
"""
    CompEval(comp::AbstractComponent, domain::AbstractDomain)

A container for a component to be evaluated on a specific domain.

# Fields:
 - `comp::AbstractComponent`: the wrapped component;
 - `domain::AbstractDomain`: the domain where the component is supposed to be evaluated;
 - `counter::Int`: the number of times the component has been evaluated since creatio of the `CompEval` object;
 - `lastparvalues::Vector{Float64}`: the parameter values used in the last evaluation.  A call to `update!()` with the same values stored in `lastparvalues` will not result in a new evaluation;
 - `deps::Vector{Vector{Float64}}`: the buffers of all dependencies;
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
    update!(ceval::CompEval, pvalues::Vector{Float64})

Evaluate a `CompEval` structure with a new component evaluation using the provided parameter values.
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
function (comp::AbstractComponent)(domain::AbstractDomain, deps=Vector{Vector{Float64}}(); kws...)
    ceval = CompEval(comp, domain)
    append!(ceval.deps, deps)
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
struct ParameterVectors
    params::PVModel{Parameter}
    values::PVModel{Float64}
    actual::PVModel{Float64}
    ifree::Vector{Int}
    ParameterVectors() = new(PVModel{Parameter}(), PVModel{Float64}(),
                             PVModel{Float64}(), Vector{Int}())
end
function empty!(pv::ParameterVectors)
    empty!(pv.params)
    empty!(pv.values)
    empty!(pv.actual)
    empty!(pv.ifree)
end
function push!(pv::ParameterVectors, cname::Symbol, pname::Symbol, par::Parameter)
    push!(pv.params, cname, pname, par)
    push!(pv.values, cname, pname, par.val)
    push!(pv.actual, cname, pname, par.actual)
    if !par.fixed
        push!(pv.ifree, length(items(pv.params)))
    end
end


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
    pv::ParameterVectors
    pvmulti::Vector{PVModel{Float64}}
    maincomp::Symbol

    function ModelEval(model::Model, domain::AbstractDomain)
        out = new(model, domain, OrderedDict{Symbol, CompEval}(),
                  ParameterVectors(), Vector{PVModel{Float64}}(),
                  find_maincomp(model))
        # update!(out)  This would cause error in the multimodel case
        return out
    end
end


free_params(meval::ModelEval) = collect(items(meval.pv.params)[meval.pv.ifree])


"""
    update!(meval::ModelEval)

Update a `ModelEval` structure by evaluating all components in the model.
"""
function update!(meval::ModelEval)
    update_init!(meval)
    update_evaluation!(meval)
    update_finalize!(meval)
    return meval
end


# Evaluation step init:
# - update internal structures before fitting
function update_init!(meval::ModelEval)
    for (cname, comp) in meval.model.comps
        (cname in keys(meval.cevals))  &&  continue
        meval.cevals[cname] = CompEval(comp, meval.domain)
    end

    empty!(meval.pv)
    for (cname, ceval) in meval.cevals
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
                @assert length(meval.pvmulti) > 0 "Parameter [:$(cname)].$pname has the mpatch field set but no other Model is being considered"
                @assert length(par.mpatch.args) in [1,2]
                if length(par.mpatch.args) == 1
                    par.fixed = true
                else
                    par.fixed = false
                end
            end
            if meval.model.fixed[cname]
                par.fixed = true
            end
            push!(meval.pv, cname, pname, par)
        end

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
end


# Set new model parameters
function update_setparvals(meval::ModelEval, pvalues::Vector{Float64})
    items(meval.pv.values)[meval.pv.ifree] .= pvalues
end


# Evaluation step fit:
# - copy all parameter values into actual;
# - update actual by invoking the patch functions;
# - evaluation of all components
function update_evaluation!(meval::ModelEval)
    # Copy pvalues into actual
    items(meval.pv.actual) .= items(meval.pv.values)
    # Patch parameter values
    for (cname, comp) in meval.pv.params
        for (pname, par) in comp
            if !isnothing(par.patch)
                @assert isnothing(par.mpatch) "Parameter [:$(cname)].$pname has both patch and mpatch fields set, while only one is allowed"
                if isa(par.patch, Symbol)  # use same param. value from a different component
                    meval.pv.actual[cname][pname] = meval.pv.values[par.patch][pname]
                else                       # invoke a patch function
                    if length(par.patch.args) == 1
                        meval.pv.actual[cname][pname] = par.patch(meval.pv.values)
                    else
                        meval.pv.actual[cname][pname] = par.patch(meval.pv.values, meval.pv.values[cname][pname])
                    end
                end
            elseif !isnothing(par.mpatch)
                @assert length(meval.pvmulti) > 0 "Parameter [:$(cname)].$pname has the mpatch field set but no other Model is being considered"
                if length(par.mpatch.args) == 1
                    meval.pv.actual[cname][pname] = par.mpatch(meval.pvmulti)
                else
                    meval.pv.actual[cname][pname] = par.mpatch(meval.pvmulti, meval.pv.values[cname][pname])
                end
            end
        end
    end

    # Evaluation of all components, starting from the main one and
    # following dependencies
    function update_compeval_recursive(meval::ModelEval, cname::Symbol)
        @batch_when_threaded per=core for d in dependencies(meval.model, cname)
            update_compeval_recursive(meval, d)
        end
        update!(meval.cevals[cname], items(meval.pv.actual[cname]))
    end
    update_compeval_recursive(meval, meval.maincomp)
end


# Evaluation step finalize:
# - copy back bestfit, actual values and uncertainties into their original Parameter structures.
function update_finalize!(meval::ModelEval, uncerts=Vector{Float64}[])
    i = 1
    for (cname, comp) in meval.pv.params
        for (pname, par) in comp
            par.val    = meval.pv.values[cname][pname]
            par.actual = meval.pv.actual[ cname][pname]
            if (length(uncerts) > 0)  &&  (!par.fixed)
                par.unc = uncerts[i]
                i += 1
            else
                par.unc = NaN
            end
        end
    end

    # Also update Model's parameters
    for (cname, ceval) in meval.cevals
        setparams!(meval.model[cname], meval.pv.params[cname])
    end
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


# Return model evaluations
(meval::ModelEval)() = meval(meval.maincomp)
(meval::ModelEval)(name::Symbol) = reshape(meval.domain, meval.cevals[name].buffer)


# ====================================================================
# Evaluate Model on the given domain
function (model::Model)(domain::AbstractDomain, cname::Union{Nothing, Symbol}=nothing)
    meval = ModelEval(model, domain)
    update!(meval)
    if isnothing(cname)
        return meval()
    else
        return meval(cname)
    end
end
