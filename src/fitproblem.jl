# ====================================================================
# Methods to deal with Vector{ModelEval}
#
function free_params_indices(mevals::Vector{ModelEval})
    out = Vector{NTuple{3, Int}}()
    i1 = 1
    for id in 1:length(mevals)
        nn = length(mevals[id].ifree)
        if nn > 0
            i2 = i1 + nn - 1
            push!(out, (id, i1, i2))
            i1 += nn
        end
    end
    return out
end


scan_model!(mevals::Vector{ModelEval}) = scan_model!.(mevals)


function free_params(mevals::Vector{ModelEval})
    out = Vector{Parameter}()
    for id in 1:length(mevals)
        append!(out, free_params(mevals[id]))
    end
    return out
end
free_params_val(mevals::Vector{ModelEval}) = getfield.(free_params(mevals), :val)
nfree(mevals::Vector{ModelEval}) = sum(nfree.(mevals))


evaluate(mevals::Vector{ModelEval}) = evaluate(mevals, free_params_val(mevals))
function evaluate(mevals::Vector{ModelEval}, pvalues::AbstractVector{T}) where T
    if length(mevals) == 1
        return [evaluate(mevals[1], pvalues)]
    else
        if length(mevals[1].tpar.pvmulti) == 0
            pv1 = [mevals[i].tpar.pvalues for i in 1:length(mevals)]
            for i in 1:length(mevals)
                append!(mevals[i].tpar.pvmulti, pv1)
            end
            pv2 = [mevals[i].tparad.pvalues for i in 1:length(mevals)]
            for i in 1:length(mevals)
                append!(mevals[i].tparad.pvmulti, pv2)
            end
        end
        for (id, i1, i2) in free_params_indices(mevals)
            set_pvalues!(mevals[id], pvalues[i1:i2])
        end
        output = Vector{Vector{T}}()
        for (id, i1, i2) in free_params_indices(mevals)
            push!(output, evaluate(mevals[id], pvalues[i1:i2]))
        end
        return output
    end
end


# ====================================================================
# Fit statistics
abstract type AbstractFitStat end
abstract type ChiSquared <: AbstractFitStat end

# ====================================================================
"""
    FitProblem{T <: AbstractFitStat}

A structure containing a `ModelEval` and a dataset, namly all inputs ton define a fit problem.

A solver can be used to reduce the distance between model and data by varying the model parameter values.

# Fields:
- `mevals::Vector{ModelEval}`: Model evaluations on the given domains;
- `data::Vector{AbstractMeasures}`: Empirical datasets to be compared to the models;
- `buffer::Vector{Float64}`: Weighted residuals for each point in the domain;
- `bestfit::Vector{PVModel{Parameter}}`: Best fit values for the parameters
"""
struct FitProblem{T <: AbstractFitStat}
    mevals::Vector{ModelEval}
    data::Vector{<: AbstractMeasures}
    bestfit::Vector{PVModel{Parameter}}
    buffer::Vector{Float64}

    FitProblem(model::Model, data::Measures{N}) where N = FitProblem([model], [data])
    FitProblem(meval::ModelEval, data::Measures{N}) where N = FitProblem([meval], [data])
    FitProblem(models::Vector{Model}, datasets::Vector{Measures{N}}) where N =
        FitProblem(ModelEval.(models, getfield.(datasets, :domain)), datasets)

    function FitProblem(mevals::Vector{ModelEval}, datasets::Vector{Measures{N}}) where N
        @assert length(mevals) == length(datasets)
        fp = new{ChiSquared}(mevals, datasets, Vector{PVModel{Parameter}}(), Vector{Float64}(undef, sum(length.(datasets))))
        return fp
    end
end

free_params(fitprob::FitProblem) = free_params(fitprob.mevals)
nfree(fitprob::FitProblem) = nfree(fitprob.mevals)
ndata(fitprob::FitProblem) = sum(length.(fitprob.data))

function set_bestfit!(fitprob::FitProblem, pvalues::Vector{Float64}, puncerts::Vector{Float64})
    for (id, i1, i2) in free_params_indices(fitprob.mevals)
        set_pvalues!(fitprob.mevals[id], pvalues[i1:i2])
    end
    evaluate(fitprob.mevals, pvalues)

    empty!(fitprob.bestfit)
    for id in 1:length(fitprob.mevals)
        push!(fitprob.bestfit, PVModel{Parameter}())
    end

    for (id, i1, i2) in free_params_indices(fitprob.mevals)
        meval = fitprob.mevals[id]
        i = 1
        for (cname, comp) in meval.model.comps
            for (pname, _par) in getparams(comp)
                par = deepcopy(_par)
                par.val    = meval.tpar.pvalues[cname][pname]
                par.actual = meval.tpar.pactual[cname][pname]
                push!(fitprob.bestfit[id], cname, pname, par)
                if length(fitprob.bestfit[id].data) in meval.ifree
                    par.unc = puncerts[i1:i2][i]
                    par.fixed = false
                    i += 1
                else
                    par.fixed = true
                end
            end
        end
    end
    nothing
end


# FitProblem{ChiSquared} specific methods
dof(fitprob::FitProblem{ChiSquared}) = ndata(fitprob) - nfree(fitprob)
fitstat(fitprob::FitProblem{ChiSquared}) = sum(abs2, fitprob.buffer) / dof(fitprob)

function evaluate_residuals!(output::Vector{T}, fitprob::FitProblem{ChiSquared}, pvalues::Vector{T}) where T
    evals = evaluate(fitprob.mevals, pvalues)
    i1 = 1
    for i in 1:length(fitprob.mevals)
        nn = length(evals[i])
        if nn > 0
            i2 = i1 + nn - 1
            output[i1:i2] .= reshape((reshape(fitprob.mevals[i].domain, evals[i]) .- values(fitprob.data[i])) ./ uncerts(fitprob.data[i]), :)
            i1 += nn
        end
    end

    (T == Float64)  &&  (fitprob.buffer .= output)
    return output
end
