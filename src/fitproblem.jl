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


update!(mevals::Vector{ModelEval}) = update!.(mevals)


function free_params(mevals::Vector{ModelEval})
    out = Vector{Parameter}()
    for id in 1:length(mevals)
        append!(out, free_params(mevals[id]))
    end
    return out
end
nfree(mevals::Vector{ModelEval}) = sum(nfree.(mevals))


function set_pvalues!(mevals::Vector{ModelEval}, pvalues::Vector{Float64})
    for (id, i1, i2) in free_params_indices(mevals)
        set_pvalues!(mevals[id], pvalues[i1:i2])
    end
end


function evaluate!(mevals::Vector{ModelEval})
    if length(mevals[1].pvmulti) == 0
        pvmulti = [mevals[i].pvalues for i in 1:length(mevals)]
        for i in 1:length(mevals)
            append!(mevals[i].pvmulti, pvmulti)
        end
    end
    return evaluate!.(mevals)
end


# ====================================================================
# Fit statistics
abstract type AbstractFitStat end
abstract type ChiSquared <: AbstractFitStat end

# ====================================================================
"""
    FitProblem{T <: AbstractFitStat}

A structure representing the distance between a `ModelEval` and a dataset. The "distance" is expressed in terms of weighted residuals.

A minimizer can be invoked via the `minimize!` function  to reduce such distance by varying the model parameter values.

# Fields:
- `mevals::Vector{ModelEval}`: Model evaluations on the given domains;
- `data::Vector{AbstractMeasures}`: Empirical datasets to be compared to the models;
- `buffer::Vector{Float64}`: Weighted residuals for each point in the domain;
"""
struct FitProblem{T <: AbstractFitStat}
    mevals::Vector{ModelEval}
    data::Vector{<: AbstractMeasures}
    buffer::Vector{Float64}
    bestfit::Vector{PVModel{Parameter}}

    FitProblem(model::Model, data::Measures{N}) where N = FitProblem([model], [data])
    FitProblem(meval::ModelEval, data::Measures{N}) where N = FitProblem([meval], [data])
    FitProblem(models::Vector{Model}, datasets::Vector{Measures{N}}) where N =
        FitProblem(ModelEval.(models, getfield.(datasets, :domain)), datasets)

    function FitProblem(mevals::Vector{ModelEval}, datasets::Vector{Measures{N}}) where N
        @assert length(mevals) == length(datasets)
        evaluate!(mevals)
        buffer = fill(NaN, sum(length.(datasets)))
        fp = new{ChiSquared}(mevals, datasets, buffer, Vector{PVModel{Parameter}}())
        populate_residuals!(fp)
        return fp
    end
end

free_params(fitprob::FitProblem) = free_params(fitprob.mevals)
nfree(fitprob::FitProblem) = nfree(fitprob.mevals)
residuals(fitprob::FitProblem) = return fitprob.buffer

function set_bestfit!(fitprob::FitProblem, pvalues::Vector{Float64}, puncerts::Vector{Float64})
    for (id, i1, i2) in free_params_indices(fitprob.mevals)
        meval = fitprob.mevals[id]
        set_pvalues!(meval, pvalues[i1:i2])
        evaluate!(meval)
    end

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
                par.val    = meval.pvalues[cname][pname]
                par.actual = meval.pactual[cname][pname]
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


function evaluate!(fitprob::FitProblem{ChiSquared}, pvalues::Vector{Float64})
    set_pvalues!(fitprob.mevals, pvalues)
    evaluate!(fitprob.mevals)
    return populate_residuals!(fitprob)
end


# FitProblem{ChiSquared} specific methods
dof(fitprob::FitProblem{ChiSquared}) = sum(length.(fitprob.data)) - nfree(fitprob)
fitstat(fitprob::FitProblem{ChiSquared}) = sum(abs2, fitprob.buffer) / dof(fitprob)
function populate_residuals!(fitprob::FitProblem{ChiSquared})
    i1 = 1
    for i in 1:length(fitprob.mevals)
        meval = fitprob.mevals[i]
        nn = length(last_evaluation(meval))
        if nn > 0
            i2 = i1 + nn - 1
            fitprob.buffer[i1:i2] .= reshape((last_evaluation(meval) .- values(fitprob.data[i])) ./ uncerts(fitprob.data[i]), :)
            i1 += nn
        end
    end
    return fitprob.buffer
end
