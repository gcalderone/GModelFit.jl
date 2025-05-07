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
    pvmulti = [mevals[i].pvalues for i in 1:length(mevals)]
    for i in 1:length(mevals)
        mevals[i].pvmulti = pvmulti
    end
end

update!(mevals::Vector{ModelEval}) = update!.(mevals)

function set_bestfit!(mevals::Vector{ModelEval}, pvalues::Vector{Float64}, uncerts::Vector{Float64})
    for (id, i1, i2) in free_params_indices(mevals)
        set_bestfit!(mevals[id], pvalues[i1:i2], uncerts[i1:i2])
    end
end


# ====================================================================
"""
    FitProblemChiSqMulti

A structure representing the distance between a `Vector{ModelEval}` and a corresponding number of datasets. The "distance" is expressed in terms of weighted residuals.

A minimizer can be invoked via the `minimize!` function  to reduce such distance by varying the parameter valuesfor each model.

# Fields:
- `mevals::Vector{ModelEval}`: Vector of model evaluations;
- `data::Vector{<: AbstractMeasures}`: Empirical datasets to be compared to the models;
- `buffer::Vector{Float64}`: Weighted residuals for each point in the domain of each dataset;
- `mzer::AbstractMinimizer`: Minimizer used to reduce the residuals.
"""
struct FitProblemChiSqMulti <: AbstractFitProblem
    mevals::Vector{ModelEval}
    data::Vector{<: AbstractMeasures}
    buffer::Vector{Float64}

    function FitProblemChiSqMulti(mevals::Vector{ModelEval}, datasets::Vector{T}) where {T <: AbstractMeasures}
        @assert length(mevals) == length(datasets)
        update!(mevals)
        buffer = fill(NaN, sum(length.(datasets)))
        return new{T,M}(mevals, datasets, buffer, mzer)
    end
end

free_params(fitprob::FitProblemChiSqMulti) = free_params(fitprob.mevals)
nfree(fitprob::FitProblemChiSqMulti) = nfree(fitprob.mevals)
dof(fitprob::FitProblemChiSqMulti) = sum(length.(fitprob.data)) - nfree(fitprob)

function update!(fitprob::FitProblemChiSqMulti, pvalues::Vector{Float64}, puncerts=Float64[]))
    # Must set pvalues on all models before any evaluation
    set_pvalues!(fitprob.mevals, pvalues)

    # Populate residuals
    i1 = 1
    for i in 1:length(fitprob.mevals)
        meval = fitprob.mevals[i]
        update!(meval)
        nn = length(last_evaluation(meval))
        if nn > 0
            i2 = i1 + nn - 1
            fitprob.buffer[i1:i2] .= (last_evaluation(meval) .- values(fitprob.data[i])) ./ uncerts(fitprob.data[i])
            i1 += nn
        end
    end

    if length(puncerts) > 0
        set_bestfit!(fitprob.meval, pvalues, puncerts)
    end

    return fitprob.buffer
end

fit_stat(fitprob::FitProblemChiSqMulti) =
    sum(abs2, fitprob.buffer) / dof(fitprob)

function finalize!(fitprob::FitProblemChiSqMulti, best::Vector{Float64}, uncerts::Vector{Float64})
    @assert nfree(fitprob) == length(best) == length(uncerts)
    residuals(fitprob, best)
    for (id, i1, i2) in free_params_indices(fitprob.mevals)
        # TODO update_finalize!(fitprob.mevals[id], uncerts[i1:i2])
    end
end


# ====================================================================
"""
    minimize!(fitprob::FitProblemChiSqMulti)

Invoke a minimizer to reduce the residuals between a set of models and a corresponding number of datasets.
"""
function minimize!(fitprob::FitProblemChiSqMulti)
    starttime = time()
    update!(fitprob.mevals)
    @assert nfree(fitprob) > 0 "No free parameter in the model"
    status = _minimize!(fitprob)
    bestfit = ModelSnapshot.(fitprob.mevals)
    stats = FitStats(fitprob, status, time() - starttime)
    # test_serialization(bestfit, stats, fitprob.data)
    return (bestfit, stats)
end


"""
    fit!(multi::Vector{Model}, data::Vector{Measures{N}}; minimizer::AbstractMinimizer=lsqfit())

Fit a multi-model to a set of empirical data sets using the specified minimizer (default: `lsqfit()`).  Upon return the parameter values in the `Model` objects are set to the best fit ones.
"""
function fit!(multi::Vector{Model}, data::Vector{Measures{N}}; minimizer::AbstractMinimizer=lsqfit()) where N
    mevals = [ModelEval(multi[i], data[i].domain) for i in 1:length(multi)]
    update!(mevals)
    fitprob = FitProblemChiSqMulti(mevals, data, minimizer)
    return minimize!(fitprob)
end


"""
    fit(multi::Vector{Model}, data::Vector{Measures{N}}; minimizer::AbstractMinimizer=lsqfit())

Fit a multi-model to a set of empirical data sets using the specified minimizer (default: `lsqfit()`).  See also `fit!`.
"""
fit(multi::Vector{Model}, data::Vector{Measures{N}}; kws...) where N = fit!(deepcopy(multi), data; kws...)


"""
    compare(multi::Vector{Model}, data::Vector{Measures{N}})

Compare a multi-model to a multi-dataset and return a `FitStats` object.
"""
compare(multi::Vector{Model}, data::Vector{Measures{N}}) where N = fit!(multi, data, dry())
