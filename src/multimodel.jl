function free_params_indices(mevals::Vector{ModelEval})
    out = Vector{NTuple{3, Int}}()
    i1 = 1
    for id in 1:length(mevals)
        nn = length(mevals[id].pv.ifree)
        if nn > 0
            i2 = i1 + nn - 1
            push!(out, (id, i1, i2))
            i1 += nn
        end
    end
    return out
end


function update!(mevals::Vector{ModelEval})
    update_init!(mevals)
    update_evaluation!(mevals)
    update_finalize!(mevals)
    return mevals
end

function update_init!(mevals::Vector{ModelEval})
    # Populate pvmulti fields in all Model structures to notify we are
    # going to perform a multi-model fitting
    for i in 1:length(mevals)
        empty!(mevals[i].pvmulti)
        for j in 1:length(mevals)
            push!(mevals[i].pvmulti, mevals[j].pv.values)
        end
    end
    update_init!.(mevals)
end

function update_setparvals(mevals::Vector{ModelEval}, pvalues::Vector{Float64})
    for (id, i1, i2) in free_params_indices(mevals)
        update_setparvals(mevals[id], pvalues[i1:i2])
    end
end

function update_evaluation!(mevals::Vector{ModelEval})
    update_evaluation!.(mevals)
end

function update_finalize!(mevals::Vector{ModelEval}, uncerts=Vector{Float64}[])
    if length(uncerts) > 0
        for (id, i1, i2) in free_params_indices(mevals)
            update_finalize!(mevals[id], pvalues[i1:i2])
        end
    else
        update_finalize!.(mevals)
    end
end

function free_params(mevals::Vector{ModelEval})
    out = Vector{Parameter}()
    for id in 1:length(mevals)
        append!(out, free_params(mevals[id]))
    end
    return out
end


# ====================================================================
struct MultiResiduals{T <: AbstractMeasures, M <: AbstractMinimizer} <: AbstractResiduals{T, M}
    mevals::Vector{ModelEval}
    resid::Vector{Residuals{T}}
    buffer::Vector{Float64}
    nfree::Int
    dof::Int
    mzer::M

    function MultiResiduals(mevals::Vector{ModelEval}, datasets::Vector{T}, mzer::M=dry()) where {T <: AbstractMeasures, M <: AbstractMinimizer}
        @assert length(mevals) == length(datasets)
        update!(mevals)
        mresid = [Residuals(mevals[id], datasets[id]) for id in 1:length(mevals)]
        buffer = fill(NaN, sum(length.(getfield.(mresid, :buffer))))
        nfree = sum(getfield.(mresid, :nfree))
        @assert nfree > 0 "No free parameter in the model"
        return new{T,M}(mevals, mresid, buffer, nfree, length(buffer) - nfree, mzer)
    end
end

free_params(mresid::MultiResiduals) = free_params(mresid.mevals)

residuals(mresid::MultiResiduals) = mresid.buffer
function residuals(mresid::MultiResiduals, pvalues::Vector{Float64})
    # Must set pvalues on all models before any evaluation
    update_setparvals(mresid.mevals, pvalues)

    # Populate residuals
    j1 = 1
    for (id, i1, i2) in free_params_indices(mresid.mevals)
        residuals(mresid.resid[id], pvalues[i1:i2])
        nn = length(mresid.resid[id].buffer)
        if nn > 0
            j2 = j1 + nn - 1
            mresid.buffer[j1:j2] .= mresid.resid[id].buffer
            j1 += nn
        end
    end
    return mresid.buffer
end

fit_stat(mresid::MultiResiduals) =
    sum(abs2, mresid.buffer) / mresid.dof

function finalize!(mresid::MultiResiduals, best::Vector{Float64}, uncerts::Vector{Float64})
    for (id, i1, i2) in free_params_indices(mresid.mevals)
        finalize!(mresid.resid[id], best[i1:i2], uncerts[i1:i2])
    end
end


# ====================================================================
function fit!(mresid::MultiResiduals)
    @assert mresid.nfree > 0 "No free parameter in the model"
    starttime = time()
    update!(mresid.mevals)
    status = minimize!(mresid)
    bestfit = ModelSnapshot.(mresid.mevals)
    stats = FitStats(mresid, status, time() - starttime)
    # test_serialization(bestfit, stats, data)
    return (bestfit, stats)
end


"""
    fit!(multi::Vector{Model}, data::Vector{Measures{N}}; minimizer::AbstractMinimizer=lsqfit())

Fit a multi-model to a set of empirical data sets using the specified minimizer (default: `lsqfit()`).  Upon return the parameter values in the `Model` objects are set to the best fit ones.
"""
function fit!(multi::Vector{Model}, data::Vector{Measures{N}}; minimizer::AbstractMinimizer=lsqfit()) where N
    mevals = [ModelEval(multi[i], data[i].domain) for i in 1:length(multi)]
    update!(mevals)
    mresid = MultiResiduals(mevals, data, minimizer)
    return fit!(mresid)
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
