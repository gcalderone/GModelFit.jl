function free_params_indices(multi::Vector{ModelEval})
    out = Vector{NTuple{3, Int}}()
    i1 = 1
    for id in 1:length(multi)
        nn = length(multi[id].pv.ifree)
        if nn > 0
            i2 = i1 + nn - 1
            push!(out, (id, i1, i2))
            i1 += nn
        end
    end
    return out
end


function update!(multi::Vector{ModelEval})
    update_init!(multi)
    update_evaluation!(multi)
    update_finalize!(multi)
    return multi
end

function update_init!(multi::Vector{ModelEval})
    # Populate pvmulti fields in all Model structures to notify we are
    # going to perform a multi-model fitting
    for i in 1:length(multi)
        empty!(multi[i].pvmulti)
        for j in 1:length(multi)
            push!(multi[i].pvmulti, multi[j].pv.values)
        end
    end
    update_init!.(multi)
end

function update_setparvals(multi::Vector{ModelEval}, pvalues::Vector{Float64})
    for (id, i1, i2) in free_params_indices(multi)
        update_setparvals(multi[id], pvalues[i1:i2])
    end
end

function update_evaluation!(multi::Vector{ModelEval})
    update_evaluation!.(multi)
end

function update_finalize!(multi::Vector{ModelEval}, uncerts=Vector{Float64}[])
    if length(uncerts) > 0
        for (id, i1, i2) in free_params_indices(multi)
            update_finalize!(multi[id], pvalues[i1:i2])
        end
    else
        update_finalize!.(multi)
    end
end

function free_params(multi::Vector{ModelEval})
    out = Vector{Parameter}()
    for id in 1:length(multi)
        append!(out, free_params(multi[id]))
    end
    return out
end


# ====================================================================
struct MultiResiduals{T <: AbstractMeasures, M <: AbstractMinimizer} <: AbstractResiduals{T, M}
    multi::Vector{ModelEval}
    resid::Vector{Residuals{T}}
    buffer::Vector{Float64}
    nfree::Int
    dof::Int
    mzer::M

    function MultiResiduals(multi::Vector{ModelEval}, datasets::Vector{T}, mzer::M=dry()) where {T <: AbstractMeasures, M <: AbstractMinimizer}
        @assert length(multi) == length(datasets)
        update!(multi)
        mresid = [Residuals(multi[id], datasets[id]) for id in 1:length(multi)]
        buffer = fill(NaN, sum(length.(getfield.(mresid, :buffer))))
        nfree = sum(getfield.(mresid, :nfree))
        @assert nfree > 0 "No free parameter in the model"
        return new{T,M}(multi, mresid, buffer, nfree, length(buffer) - nfree, mzer)
    end
end

free_params(mresid::MultiResiduals) = free_params(mresid.multi)

residuals(mresid::MultiResiduals) = mresid.buffer
function residuals(mresid::MultiResiduals, pvalues::Vector{Float64})
    # Must set pvalues on all models before any evaluation
    update_setparvals(mresid.multi, pvalues)

    # Populate residuals
    j1 = 1
    for (id, i1, i2) in free_params_indices(mresid.multi)
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
    for (id, i1, i2) in free_params_indices(mresid.multi)
        finalize!(mresid.resid[id], best[i1:i2], uncerts[i1:i2])
    end
end


"""
    fit(multi::Vector{Model}, data::Vector{Measures{N}}; minimizer::AbstractMinimizer=lsqfit())

Fit a multi-model to a set of empirical data sets using the specified minimizer (default: `lsqfit()`).
"""
function fit!(multi::Vector{ModelEval}, data::Vector{Measures{N}}; minimizer::AbstractMinimizer=lsqfit()) where N
    starttime = time()
    update!(multi)
    mresid = MultiResiduals(multi, data, minimizer)
    status = minimize!(mresid)
    bestfit = ModelSnapshot.(mresid.multi)
    stats = FitStats(mresid, status, time() - starttime)
    # test_serialization(bestfit, stats, data)
    return (bestfit, stats)
end
fit!(multi::Vector{Model}, data::Vector{Measures{N}}; kws...) where N = fit!([ModelEval(multi[i], data[i].domain) for i in 1:length(multi)], data; kws...)
fit( multi::Vector{Model}, data::Vector{Measures{N}}; kws...) where N = fit!(deepcopy(multi), data; kws...)


"""
    compare(multi::Vector{Model}, data::Vector{Measures{N}})

Compare a multi-model to a multi-dataset and return a `FitStats` object.
"""
function compare(multi::Vector{ModelEval}, data::Vector{Measures{N}}) where N
    mresid = MultiResiduals(multi, data, dry())
    status = minimize!(mresid)
    return FitStats(mresid, status)
end
compare(multi::Vector{Model}, data::Vector{Measures{N}}) where N = compare([ModelEval(multi[i], data[i].domain) for i in 1:length(multi)], data)
