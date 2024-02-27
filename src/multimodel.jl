function free_params_indices(multi::Vector{Model})
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


function update!(multi::Vector{Model})
    update_step_init(multi)
    update_step_evaluation(multi)
    update_step_finalize(multi)
    return multi
end

function update_step_init(multi::Vector{Model})
    # Populate pvmulti fields in all Model structures to notify we are
    # going to perform a multi-model fitting
    for i in 1:length(multi)
        empty!(multi[i].pvmulti)
        for j in 1:length(multi)
            push!(multi[i].pvmulti, multi[j].pv.values)
        end
    end
    update_step_init.(multi)
end

function update_step_newpvalues(multi::Vector{Model}, pvalues::Vector{Float64})
    for (id, i1, i2) in free_params_indices(multi)
        update_step_newpvalues(multi[id], pvalues[i1:i2])
    end
end

function update_step_evaluation(multi::Vector{Model})
    update_step_evaluation.(multi)
end

function update_step_finalize(multi::Vector{Model}, uncerts=Vector{Float64}[])
    if length(uncerts) > 0
        for (id, i1, i2) in free_params_indices(multi)
            update_step_finalize(multi[id], pvalues[i1:i2])
        end
    else
        update_step_finalize.(multi)
    end
end

function free_params(multi::Vector{Model})
    out = Vector{Parameter}()
    for id in 1:length(multi)
        append!(out, free_params(multi[id]))
    end
    return out
end


# ====================================================================
struct MultiFitProblem <: AbstractFitProblem
    timestamp::DateTime
    multi::Vector{Model}
    fp::Vector{FitProblem}
    resid::Vector{Float64}
    nfree::Int
    dof::Int

    function MultiFitProblem(multi::Vector{Model}, datasets::Vector{T}) where T <: AbstractMeasures
        @assert length(multi) == length(datasets)
        update!(multi)
        fp = [FitProblem(multi[id], datasets[id]) for id in 1:length(multi)]
        resid = fill(NaN, sum(length.(getfield.(fp, :resid))))
        nfree = sum(getfield.(fp, :nfree))
        @assert nfree > 0 "No free parameter in the model"
        return new(now(), multi, fp, resid, nfree, length(resid) - nfree)
    end
end

free_params(fp::MultiFitProblem) = free_params(fp.multi)

residuals(fp::MultiFitProblem) = fp.resid
function residuals(fp::MultiFitProblem, pvalues::Vector{Float64})
    # Must set pvalues on all models before any evaluation
    update_step_newpvalues(fp.multi, pvalues)

    # Populate residuals
    j1 = 1
    for (id, i1, i2) in free_params_indices(fp.multi)
        residuals(fp.fp[id], pvalues[i1:i2])
        nn = length(fp.fp[id].resid)
        if nn > 0
            j2 = j1 + nn - 1
            fp.resid[j1:j2] .= fp.fp[id].resid
            j1 += nn
        end
    end
    return fp.resid
end

fit_stat(fp::MultiFitProblem) =
    sum(abs2, fp.resid) / fp.dof

function finalize!(fp::MultiFitProblem, best::Vector{Float64}, uncerts::Vector{Float64})
    for (id, i1, i2) in free_params_indices(fp.multi)
        finalize!(fp.fp[id], best[i1:i2], uncerts[i1:i2])
    end
end


"""
    fit(multi::Vector{Model}, data::Vector{Measures{N}}; minimizer::AbstractMinimizer=lsqfit())

Fit a multi-model to a set of empirical data sets using the specified minimizer (default: `lsqfit()`).
"""
function fit(multi::Vector{Model}, data::Vector{Measures{N}}; minimizer::AbstractMinimizer=lsqfit()) where N
    # Ensure Model.maincomp has a fixed value
    for i in 1:length(multi)
        push!(multi[i].maincomp, find_maincomp(multi[i]))
    end
    fp = MultiFitProblem(multi, data)
    status = fit(minimizer, fp)
    for i in 1:length(multi)
        deleteat!(multi[i].maincomp, length(multi[i].maincomp))  # restore original empty value
    end
    return ModelSnapshot.(fp.multi), FitStats(fp, status)
end


"""
    compare(multi::Vector{Model}, data::Vector{Measures{N}})

Compare a multi-model to a multi-dataset and return a `FitStats` object.
"""
function compare(multi::Vector{Model}, data::Vector{Measures{N}}) where N
    fp = MultiFitProblem(multi, data)
    status = fit(dry(), fp)
    return FitStats(fp, MinimizerStatus(MinDRY))
end
