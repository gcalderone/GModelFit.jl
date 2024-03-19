
# ====================================================================
struct Residuals{T <: AbstractMeasures, M <: AbstractMinimizer} <: AbstractResiduals{T, M}
    meval::ModelEval
    data::T
    buffer::Vector{Float64}
    mzer::M

    function Residuals(meval::ModelEval, data::T, mzer::M=dry()) where {T <: AbstractMeasures, M <: AbstractMinimizer}
        update!(meval)
        buffer = fill(NaN, length(last_evaluation(meval)))
        return new{T,M}(meval, data, buffer, mzer)
    end
end

free_params(resid::Residuals) = free_params(resid.meval)
nfree(resid::Residuals) = nfree(resid.meval)
dof(resid::Residuals) = length(resid.data) - nfree(resid)

residuals(resid::Residuals) = resid.buffer
function residuals(resid::Residuals, pvalues::Vector{Float64})
    update_setparvals(resid.meval, pvalues)
    update_evaluation!(resid.meval)
    resid.buffer .= reshape((last_evaluation(resid.meval) .- values(resid.data)) ./ uncerts(resid.data), :)
    return resid.buffer
end

fit_stat(resid::Residuals{Measures{N}}) where N =
    sum(abs2, resid.buffer) / dof(resid)

function finalize!(resid::Residuals, best::Vector{Float64}, uncerts::Vector{Float64})
    @assert nfree(resid) == length(best) == length(uncerts)
    residuals(resid, best)
    update_finalize!(resid.meval, uncerts)
end


# ====================================================================
"""
    FitStats

A structure representing the results of a fitting process.

# Fields:
- `timestamp::DateTime`: time at which the fitting process has started;
- `elapsed::Float64`: elapsed time (in seconds);
- `ndata::Int`: number of data empirical points;
- `nfree::Int`: number of free parameters;
- `dof::Int`: ndata - nfree;
- `fitstat::Float64`: fit statistics (equivalent ro reduced χ^2 for `Measures` objects);
- `status`: minimizer exit status (tells whether convergence criterion has been satisfied, or if an error has occurred during fitting);

Note: the `FitStats` fields are supposed to be accessed directly by the user.
"""
struct FitStats
    elapsed::Float64
    ndata::Int
    nfree::Int
    dof::Int
    fitstat::Float64
    # gofstat::Float64
    # log10testprob::Float64
    status::AbstractMinimizerStatus
end

function FitStats(resid::AbstractResiduals, status::AbstractMinimizerStatus, elapsed=NaN)
    # gof_stat = sum(abs2, residuals(resid))
    # tp = logccdf(Chisq(resid.dof), gof_stat) * log10(exp(1))
    ndata = length(residuals(resid))
    nf = nfree(resid)
    FitStats(elapsed,
             ndata, nf, ndata - nf, fit_stat(resid), # tp,
             status)
end



# ====================================================================
function fit!(resid::Residuals)
    starttime = time()
    update!(resid.meval)
    @assert nfree(resid) > 0 "No free parameter in the model"
    status = minimize!(resid)
    bestfit = ModelSnapshot(resid.meval)
    stats = FitStats(resid, status, time() - starttime)
    # test_serialization(bestfit, stats, data)
    return (bestfit, stats)
end


"""
    fit!(model::Model, data::Measures; minimizer::AbstractMinimizer=lsqfit())

Fit a model to an empirical data set using the specified minimizer (default: `lsqfit()`).  Upon return the parameter values in the `Model` object are set to the best fit ones.
"""
function fit!(model::Model, data::Measures; minimizer::AbstractMinimizer=lsqfit())
    meval = ModelEval(model, data.domain)
    update!(meval)
    resid = Residuals(meval, data, minimizer)
    return fit!(resid)
end


"""
    fit(model::Model, data::Measures; minimizer::AbstractMinimizer=lsqfit())

Fit a model to an empirical data set using the specified minimizer (default: `lsqfit()`).  See also `fit!`.
"""
fit(model::Model, data::Measures; kws...) = fit!(deepcopy(model), data; kws...)


"""
    compare(model::Model, data::Measures)

Compare a model to a dataset and return a `FitStats` object.
"""
compare(model::Model, data::Measures) = fit!(model, data, minimizer=dry())
