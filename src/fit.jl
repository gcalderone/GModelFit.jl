
# ====================================================================
struct FitProblem{T <: AbstractMeasures} <: AbstractFitProblem
    timestamp::DateTime
    meval::ModelEval
    measures::AbstractMeasures
    resid::Vector{Float64}
    nfree::Int
    dof::Int

    function FitProblem(meval::ModelEval, data::T) where T <: AbstractMeasures
        update!(meval)
        resid = fill(NaN, length(data))
        nfree = length(free_params(meval))
        @assert nfree > 0 "No free parameter in the model"
        return new{T}(now(), meval, data, resid, nfree, length(resid) - nfree)
    end
end

free_params(fp::FitProblem) = free_params(fp.meval)

residuals(fp::FitProblem) = fp.resid
function residuals(fp::FitProblem, pvalues::Vector{Float64})
    update_setparvals(fp.meval, pvalues)
    update_evaluation!(fp.meval)
    fp.resid .= reshape((fp.meval() .- values(fp.measures)) ./ uncerts(fp.measures), :)
    return fp.resid
end

fit_stat(fp::FitProblem{Measures{N}}) where N =
    sum(abs2, fp.resid) / fp.dof

function finalize!(fp::FitProblem, best::Vector{Float64}, uncerts::Vector{Float64})
    @assert fp.nfree == length(best) == length(uncerts)
    residuals(fp, best)
    update_finalize!(fp.meval, uncerts)
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

Note: the `FitStats` fields are supposed to be accessed directly by the user, without invoking any get() method.
"""
struct FitStats
    timestamp::DateTime
    elapsed::Float64
    ndata::Int
    nfree::Int
    dof::Int
    fitstat::Float64
    # gofstat::Float64
    # log10testprob::Float64
    status::MinimizerStatus
end

function FitStats(fp::AbstractFitProblem, status::MinimizerStatus)
    # gof_stat = sum(abs2, residuals(fp))
    # tp = logccdf(Chisq(fp.dof), gof_stat) * log10(exp(1))
    FitStats(fp.timestamp, (now() - fp.timestamp).value / 1e3,
             length(residuals(fp)), fp.nfree, fp.dof, fit_stat(fp), # tp,
             status)
end



# ====================================================================
"""
    fit(model::Model, data::Measures; minimizer::AbstractMinimizer=lsqfit())

Fit a model to an empirical data set using the specified minimizer (default: `lsqfit()`).
"""
function fit!(meval::ModelEval, data::Measures; minimizer::AbstractMinimizer=lsqfit())
    update!(meval)
    fp = FitProblem(meval, data)
    status = fit(minimizer, fp)
    return ModelSnapshot(fp.meval), FitStats(fp, status)
end
fit!(model::Model, data::Measures; kws...) = fit!(ModelEval(model, data.domain), data; kws...)
fit( model::Model, data::Measures; kws...) = fit!(deepcopy(model), data; kws...)


"""
    compare(model::Model, data::Measures)

Compare a model to a dataset and return a `FitStats` object.
"""
function compare(meval::ModelEval, data::Measures)
    fp = FitProblem(meval, data)
    status = fit(dry(), fp)
    return FitStats(fp, MinimizerStatus(MinDRY))
end
compare(model::Model, data::Measures) = compare(ModelEval(model, data.domain), data)
