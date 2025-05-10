# ====================================================================
"""
    FitSummary

A structure summarizing the results of a fitting process.

# Fields:
- `elapsed::Float64`: elapsed time (in seconds);
- `ndata::Int`: number of data empirical points;
- `nfree::Int`: number of free parameters;
- `dof::Int`: ndata - nfree;
- `fitstat::Float64`: fit statistics (equivalent ro reduced χ^2 for `Measures` objects);
- `status`: solver exit status (tells whether convergence criterion has been satisfied, or if an error has occurred during fitting);

Note: the `FitSummary` fields are supposed to be accessed directly by the user.
"""
struct FitSummary
    elapsed::Float64
    ndata::Int
    nfree::Int
    dof::Int
    fitstat::Float64
    status::AbstractSolverStatus
end

FitSummary(fitprob::FitProblem, status::AbstractSolverStatus, elapsed::Float64) =
    FitSummary(elapsed, ndata(fitprob), nfree(fitprob), dof(fitprob),
               fitstat(fitprob), status)


# ====================================================================
# function fit(fitprob::FitProblem, solver::Union{AbstractSolver, NonlinearSolveBase.AbstractNonlinearSolveAlgorithm}=lsqfit())
function fit(fitprob::FitProblem, solver=Solvers.lsqfit())
    starttime = time()
    @assert nfree(fitprob) > 0 "No free parameter in the model"
    status = Solvers.solve!(fitprob, solver)
    bestfit = [ModelSnapshot(fitprob.mevals[i], fitprob.bestfit[i]) for i in 1:length(fitprob.mevals)]
    stats = FitSummary(fitprob, status, time() - starttime)
    return bestfit, stats
end

function fit!(fitprob::FitProblem, args...; kws...)
    bestfit, stats = fit(fitprob, args...; kws...) # invoke non-modifying fit() function
    for i in 1:length(fitprob.mevals)
        for (cname, comp) in fitprob.mevals[i].model
            for (pname, par) in getparams(comp)
                par.val = bestfit[i][cname][pname].val
                par.unc = bestfit[i][cname][pname].unc
            end
        end
    end
    return bestfit, stats
end


"""
    fitstat(model::Model, data::AbstractMeasures)

Compare a model to a dataset and return the fit statistic.
"""
fitstat(model::Model, data::AbstractMeasures) = fitstat(FitProblem(model, data))

"""
    fitstat(models::Vector{Model}, data::Vector{<: AbstractMeasures})

Compare a multi-model to a multi-dataset and return fit statistic.
"""
fitstat(models::Vector{Model}, data::Vector{<: AbstractMeasures}) = fitstat(FitProblem(models, data))


"""
    fit(model::Model, data::Measures, solver=lsqfit())

Fit a model to an empirical data set using the specified solver (default: `lsqfit()`).  See also `fit!`.
"""
function fit(model::Model, data::AbstractMeasures, args...; kws...)
    bestfit, stats = fit(FitProblem(model, data), args...; kws...)
    return bestfit[1], stats
end


"""
    fit!(model::Model, data::Measures, solver=lsqfit())

Fit a model to an empirical data set using the specified solver (default: `lsqfit()`).  Upon return the parameter values in the `Model` object are set to the best fit ones.  See also `fit`.
"""
function fit!(model::Model, data::AbstractMeasures, args...; kws...)
    bestfit, stats = fit!(FitProblem(model, data), args...; kws...)
    return bestfit[1], stats
end


"""
    fit(multi::Vector{Model}, data::Vector{Measures{N}}, solver=lsqfit())

Fit a multi-model to a set of empirical data sets using the specified solver (default: `lsqfit()`).  See also `fit!`.
"""
fit(models::Vector{Model}, datasets::Vector{Measures{N}}, args...; kws...) where N =
    fit(FitProblem(models, datasets), args...; kws...)


"""
    fit!(multi::Vector{Model}, data::Vector{Measures{N}}, solver=lsqfit())

Fit a multi-model to a set of empirical data sets using the specified solver (default: `lsqfit()`).  Upon return the parameter values in the `Model` objects are set to the best fit ones.
"""
fit!(models::Vector{Model}, datasets::Vector{Measures{N}}, args...; kws...) where N =
    fit!(FitProblem(models, datasets), args...; kws...)
