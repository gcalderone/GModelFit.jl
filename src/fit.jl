# ====================================================================
"""
    FitStats

A structure representing the results of a fitting process.

# Fields:
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
    status::AbstractMinimizerStatus
end

function FitStats(fitprob::FitProblem, status::AbstractMinimizerStatus, elapsed::Float64)
    ndata = length(residuals(fitprob))
    nf = nfree(fitprob)
    return FitStats(elapsed,
                    ndata, nf, ndata - nf,
                    fitstat(fitprob), status)
end

# ====================================================================
function fit(fitprob::FitProblem, mzer::AbstractMinimizer=lsqfit())
    starttime = time()
    status = minimize!(fitprob, mzer)
    bestfit = [ModelSnapshot(meval) for meval in fitprob.mevals]
    stats = FitStats(fitprob, status, time() - starttime)
    return bestfit, stats
end

function fit!(fitprob::FitProblem, args...; kws...)
    bestfit, stats = fit(fitprob, args...; kws...) # invoke non-modifying fit() function
    for i in 1:length(fitprob.mevals)
        for (cname, comp) in fitprob.mevals[i].model
            for (pname, par) in getparams(comp)
                par.val = bestfit[i][cname][pname].val
                par.unc = NaN
            end
        end
    end
    return bestfit, stats
end


"""
    compare(model::Model, data::AbstractMeasures)

Compare a model to a dataset and return a `FitStats` object.
"""
compare(model::Model, data::AbstractMeasures) = fit(model, data, dry())

"""
    compare(models::Vector{Model}, data::Vector{<: AbstractMeasures})

Compare a multi-model to a multi-dataset and return a `FitStats` object.
"""
compare(models::Vector{Model}, data::Vector{<: AbstractMeasures}) = fit(models, data, dry())


"""
    fit(model::Model, data::Measures, minimizer::AbstractMinimizer=lsqfit())

Fit a model to an empirical data set using the specified minimizer (default: `lsqfit()`).  See also `fit!`.
"""
function fit(model::Model, data::AbstractMeasures, args...; kws...)
    bestfit, stats = fit([model], [data], args...; kws...)
    return bestfit[1], stats
end


"""
    fit!(model::Model, data::Measures, minimizer::AbstractMinimizer=lsqfit())

Fit a model to an empirical data set using the specified minimizer (default: `lsqfit()`).  Upon return the parameter values in the `Model` object are set to the best fit ones.  See also `fit`.
"""
function fit!(model::Model, data::AbstractMeasures, args...; kws...)
    bestfit, stats = fit!([model], [data], args...; kws...)
    return bestfit[1], stats
end


"""
    fit(multi::Vector{Model}, data::Vector{Measures{N}}, minimizer::AbstractMinimizer=lsqfit())

Fit a multi-model to a set of empirical data sets using the specified minimizer (default: `lsqfit()`).  See also `fit!`.
"""
fit(models::Vector{Model}, datasets::Vector{Measures{N}}, args...; kws...) where N =
    fit(FitProblem(models, datasets)               , args...; kws...)


"""
    fit!(multi::Vector{Model}, data::Vector{Measures{N}}, minimizer::AbstractMinimizer=lsqfit())

Fit a multi-model to a set of empirical data sets using the specified minimizer (default: `lsqfit()`).  Upon return the parameter values in the `Model` objects are set to the best fit ones.
"""
fit!(models::Vector{Model}, datasets::Vector{Measures{N}}, args...; kws...) where N =
    fit!(FitProblem(models, datasets)               , args...; kws...)
