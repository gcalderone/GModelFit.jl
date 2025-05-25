# ====================================================================
function fit(fitprob::FitProblem, solver=Solvers.lsqfit(); kws...)
    @assert nfree(fitprob) > 0 "No free parameter in the model"
    fsumm = Solvers.solve!(fitprob, solver; kws...)
    bestfit = [ModelSnapshot(fitprob.mevals[i], fitprob.bestfit[i]) for i in 1:length(fitprob.mevals)]
    return bestfit, fsumm
end

function fit!(fitprob::FitProblem, args...; kws...)
    bestfit, fsumm = fit(fitprob, args...; kws...) # invoke non-modifying fit() function
    for i in 1:length(fitprob.mevals)
        for (cname, comp) in fitprob.mevals[i].model
            for (pname, par) in getparams(comp)
                par.val = bestfit[i][cname][pname].val
                par.unc = bestfit[i][cname][pname].unc
            end
        end
    end
    return bestfit, fsumm
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
    bestfit, fsumm = fit(FitProblem([model], [data]), args...; kws...)
    return bestfit[1], fsumm
end


"""
    fit!(model::Model, data::Measures, solver=lsqfit())

Fit a model to an empirical data set using the specified solver (default: `lsqfit()`).  Upon return the parameter values in the `Model` object are set to the best fit ones.  See also `fit`.
"""
function fit!(model::Model, data::AbstractMeasures, args...; kws...)
    bestfit, fsumm = fit!(FitProblem(model, data), args...; kws...)
    return bestfit[1], fsumm
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
