# ====================================================================
"""
    FitProblemChiSq{T <: AbstractMeasures, M <: AbstractMinimizer}

A structure representing the distance between a `ModelEval` and a dataset. The "distance" is expressed in terms of weighted residuals.

A minimizer can be invoked via the `minimize!` function  to reduce such distance by varying the model parameter values.

# Fields:
- `mevals::Vector{ModelEval}`: Model evaluations on the given domains;
- `data::Vector{AbstractMeasures}`: Empirical datasets to be compared to the models;
- `buffer::Vector{Float64}`: Weighted residuals for each point in the domain;
"""
struct FitProblemChiSq <: AbstractFitProblem
    mevals::Vector{ModelEval}
    data::Vector{<: AbstractMeasures}
    buffer::Vector{Float64}

    FitProblemChiSq(model::Model, data::AbstractMeasures) = FitProblemChiSq([model], [data])

    function FitProblemChiSq(models::Vector{Model}, datasets::Vector{<: AbstractMeasures})
        @assert length(models) == length(datasets)
        mevals = ModelEval.(models, getfield.(datasets, :domain))
        update!(mevals)
        @assert nfree(mevals) > 0 "No free parameter in the model"
        buffer = fill(NaN, sum(length.(datasets)))
        return new(mevals, datasets, buffer)
    end
end

free_params(fitprob::FitProblemChiSq) = free_params(fitprob.mevals)
nfree(fitprob::FitProblemChiSq) = nfree(fitprob.mevals)
dof(fitprob::FitProblemChiSq) = sum(length.(fitprob.data)) - nfree(fitprob)
residuals(fitprob::FitProblemChiSq) = return fitprob.buffer
fitstat(fitprob::FitProblemChiSq) = sum(abs2, fitprob.buffer) / dof(fitprob)

function update!(fitprob::FitProblemChiSq, pvalues::Vector{Float64}, puncerts=Float64[])
    set_pvalues!(fitprob.mevals, pvalues)
    update!(fitprob.mevals)

    # Populate residuals
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

    if length(puncerts) > 0
        set_bestfit!(fitprob.mevals, pvalues, puncerts)
    end
    return fitprob.buffer
end


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

    function FitStats(fitprob::T, status::AbstractMinimizerStatus, elapsed::Float64) where T <: AbstractFitProblem
        ndata = length(residuals(fitprob))
        nf = nfree(fitprob)
        return new(elapsed,
                   ndata, nf, ndata - nf,
                   fitstat(fitprob), status)
    end
end


# ====================================================================
function fit(fitprob::AbstractFitProblem, mzer::AbstractMinimizer=lsqfit())
    starttime = time()
    status = minimize!(fitprob, mzer)
    bestfit = [ModelSnapshot(meval) for meval in fitprob.mevals]
    stats = FitStats(fitprob, status, time() - starttime)
    return bestfit, stats
end

function fit!(fitprob::AbstractFitProblem, args...; kws...)
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
    fit(FitProblemChiSq(models, datasets)               , args...; kws...)


"""
    fit!(multi::Vector{Model}, data::Vector{Measures{N}}, minimizer::AbstractMinimizer=lsqfit())

Fit a multi-model to a set of empirical data sets using the specified minimizer (default: `lsqfit()`).  Upon return the parameter values in the `Model` objects are set to the best fit ones.
"""
fit!(models::Vector{Model}, datasets::Vector{Measures{N}}, args...; kws...) where N =
    fit!(FitProblemChiSq(models, datasets)               , args...; kws...)
