# ====================================================================
"""
    FitProblemChiSq{T <: AbstractMeasures, M <: AbstractMinimizer}

A structure representing the distance between a `ModelEval` and a dataset. The "distance" is expressed in terms of weighted residuals.

A minimizer can be invoked via the `minimize!` function  to reduce such distance by varying the model parameter values.

# Fields:
- `meval::ModelEval`: Model evaluation on a given domain;
- `data::AbstractMeasures`: Empirical dataset to be compared to the model;
- `buffer::Vector{Float64}`: Weighted residuals for each point in the domain;
"""
struct FitProblemChiSq <: AbstractFitProblem
    meval::ModelEval
    data::AbstractMeasures
    buffer::Vector{Float64}

    function FitProblemChiSq(model::Model, data::AbstractMeasures)
        meval = ModelEval(model, data.domain)
        update!(meval)
        @assert nfree(meval) > 0 "No free parameter in the model"
        buffer = fill(NaN, length(data))
        return new(meval, data, buffer)
    end
end

free_params(fitprob::FitProblemChiSq) = free_params(fitprob.meval)
nfree(fitprob::FitProblemChiSq) = nfree(fitprob.meval)
dof(fitprob::FitProblemChiSq) = length(fitprob.data) - nfree(fitprob)
residuals(fitprob::FitProblemChiSq) = return fitprob.buffer
fitstat(fitprob::FitProblemChiSq) = sum(abs2, fitprob.buffer) / dof(fitprob)

function update!(fitprob::FitProblemChiSq, pvalues::Vector{Float64}, puncerts=Float64[])
    set_pvalues!(fitprob.meval, pvalues)
    update!(fitprob.meval)
    fitprob.buffer .= reshape((last_evaluation(fitprob.meval) .- values(fitprob.data)) ./ uncerts(fitprob.data), :)
    if length(puncerts) > 0
        set_bestfit!(fitprob.meval, pvalues, puncerts)
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
        nf = nfree(fitprob.meval)
        return new(elapsed,
                   ndata, nf, ndata - nf,
                   fitstat(fitprob), status)
    end
end


# ====================================================================
"""
    fit!(model::Model, data::Measures; minimizer::AbstractMinimizer=lsqfit())

Fit a model to an empirical data set using the specified minimizer (default: `lsqfit()`).  Upon return the parameter values in the `Model` object are set to the best fit ones.  See also `fit`.
"""
function fit!(model::Model, args...; kws...)
    bestfit, stats = fit(model, args...; kws...)
    for (cname, comp) in model
        for (pname, par) in getparams(comp)
            par.val = bestfit[cname][pname].val
            par.unc = NaN
        end
    end
    return bestfit, stats
end


"""
    fit(model::Model, data::Measures; minimizer::AbstractMinimizer=lsqfit())

Fit a model to an empirical data set using the specified minimizer (default: `lsqfit()`).  See also `fit!`.
"""
function fit(model::Model, data::Measures, mzer::AbstractMinimizer=lsqfit())
    starttime = time()
    fitprob = FitProblemChiSq(model, data)
    status = minimize!(fitprob, mzer)
    bestfit = ModelSnapshot(fitprob.meval)
    stats = FitStats(fitprob, status, time() - starttime)
    return bestfit, stats
end


"""
    compare(model::Model, data::Measures)

Compare a model to a dataset and return a `FitStats` object.
"""
compare(model::Model, data::Measures) = fit!(model, data, dry())
