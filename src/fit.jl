# ====================================================================
# Fit statistics
abstract type AbstractFitStat end

# ====================================================================
struct FitProblem{M <: AbstractMeasures, FS <: AbstractFitStat}
    meval::MEval
    data::Vector{M}
    bestfit::Vector{PVModel{Parameter}}
    buffer::Vector{Float64}  # local buffer used to calculate fit statistic
    fitstat::FS

    FitProblem(meval::MEval, datasets::Vector{<: AbstractMeasures}) = FitProblem(meval, datasets, default_fitstat(datasets[1]))
    function FitProblem(meval::MEval, datasets::Vector{M}, fitstat::FS) where {M <: AbstractMeasures, FS <: AbstractFitStat}
        @assert length(meval) == length(datasets)
        fp = new{M, FS}(meval, datasets, Vector{PVModel{Parameter}}(), Vector{Float64}(undef, sum(length.(datasets))), fitstat)
        update_eval!(fp, fp.buffer, free_params_val(fp.meval))
        return fp
    end
end

FitProblem(models::Vector{Model}, datasets::Vector{T}, fitstat=default_fitstat(datasets[1])) where T =
    FitProblem(MEval(models, getfield.(datasets, :domain)), datasets, fitstat)


free_params(fitprob::FitProblem) = free_params(fitprob.meval)
nfree(fitprob::FitProblem) = nfree(fitprob.meval)
ndata(fitprob::FitProblem) = sum(length.(fitprob.data))

function set_bestfit!(fitprob::FitProblem, pvalues::Vector{Float64}, puncerts::Vector{Float64})
    update_eval!(fitprob.meval, pvalues)

    empty!(fitprob.bestfit)
    for id in 1:length(fitprob.meval)
        push!(fitprob.bestfit, PVModel{Parameter}())
    end

    for (id, i1, i2) in free_params_indices(fitprob.meval)
        meval = fitprob.meval.v[id]
        i = 1
        for (cname, comp) in meval.model.comps
            for (pname, _par) in getparams(comp)
                par = deepcopy(_par)
                par.val    = meval.tpar.pvalues[cname][pname]
                par.actual = meval.tpar.pactual[cname][pname]
                push!(fitprob.bestfit[id], cname, pname, par)
                if length(fitprob.bestfit[id].data) in meval.ifree
                    par.unc = puncerts[i1:i2][i]
                    par.fixed = false
                    i += 1
                else
                    par.fixed = true
                end
            end
        end
    end
    nothing
end


# ChiSquared fit statistic
struct ChiSquared <: AbstractFitStat; end

default_fitstat(::Measures) = ChiSquared()
dof(fitprob::FitProblem{M, ChiSquared}) where M = ndata(fitprob) - nfree(fitprob)
fitstat(fitprob::FitProblem{M, ChiSquared}) where M = sum(abs2, fitprob.buffer) / dof(fitprob)

function update_eval!(fitprob::FitProblem{M, ChiSquared}, output::Vector{T}, pvalues::Vector{T}) where {M,T}
    evals = update_eval!(fitprob.meval, pvalues)
    i1 = 1
    for i in 1:length(fitprob.meval)
        nn = length(evals[i])
        if nn > 0
            i2 = i1 + nn - 1
            output[i1:i2] .= reshape((reshape(fitprob.meval.v[i].domain, evals[i]) .- values(fitprob.data[i])) ./ uncerts(fitprob.data[i]), :)
            i1 += nn
        end
    end

    (T == Float64)  &&  (fitprob.buffer .= output)  # update local buffer
    return output
end


# ====================================================================
include("solvers.jl")
using .Solvers
export cmpfit, lsqfit  # export symbols from .Solvers


# ====================================================================
function fit(fitprob::FitProblem, solver=Solvers.lsqfit(); kws...)
    @assert nfree(fitprob) > 0 "No free parameter in the model"
    fsumm = Solvers.solve!(fitprob, solver; kws...)
    bestfit = [ModelSnapshot(fitprob.meval.v[i], fitprob.bestfit[i]) for i in 1:length(fitprob.meval)]
    return bestfit, fsumm
end

function fit!(fitprob::FitProblem, args...; kws...)
    bestfit, fsumm = fit(fitprob, args...; kws...) # invoke non-modifying fit() function
    for i in 1:length(fitprob.meval)
        for (cname, comp) in fitprob.meval.v[i].model
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
    bestfit, fsumm = fit!(FitProblem([model], [data]), args...; kws...)
    return bestfit[1], fsumm
end


"""
    fit(multi::Vector{Model}, data::Vector{Measures{N}}, solver=lsqfit())

Fit a multi-model to a set of empirical data sets using the specified solver (default: `lsqfit()`).  See also `fit!`.
"""
fit(models::Vector{Model}, datasets::Vector{<: AbstractMeasures}, args...; kws...) =
    fit(FitProblem(models, datasets), args...; kws...)


"""
    fit!(multi::Vector{Model}, data::Vector{Measures{N}}, solver=lsqfit())

Fit a multi-model to a set of empirical data sets using the specified solver (default: `lsqfit()`).  Upon return the parameter values in the `Model` objects are set to the best fit ones.
"""
fit!(models::Vector{Model}, datasets::Vector{<: AbstractMeasures}, args...; kws...) =
    fit!(FitProblem(models, datasets), args...; kws...)
