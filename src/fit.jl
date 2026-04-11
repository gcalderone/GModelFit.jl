# ====================================================================
# Fit statistics
abstract type AbstractFitStat end

# ====================================================================
struct Likelihood{M <: AbstractMeasures, FS <: AbstractFitStat}
    mseval::ModelSetEval
    data::Vector{M}
    buffer::Vector  # local buffer used to calculate fit statistic
    fitstat::FS

    function Likelihood(datasets::Vector{M}, ms::ModelSet, fitstat::AbstractFitStat=default_fitstat(datasets[1]); use_AD=false) where {M <: AbstractMeasures}
        @assert length(ms) == length(datasets)
        T = use_AD  ?  Union{Dual, Float64}  :  Float64
        mseval = ModelSetEval{T}(ms, getfield.(datasets, :domain))

        out = new{M, typeof(fitstat)}(mseval, datasets,
                                      Vector{T}(undef, sum(length.(datasets))),
                                      fitstat)
        out.buffer .= 0.  # needed to avoid errors with CMPFit
        return out
    end
end


getparams(lh::Likelihood; kws...) = getparams(lh.mseval; kws...)
nfree(lh::Likelihood) = nfree(lh.mseval)
ndata(lh::Likelihood) = length(lh.buffer)
set_bestfit!(lh::Likelihood, pvalues::Vector{Float64}, puncerts::Vector{Float64}) = set_bestfit!(lh.mseval, pvalues, puncerts)


# ChiSquared fit statistic ===========================================
struct ChiSquared <: AbstractFitStat; end

default_fitstat(::Measures) = ChiSquared()
dof(lh::Likelihood{M, ChiSquared}) where M = ndata(lh) - nfree(lh)
fitstat(lh::Likelihood{M, ChiSquared}) where M = loglikelihood(lh) / (-0.5) / dof(lh)
loglikelihood(lh::Likelihood{M, ChiSquared}) where M = -0.5 * sum(abs2, lh.buffer)

function update_eval!(lh::Likelihood{M, ChiSquared}, pvalues::Vector{T}) where {M,T}
    update_eval!(lh.mseval, pvalues)
    i1 = 1
    for i in 1:length(lh.mseval.vec)
        model = lh.mseval.vec[i].folded
        nn = length(model)
        if nn > 0
            i2 = i1 + nn - 1
            lh.buffer[i1:i2] .= view((model .- values(lh.data[i])) ./ uncerts(lh.data[i]), :)
            i1 += nn
        end
    end
    return convert(Vector{T}, lh.buffer)
end


# ====================================================================
include("solvers.jl")
using .Solvers


# ====================================================================
"""
    fitstat(model::Model, data::AbstractMeasures)

Compare a model to a dataset and return the fit statistic.
"""
fitstat(model::Model, data::AbstractMeasures) = fitstat(Likelihood(data, model))

"""
    fitstat(models::Vector{Model}, data::Vector{<: AbstractMeasures})

Compare a multi-model to a multi-dataset and return fit statistic.
"""
fitstat(ms::ModelSet, data::Vector{<: AbstractMeasures}) = fitstat(Likelihood(data, ms))


function fit(lh::Likelihood, solver::AbstractSolver)
    @assert nfree(lh) > 0 "No free parameter in the model"
    fsumm = Solvers.solve!(lh, solver)
    bestfit = ModelSetSnapshot(lh.mseval)
    return bestfit, fsumm
end

function fit!(lh::Likelihood, solver::AbstractSolver)
    bestfit, fsumm = fit(lh, solver) # invoke non-modifying fit() function
    newpars = getparams(lh)
    for (key, par) in getparams(lh.mseval.ms)
        par.val = newpars[key].val
        par.unc = newpars[key].unc
    end
    return bestfit, fsumm
end

"""
    fit(multi::Vector{Model}, data::Vector{Measures{N}}, solver=lsqfit())

Fit a multi-model to a set of empirical data sets using the specified solver (default: `lsqfit()`).  See also `fit!`.
"""
fit(ms::ModelSet, datasets::Vector{<: AbstractMeasures}, solver::AbstractSolver=Solvers.lsqfit()) =
    fit(Likelihood(datasets, ms, use_AD=use_AD(solver)), solver)

"""
    fit!(multi::Vector{Model}, data::Vector{Measures{N}}, solver=lsqfit())

Fit a multi-model to a set of empirical data sets using the specified solver (default: `lsqfit()`).  Upon return the parameter values in the `Model` objects are set to the best fit ones.
"""
fit!(ms::ModelSet, datasets::Vector{<: AbstractMeasures}, solver::AbstractSolver=Solvers.lsqfit()) =
    fit!(Likelihood(datasets, ms, use_AD=use_AD(solver)), solver)

"""
    fit(model::Model, data::Measures, solver=lsqfit())

Fit a model to an empirical data set using the specified solver (default: `lsqfit()`).  See also `fit!`.
"""
function fit(model::Model, data::AbstractMeasures, args...)
    bestfit, fsumm = fit(ModelSet(:_ => model), [data], args...)
    return bestfit[:_], fsumm
end

"""
    fit!(model::Model, data::Measures, solver=lsqfit())

Fit a model to an empirical data set using the specified solver (default: `lsqfit()`).  Upon return the parameter values in the `Model` object are set to the best fit ones.  See also `fit`.
"""
function fit!(model::Model, data::AbstractMeasures, args...)
    bestfit, fsumm = fit!(ModelSet(:_ => model), [data], args...)
    return bestfit[:_], fsumm
end
