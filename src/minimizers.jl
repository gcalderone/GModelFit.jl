# ====================================================================
# Minimizers
#
abstract type AbstractMinimizer end
abstract type AbstractMinimizerStatus end


# --------------------------------------------------------------------
struct MinimizerStatusUnknown <: AbstractMinimizerStatus
end

struct MinimizerStatusDry <: AbstractMinimizerStatus
end

struct MinimizerStatusOK <: AbstractMinimizerStatus
    internal
end

struct MinimizerStatusWarn <: AbstractMinimizerStatus
    message::String
    internal
end

struct MinimizerStatusError <: AbstractMinimizerStatus
    message::String
    internal
end

function as_string(status::AbstractMinimizerStatus)
    if isa(status, MinimizerStatusOK)
        return (crayon"green", "OK", "")
    elseif isa(status, MinimizerStatusWarn)
        return (crayon"bold yellow", "WARNING", status.message)
    elseif isa(status, MinimizerStatusError)
        return (crayon"bold red", "ERROR", status.message)
    elseif isa(status, MinimizerStatusDry)
        return (crayon"bold red", "DRY", "")
    else
        error("Unsupported type: $(typeof(status))")
    end
end


# --------------------------------------------------------------------
struct dry <: AbstractMinimizer; end
function fit!(minimizer::dry, fp::AbstractFitProblem)
    params = free_params(fp)
    println(params)
    finalize!(fp,
              getfield.(params, :val),
              fill(NaN, length(params)))
    return MinimizerStatusDry()
end


# --------------------------------------------------------------------
import LsqFit

struct lsqfit <: AbstractMinimizer; end

function fit!(minimizer::lsqfit, fp::AbstractFitProblem)
    params = free_params(fp)
    ndata = length(residuals(fp))
    prog = ProgressUnknown("Model (dof=$(fp.dof)) evaluations:", dt=0.5, showspeed=true)
    res = LsqFit.curve_fit((dummy, pvalues) -> begin
                           ProgressMeter.next!(prog; showvalues=() -> [(:fit_stat, fit_stat(fp))])
                           update!(fp, pvalues)
                           end,
                           1.:ndata, fill(0., ndata),
                           getfield.(params, :val),
                           lower=getfield.(params, :low),
                           upper=getfield.(params, :high))
    ProgressMeter.finish!(prog)
    if !res.converged
        error!(fp)
        return MinimizerStatusError("Not converged", res)
    end

    finalize!(fp, getfield.(Ref(res), :param), LsqFit.stderror(res))
    return MinimizerStatusOK(res)
end


# --------------------------------------------------------------------
import CMPFit

#=
NOTE: using custom thresholds for ftol, gtol and xtol may lead to
unexpected behaviour.  E.g. settings ftol = 1.e-6 may lead to a
non-optimal fit with exit status 2 (possibly because in a single
iteration the improvement is particularly small).

The best approach is probably to use default tolerance values and
either increase the maximum allowed number of iterations
(config.maxiter) or set a threshold for relative fit statistic
improvements (Δfitstat_threshold) to be checked after the minimizer
iterated for the maximum allowed number of times.
=#

mutable struct cmpfit <: AbstractMinimizer
    config::CMPFit.Config
    Δfitstat_threshold::Float64
    cmpfit() = new(CMPFit.Config(), NaN)
end

function fit!(minimizer::cmpfit, fp::AbstractFitProblem)
    params = free_params(fp)
    guess = getfield.(params, :val)
    low   = getfield.(params, :low)
    high  = getfield.(params, :high)
    parinfo = CMPFit.Parinfo(length(guess))
    for i in 1:length(guess)
        llow  = isfinite(low[i])   ?  1  :  0
        lhigh = isfinite(high[i])  ?  1  :  0
        parinfo[i].limited = (llow, lhigh)
        parinfo[i].limits  = (low[i], high[i])
    end

    update!(fp, guess)
    last_fitstat = sum(abs2, residuals(fp))
    while true
        prog = ProgressUnknown("Model (dof=$(fp.dof)) evaluations:", dt=0.5, showspeed=true)
        res = CMPFit.cmpfit((pvalues) -> begin
                            ProgressMeter.next!(prog; showvalues=() -> [(:fit_stat, fit_stat(fp))])
                            update!(fp, pvalues)
                            end,
                            guess, parinfo=parinfo, config=minimizer.config)
        ProgressMeter.finish!(prog)

        if res.status <= 0
            error!(fp)
            return MinimizerStatusError("Status = $(res.status)", res)
        end

        if (res.status == 5)  &&  isfinite(minimizer.Δfitstat_threshold)
            Δfitstat = (last_fitstat - res.bestnorm) / last_fitstat
            if Δfitstat > minimizer.Δfitstat_threshold
                println("\nReached max. number of iteration but relative Δfitstat = $(Δfitstat) > $(minimizer.Δfitstat_threshold), continue minimization...\n")
                last_fitstat = res.bestnorm
                guess = getfield.(Ref(res), :param)
                continue
            end
        end

        finalize!(fp,
                  getfield.(Ref(res), :param),
                  getfield.(Ref(res), :perror))

        if res.status == 2
            return MinimizerStatusWarn("CMPFit status = 2 may imply one (or more) guess values are too far from optimum", res)
        end
        return MinimizerStatusOK(res)
    end
end
