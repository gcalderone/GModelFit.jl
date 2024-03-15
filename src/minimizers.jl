# ====================================================================
# Minimizers
#

# --------------------------------------------------------------------
abstract type AbstractMinimizerStatus end

struct MinimizerStatusOK <: AbstractMinimizerStatus
end

struct MinimizerStatusDry <: AbstractMinimizerStatus
end

struct MinimizerStatusWarn <: AbstractMinimizerStatus
    message::String
end

struct MinimizerStatusError <: AbstractMinimizerStatus
    message::String
end


# --------------------------------------------------------------------
abstract type AbstractMinimizer end


# --------------------------------------------------------------------
struct dry <: AbstractMinimizer; end
function fit(minimizer::dry, fp::AbstractFitProblem)
    params = free_params(fp)
    residuals(fp, getfield.(params, :val))
    finalize!(fp,
              getfield.(params, :val),
              fill(NaN, length(params)))
    return MinimizerStatusDry
end


# --------------------------------------------------------------------
import LsqFit

mutable struct lsqfit <: AbstractMinimizer;
    result
    lsqfit() = new(nothing)
end

function fit(mzer::lsqfit, fp::AbstractFitProblem)
    params = free_params(fp)
    ndata = length(residuals(fp))
    prog = ProgressUnknown(desc="Model (dof=$(fp.dof)) evaluations:", dt=0.5, showspeed=true, color=:light_black)
    mzer.result = LsqFit.curve_fit((dummy, pvalues) -> begin
                                       ProgressMeter.next!(prog; showvalues=() -> [(:fit_stat, fit_stat(fp))])
                                       residuals(fp, pvalues)
                                   end,
                                   1.:ndata, fill(0., ndata),
                                   getfield.(params, :val),
                                   lower=getfield.(params, :low),
                                   upper=getfield.(params, :high))
    ProgressMeter.finish!(prog)
    if !mzer.result.converged
        return MnimizerStatusError("Not converged")
    end

    finalize!(fp, getfield.(Ref(mzer.result), :param), LsqFit.stderror(mzer.result))
    return MinimizerStatusOK()
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
improvements (ftol_after_maxiter) to be checked after the minimizer
iterated for the maximum allowed number of times.
=#

mutable struct cmpfit <: AbstractMinimizer
    config::CMPFit.Config
    ftol_after_maxiter::Float64
    result
    function cmpfit()
        out = new(CMPFit.Config(), 1e-4, nothing)
        out.config.maxiter = 1000
        return out
    end
end

function fit(mzer::cmpfit, fp::AbstractFitProblem)
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

    residuals(fp, guess)
    last_fitstat = sum(abs2, residuals(fp))
    prog = ProgressUnknown(desc="Model (dof=$(fp.dof)) evaluations:", dt=0.5, showspeed=true, color=:light_black)
    while true
        mzer.result = CMPFit.cmpfit((pvalues) -> begin
                                        ProgressMeter.next!(prog; showvalues=() -> [(:fit_stat, fit_stat(fp))])
                                        residuals(fp, pvalues)
                                    end,
                                    guess, parinfo=parinfo, config=mzer.config)
        if mzer.result.status <= 0
            return MinimizerStatusError("CMPFit status = $(mzer.result.status)")
        end

        if (mzer.result.status == 5)
            Δfitstat = (last_fitstat - mzer.result.bestnorm) / last_fitstat
            if Δfitstat > mzer.ftol_after_maxiter
                println("Reached max. number of iteration but relative Δfitstat = $(Δfitstat) > $(mzer.ftol_after_maxiter), continue minimization...\n")
                last_fitstat = mzer.result.bestnorm
                guess = getfield.(Ref(mzer.result), :param)
                continue
            end
        end

        ProgressMeter.finish!(prog)
        finalize!(fp,
                  getfield.(Ref(mzer.result), :param),
                  getfield.(Ref(mzer.result), :perror))

        if mzer.result.status == 2
            return MinimizerStatusWarn("CMPFit status = 2 may imply one (or more) guess values are too far from optimum")
        elseif mzer.result.status == 5
            return MinimizerStatusWarn("CMPFit status = 5, reached maximum allowed number of iteration.")
        end
        return MinimizerStatusOK()
    end
end
