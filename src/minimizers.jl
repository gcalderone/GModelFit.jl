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
abstract type AbstractResiduals{T <: AbstractMeasures, M <: AbstractMinimizer} end


# --------------------------------------------------------------------
struct dry <: AbstractMinimizer; end
function minimize!(resid::AbstractResiduals{Measures{N}, dry}) where N
    params = free_params(resid)
    residuals(resid, getfield.(params, :val))
    finalize!(resid,
              getfield.(params, :val),
              fill(NaN, length(params)))
    return MinimizerStatusDry()
end


# --------------------------------------------------------------------
import LsqFit

mutable struct lsqfit <: AbstractMinimizer
    result
    lsqfit() = new(nothing)
end

function minimize!(resid::AbstractResiduals{Measures{N}, lsqfit}) where N
    params = free_params(resid)
    ndata = length(residuals(resid))
    prog = ProgressUnknown(desc="Model (dof=$(resid.dof)) evaluations:", dt=0.5, showspeed=true, color=:light_black)
    resid.mzer.result = LsqFit.curve_fit((dummy, pvalues) -> begin
                                             ProgressMeter.next!(prog; showvalues=() -> [(:fit_stat, fit_stat(resid))])
                                             residuals(resid, pvalues)
                                         end,
                                         1.:ndata, fill(0., ndata),
                                         getfield.(params, :val),
                                         lower=getfield.(params, :low),
                                         upper=getfield.(params, :high))
    ProgressMeter.finish!(prog)
    if !resid.mzer.result.converged
        return MnimizerStatusError("Not converged")
    end

    finalize!(resid, getfield.(Ref(resid.mzer.result), :param), LsqFit.stderror(resid.mzer.result))
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

function minimize!(resid::AbstractResiduals{Measures{N}, cmpfit}) where N
    params = free_params(resid)
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

    residuals(resid, guess)
    last_fitstat = sum(abs2, residuals(resid))
    prog = ProgressUnknown(desc="Model (dof=$(resid.dof)) evaluations:", dt=0.5, showspeed=true, color=:light_black)
    while true
        resid.mzer.result = CMPFit.cmpfit((pvalues) -> begin
                                              ProgressMeter.next!(prog; showvalues=() -> [(:fit_stat, fit_stat(resid))])
                                              residuals(resid, pvalues)
                                          end,
                                          guess, parinfo=parinfo, config=resid.mzer.config)
        if resid.mzer.result.status <= 0
            return MinimizerStatusError("CMPFit status = $(resid.mzer.result.status)")
        end

        if (resid.mzer.result.status == 5)
            Δfitstat = (last_fitstat - resid.mzer.result.bestnorm) / last_fitstat
            if Δfitstat > resid.mzer.ftol_after_maxiter
                println("Reached max. number of iteration but relative Δfitstat = $(Δfitstat) > $(resid.mzer.ftol_after_maxiter), continue minimization...\n")
                last_fitstat = resid.mzer.result.bestnorm
                guess = getfield.(Ref(resid.mzer.result), :param)
                continue
            end
        end

        ProgressMeter.finish!(prog)
        finalize!(resid,
                  getfield.(Ref(resid.mzer.result), :param),
                  getfield.(Ref(resid.mzer.result), :perror))

        if resid.mzer.result.status == 2
            return MinimizerStatusWarn("CMPFit status = 2 may imply one (or more) guess values are too far from optimum")
        elseif resid.mzer.result.status == 5
            return MinimizerStatusWarn("CMPFit status = 5, reached maximum allowed number of iteration.")
        end
        return MinimizerStatusOK()
    end
end
