# ====================================================================
# Minimizers
#
abstract type AbstractMinimizer end
abstract type AbstractMinimizerStatus end


# --------------------------------------------------------------------
struct MinimizerStatusOK <: AbstractMinimizerStatus
    best::Vector{Float64}
    unc::Vector{Float64}
    specific
end

struct MinimizerStatusWarn <: AbstractMinimizerStatus
    best::Vector{Float64}
    unc::Vector{Float64}
    message::String
    specific
end

struct MinimizerStatusError <: AbstractMinimizerStatus
    message::String
    specific
end

function as_string(status::AbstractMinimizerStatus)
    if isa(status, MinimizerStatusOK)
        return (crayon"green", "OK", "")
    elseif isa(status, MinimizerStatusWarn)
        return (crayon"bold yellow", "WARNING!", status.message)
    elseif isa(status, MinimizerStatusWarn)
        return (crayon"bold red", "ERROR!", status.message)
    else
        error("Unsupported type: $(typeof(status))")
    end
end



# --------------------------------------------------------------------
import LsqFit

struct lsqfit <: AbstractMinimizer; end

function minimize(minimizer::lsqfit, func::Function, params::Vector{Parameter})
    ndata = length(func(getfield.(params, :val)))
    res = LsqFit.curve_fit((dummy, pvalues) -> func(pvalues),
                           1.:ndata, fill(0., ndata),
                           getfield.(params, :val),
                           lower=getfield.(params, :low),
                           upper=getfield.(params, :high))

    if !res.converged
        return MinimizerStatusError("Not converged", res)
    end
    return MinimizerStatusOK(getfield.(Ref(res), :param),
                             LsqFit.stderror(res),
                             res)
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
improvements (Δfitstat_theshold) to be checked after the minimizer
iterated for the maximum allowed number of times.
=#

mutable struct cmpfit <: AbstractMinimizer
    config::CMPFit.Config
    Δfitstat_theshold::Float64
    cmpfit() = new(CMPFit.Config(), NaN)
end

function minimize(minimizer::cmpfit, func::Function, params::Vector{Parameter})
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

    last_fitstat = sum(abs2, func(guess))
    while true
        res = CMPFit.cmpfit(func,
                            guess, parinfo=parinfo, config=minimizer.config)

        if res.status <= 0
            return MinimizerStatusError("Status = $(res.status)", res)
        end

        if (res.status == 5)  &&  isfinite(minimizer.Δfitstat_theshold)
            Δfitstat = (last_fitstat - res.bestnorm) / last_fitstat
            if Δfitstat > minimizer.Δfitstat_theshold
                println("\nReached max. number of iteration but relative Δfitstat = $(Δfitstat) > $(minimizer.Δfitstat_theshold), continue minimization...\n")
                last_fitstat = res.bestnorm
                guess = getfield.(Ref(res), :param)
                continue
            end
        end

        if res.status == 2
            return MinimizerStatusWarn(
                getfield.(Ref(res), :param),
                getfield.(Ref(res), :perror),
                "CMPFit status = 2 may imply one (or more) guess values are too far from optimum",
                res)
        end

        return MinimizerStatusOK(getfield.(Ref(res), :param),
                                 getfield.(Ref(res), :perror),
                                 res)
    end
end

