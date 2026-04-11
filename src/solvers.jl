module Solvers

using ProgressMeter, Dates

export AbstractSolverStatus, SolverStatusOK, SolverStatusWarn, SolverStatusError, AbstractSolver, need_autodiff, maximize!

import ..GModelFit: AbstractData, Likelihood, getfreepars, nfree, ndata, gofstat, evaluate_resid!, _set_bestfit!
import ForwardDiff

# --------------------------------------------------------------------
abstract type AbstractSolverStatus{T} end

struct SolverStatusOK{T} <: AbstractSolverStatus{T}
    tstart::DateTime
    elapsed::Float64
    status::T
    SolverStatusOK(shared, status::T) where T = new{typeof(status)}(shared.tstart, Dates.value(convert(Millisecond, now() - shared.tstart)) / 1000., status)
end

struct SolverStatusWarn{T} <: AbstractSolverStatus{T}
    tstart::DateTime
    elapsed::Float64
    status::T
    message::String
    SolverStatusWarn(shared, message::String, status::T) where T = new{typeof(status)}(shared.tstart, Dates.value(convert(Millisecond, now() - shared.tstart)) / 1000., status, message)
end

struct SolverStatusError{T} <: AbstractSolverStatus{T}
    tstart::DateTime
    elapsed::Float64
    status::T
    message::String
    SolverStatusError(shared, message::String, status::T) where T = new{typeof(status)}(shared.tstart, Dates.value(convert(Millisecond, now() - shared.tstart)) / 1000., status, message)
end

# --------------------------------------------------------------------
"""
    FitSummary

A structure summarizing the results of a fitting process.

# Fields:
- `tstart::DateTime`: timestamp at the beginning of the fitting process;
- `elapsed::Float64`: elapsed time (in seconds);
- `status`: minimization process status (tells whether convergence criterion has been satisfied, or if an error has occurred during fitting);
- `solver_retval`: solver return value.

Note: the `solver_retval` field can not be serialized, will contain `nothing` when deserialized.
"""

# --------------------------------------------------------------------
abstract type AbstractSolver end

# --------------------------------------------------------------------
function eval_funct(lh::Likelihood)
    params = collect(values(getfreepars(lh)))
    guess  = getfield.(params, :val)
    lowb   = getfield.(params, :low)
    highb  = getfield.(params, :high)

    prog = ProgressUnknown(desc="Nfree=$(nfree(lh)), evaluations:", dt=0.5, showspeed=true, color=:light_black)

    shared = (lh=lh, tstart=now(), guess=guess, lowb=lowb, highb=highb)
    funct = let prog=prog, shared=shared
        pvalues -> begin
            ProgressMeter.next!(prog; showvalues=() -> begin
                                    return [(:gofstat, ForwardDiff.value(gofstat(shared.lh)))]
                                end)
            return evaluate_resid!(shared.lh, pvalues)
        end
    end
    return prog, shared, funct
end


# --------------------------------------------------------------------
struct dry <: AbstractSolver end
need_autodiff(::dry) = false

function maximize!(lh::Likelihood, solver::dry)
    prog, shared, funct = eval_funct(lh)
    funct(shared.guess)
    ProgressMeter.finish!(prog)
    _set_bestfit!(lh, shared.guess, shared.guess .* 0.)
    return SolverStatusWarn(shared, "dry solver", nothing)
end


# --------------------------------------------------------------------
import LsqFit

struct lsqfit <: AbstractSolver end
need_autodiff(::lsqfit) = false

function maximize!(lh::Likelihood, solver::lsqfit)
    prog, shared, funct = eval_funct(lh)
    solver_retval = LsqFit.curve_fit((dummy, pvalues) -> funct(pvalues),
                                     1.:ndata(lh), fill(0., ndata(lh)),
                                     shared.guess, lower=shared.lowb, upper=shared.highb,
                                     autodiff=LsqFit.AutoFiniteDiff())
    ProgressMeter.finish!(prog)
    if solver_retval.converged
        status = SolverStatusOK(shared, solver_retval)
    else
        status = SolverStatusError(shared, "Not converged", solver_retval)
    end
    _set_bestfit!(lh, getfield.(Ref(solver_retval), :param), LsqFit.stderror(solver_retval))
    return status
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
improvements (ftol_after_maxiter) to be checked after the solver
iterated for the maximum allowed number of times.
=#

mutable struct cmpfit <: AbstractSolver
    config::CMPFit.Config
    ftol_after_maxiter::Float64

    function cmpfit()
        out = new(CMPFit.Config(), 1e-4)
        out.config.maxiter = 1000
        return out
    end
end
need_autodiff(::cmpfit) = false

function maximize!(lh::Likelihood, solver::cmpfit)
    prog, shared, funct = eval_funct(lh)
    parinfo = CMPFit.Parinfo(length(shared.guess))
    for i in 1:length(shared.guess)
        llow  = isfinite(shared.lowb[i])   ?  1  :  0
        lhigh = isfinite(shared.highb[i])  ?  1  :  0
        parinfo[i].limited = (llow, lhigh)
        parinfo[i].limits  = (shared.lowb[i], shared.highb[i])
    end

    last_gofstat = gofstat(lh)
    guess = shared.guess
    solver_retval = nothing
    while true
        solver_retval = CMPFit.cmpfit(funct, guess, parinfo=parinfo, config=solver.config)
        (solver_retval.status <= 0)  &&  break

        if (solver_retval.status == 5)
            Δgofstat = (last_gofstat - solver_retval.bestnorm) / last_gofstat
            if Δgofstat > solver.ftol_after_maxiter
                println("Reached max. number of iteration but relative Δgofstat = $(Δgofstat) > $(solver.ftol_after_maxiter), continue minimization...\n")
                last_gofstat = solver_retval.bestnorm
                guess = getfield.(Ref(solver_retval), :param)
                continue
            end
        end
        break
    end
    ProgressMeter.finish!(prog)

    _set_bestfit!(lh,
                  getfield.(Ref(solver_retval), :param),
                  getfield.(Ref(solver_retval), :perror))

    if solver_retval.status <= 0
        return SolverStatusError(shared, "CMPFit status = $(solver_retval.status)", solver_retval)
    elseif solver_retval.status == 2
        return SolverStatusWarn(shared, "CMPFit status = 2 may imply one (or more) guess values are too far from optimum", solver_retval)
    elseif solver_retval.status == 5
        return SolverStatusWarn(shared, "CMPFit status = 5, reached maximum allowed number of iteration.", solver_retval)
    end
    return SolverStatusOK(shared, solver_retval)
end


# --------------------------------------------------------------------
import CurveFit

struct curvefit <: AbstractSolver
    alg
    curvefit(alg=nothing) = new(alg)
end
need_autodiff(::curvefit) = true

function maximize!(lh::Likelihood, solver::curvefit)
    prog, shared, funct = eval_funct(lh)

    p = CurveFit.NonlinearCurveFitProblem((pvalues, dummy) -> funct(pvalues),
                                          shared.guess, 1.:ndata(lh), fill(0., ndata(lh)))
    solver_retval = isnothing(solver.alg)  ?  CurveFit.solve(p)  :  CurveFit.solve(p, solver.alg)
    ProgressMeter.finish!(prog)

    _set_bestfit!(lh, solver_retval.u, CurveFit.stderror(solver_retval))
    funct(solver_retval.u) # update evaluation to best fit values (these may have been modified when invoking CurveFit.stderror)
    if CurveFit.isconverged(solver_retval)
        return SolverStatusOK(shared, solver_retval)
    end
    return SolverStatusError(shared, "Not converged", solver_retval)
end

end
