module Solvers

using ProgressMeter, Dates

export AbstractSolverStatus, SolverStatusOK, SolverStatusWarn, SolverStatusError, FitSummary, AbstractSolver, solve!

import ..GModelFit: FitProblem, free_params, nfree, ndata, fitstat, update_eval!, set_bestfit!
import ForwardDiff

# --------------------------------------------------------------------
abstract type AbstractSolverStatus end

struct SolverStatusOK <: AbstractSolverStatus
end

struct SolverStatusWarn <: AbstractSolverStatus
    message::String
end

struct SolverStatusError <: AbstractSolverStatus
    message::String
end

# --------------------------------------------------------------------
"""
    FitSummary

A structure summarizing the results of a fitting process.

# Fields:
- `start::DateTime`: timestamp at the beginning of the fitting process;
- `elapsed::Float64`: elapsed time (in seconds);
- `ndata::Int`: number of data empirical points;
- `nfree::Int`: number of free parameters;
- `fitstat::Float64`: fit statistics (equivalent ro reduced χ^2 for `Measures` objects);
- `status`: minimization process status (tells whether convergence criterion has been satisfied, or if an error has occurred during fitting);
- `solver_retval`: solver return value.

Note: the `solver_retval` field can not be serialized, will contain `nothing` when deserialized.
"""
struct FitSummary
    start::DateTime
    elapsed::Float64
    ndata::Int
    nfree::Int
    fitstat::Float64
    status::AbstractSolverStatus
    solver_retval
end

FitSummary(fitprob::FitProblem, status::AbstractSolverStatus, start::DateTime, elapsed::TimePeriod, solver_retval=nothing) =
    FitSummary(start, Dates.value(convert(Millisecond, elapsed)) / 1000., ndata(fitprob), nfree(fitprob), fitstat(fitprob), status, solver_retval)


# --------------------------------------------------------------------
abstract type AbstractSolver end


# --------------------------------------------------------------------
function eval_funct(fitprob::FitProblem)
    params = free_params(fitprob)
    guess  = getfield.(params, :val)
    lowb   = getfield.(params, :low)
    highb  = getfield.(params, :high)

    prog = ProgressUnknown(desc="Nfree=$(nfree(fitprob)), evaluations:", dt=0.5, showspeed=true, color=:light_black)

    shared = (fp=fitprob, start=now(), guess=guess, lowb=lowb, highb=highb)
    funct = let prog=prog, shared=shared
        pvalues -> begin
            ProgressMeter.next!(prog; showvalues=() -> [(:fitstat, fitstat(shared.fp))])
            return update_eval!(shared.fp, pvalues)
        end
    end
    return prog, shared, funct
end


# --------------------------------------------------------------------
struct dry <: AbstractSolver end

function solve!(fitprob::FitProblem, solver::dry)
    prog, shared, funct = eval_funct(fitprob)
    funct(shared.guess)
    ProgressMeter.finish!(prog)
    set_bestfit!(fitprob, shared.guess, shared.guess .* 0.)
    return FitSummary(fitprob, SolverStatusWarn("dry solver"), shared.start, now() - shared.start, nothing)
end


# --------------------------------------------------------------------
import LsqFit

struct lsqfit <: AbstractSolver end

function solve!(fitprob::FitProblem, solver::lsqfit)
    prog, shared, funct = eval_funct(fitprob)
    solver_retval = LsqFit.curve_fit((dummy, pvalues) -> funct(pvalues),
                              1.:ndata(fitprob), fill(0., ndata(fitprob)),
                              shared.guess, lower=shared.lowb, upper=shared.highb)
    ProgressMeter.finish!(prog)

    status = SolverStatusOK()
    if !solver_retval.converged
        status = SolverStatusError("Not converged")
    end

    set_bestfit!(fitprob, getfield.(Ref(solver_retval), :param), LsqFit.stderror(solver_retval))
    return FitSummary(fitprob, status, shared.start, now() - shared.start, solver_retval)
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

function solve!(fitprob::FitProblem, solver::cmpfit)
    prog, shared, funct = eval_funct(fitprob)
    parinfo = CMPFit.Parinfo(length(shared.guess))
    for i in 1:length(shared.guess)
        llow  = isfinite(shared.lowb[i])   ?  1  :  0
        lhigh = isfinite(shared.highb[i])  ?  1  :  0
        parinfo[i].limited = (llow, lhigh)
        parinfo[i].limits  = (shared.lowb[i], shared.highb[i])
    end

    last_fitstat = fitstat(fitprob)
    guess = shared.guess
    solver_retval = nothing
    status = SolverStatusOK()
    while true
        solver_retval = CMPFit.cmpfit(funct, guess, parinfo=parinfo, config=solver.config)
        if solver_retval.status <= 0
            status = SolverStatusError("CMPFit status = $(solver_retval.status)")
            break
        end

        if (solver_retval.status == 5)
            Δfitstat = (last_fitstat - solver_retval.bestnorm) / last_fitstat
            if Δfitstat > solver.ftol_after_maxiter
                println("Reached max. number of iteration but relative Δfitstat = $(Δfitstat) > $(solver.ftol_after_maxiter), continue minimization...\n")
                last_fitstat = solver_retval.bestnorm
                guess = getfield.(Ref(solver_retval), :param)
                continue
            end
        end
        break
    end
    ProgressMeter.finish!(prog)

    set_bestfit!(fitprob,
                 getfield.(Ref(solver_retval), :param),
                 getfield.(Ref(solver_retval), :perror))

    if solver_retval.status == 2
        status = SolverStatusWarn("CMPFit status = 2 may imply one (or more) guess values are too far from optimum")
    elseif solver_retval.status == 5
        status = SolverStatusWarn("CMPFit status = 5, reached maximum allowed number of iteration.")
    end
    return FitSummary(fitprob, status, shared.start, now() - shared.start, solver_retval)
end


# --------------------------------------------------------------------
import CurveFit

struct curvefit <: AbstractSolver
    alg
    curvefit(alg=nothing) = new(alg)
end

function solve!(fitprob::FitProblem, solver::curvefit)
    prog, shared, funct = eval_funct(fitprob)

    p = CurveFit.NonlinearCurveFitProblem((pvalues, dummy) -> funct(pvalues),
                                          shared.guess, 1.:ndata(fitprob))
    solver_retval = isnothing(solver.alg)  ?  CurveFit.solve(p)  :  CurveFit.solve(p, solver.alg)
    ProgressMeter.finish!(prog)

    status = CurveFit.isconverged(solver_retval)  ?  SolverStatusOK()  :  SolverStatusError("Not converged")
    set_bestfit!(fitprob, solver_retval.u, solver_retval.u, CurveFit.stderror(solver_retval))
    return FitSummary(fitprob, status, shared.start, now() - shared.start, solver_retval)
end

end
