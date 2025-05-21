module Solvers

using ProgressMeter

export AbstractSolverStatus, SolverStatusOK, SolverStatusWarn, SolverStatusError, AbstractSolver, CaptureSolution, solve!, lsqfit, cmpfit

import ..GModelFit: FitProblem, free_params, nfree, ndata, fitstat, update_eval!, set_bestfit!
import NonlinearSolve

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
abstract type AbstractSolver end

mutable struct CaptureSolution{T <: Union{AbstractSolver, NonlinearSolve.NonlinearSolveBase.AbstractNonlinearSolveAlgorithm}}
    solver::T
    sol
    CaptureSolution(solver::T) where T = new{T}(solver, nothing)
end

solve!(fitprob::FitProblem, solver::Union{AbstractSolver, NonlinearSolve.NonlinearSolveBase.AbstractNonlinearSolveAlgorithm}; kws...) =
    solve!(fitprob, CaptureSolution(solver); kws...)

# --------------------------------------------------------------------
function eval_funct(fitprob::FitProblem; nonlinearsolve=false)
    params = free_params(fitprob)
    guess  = getfield.(params, :val)
    lowb   = getfield.(params, :low)
    highb  = getfield.(params, :high)

    prog = ProgressUnknown(desc="Nfree=$(nfree(fitprob)), evaluations:", dt=0.5, showspeed=true, color=:light_black)

    shared = (fp=fitprob, guess=guess, lowb=lowb, highb=highb, output=Vector{Float64}(undef, ndata(fitprob)))
    if nonlinearsolve
        funct = let prog=prog
            (du, pvalues, shared) -> begin
                ProgressMeter.next!(prog; showvalues=() -> [(:fitstat, fitstat(shared.fp))])
                update_eval!(shared.fp, du, pvalues)
            end
        end
    else
        funct = let prog=prog, shared=shared
            pvalues -> begin
                ProgressMeter.next!(prog; showvalues=() -> [(:fitstat, fitstat(shared.fp))])
                return update_eval!(shared.fp, shared.output, pvalues)
            end
        end
    end
    return prog, shared, funct
end


# --------------------------------------------------------------------
import LsqFit

struct lsqfit <: AbstractSolver end

function solve!(fitprob::FitProblem, capture::CaptureSolution{lsqfit})
    prog, shared, funct = eval_funct(fitprob)
    capture.sol = LsqFit.curve_fit((dummy, pvalues) -> funct(pvalues),
                                   1.:ndata(fitprob), fill(0., ndata(fitprob)),
                                   shared.guess, lower=shared.lowb, upper=shared.highb)
    ProgressMeter.finish!(prog)
    if !capture.sol.converged
        return SolverStatusError("Not converged")
    end

    set_bestfit!(fitprob, getfield.(Ref(capture.sol), :param), LsqFit.stderror(capture.sol))
    return SolverStatusOK()
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

function solve!(fitprob::FitProblem, capture::CaptureSolution{cmpfit})
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
    while true
        capture.sol = CMPFit.cmpfit(funct, guess, parinfo=parinfo, config=capture.solver.config)
        if capture.sol.status <= 0
            return SolverStatusError("CMPFit status = $(capture.sol.status)")
        end

        if (capture.sol.status == 5)
            Δfitstat = (last_fitstat - capture.sol.bestnorm) / last_fitstat
            if Δfitstat > capture.solver.ftol_after_maxiter
                println("Reached max. number of iteration but relative Δfitstat = $(Δfitstat) > $(capture.ftol_after_maxiter), continue minimization...\n")
                last_fitstat = capture.sol.bestnorm
                guess = getfield.(Ref(capture.sol), :param)
                continue
            end
        end
        break
    end
    ProgressMeter.finish!(prog)

    set_bestfit!(fitprob,
                 getfield.(Ref(capture.sol), :param),
                 getfield.(Ref(capture.sol), :perror))

    if capture.sol.status == 2
        return SolverStatusWarn("CMPFit status = 2 may imply one (or more) guess values are too far from optimum")
    elseif capture.sol.status == 5
        return SolverStatusWarn("CMPFit status = 5, reached maximum allowed number of iteration.")
    end
    return SolverStatusOK()
end


# --------------------------------------------------------------------
function solve!(fitprob::FitProblem, capture::CaptureSolution{T}) where T <: NonlinearSolve.NonlinearSolveBase.AbstractNonlinearSolveAlgorithm
    prog, shared, funct = eval_funct(fitprob, nonlinearsolve=true)

    capture.sol = NonlinearSolve.solve(NonlinearSolve.NonlinearLeastSquaresProblem(
        NonlinearSolve.NonlinearFunction(funct, resid_prototype = zeros(ndata(fitprob))),
        shared.guess, shared), capture.solver)
    ProgressMeter.finish!(prog)

    set_bestfit!(fitprob, capture.sol.u, capture.sol.u .* NaN)
    if NonlinearSolve.SciMLBase.successful_retcode(capture.sol.retcode)
        return SolverStatusOK()
    end
    return SolverStatusError(string(capture.sol.retcode))
end

end
