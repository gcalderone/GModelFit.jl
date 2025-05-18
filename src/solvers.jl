module Solvers

using ProgressMeter

export AbstractSolverStatus, SolverStatusOK, SolverStatusWarn, SolverStatusError, AbstractSolver, WrapSolver, solve!, lsqfit, cmpfit

import ..GModelFit: FitProblem, free_params, nfree, ndata, fitstat, evaluate_residuals!, set_bestfit!
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

mutable struct WrapSolver{T <: Union{AbstractSolver, NonlinearSolve.NonlinearSolveBase.AbstractNonlinearSolveAlgorithm}}
    solver::T
    result
    WrapSolver(solver::T) where T = new{T}(solver, nothing)
end

solve!(fitprob::FitProblem, solver::AbstractSolver; kws...) = solve!(fitprob, WrapSolver(solver); kws...)
solve!(fitprob::FitProblem, solver::NonlinearSolve.NonlinearSolveBase.AbstractNonlinearSolveAlgorithm; kws...) = solve!(fitprob, WrapSolver(solver); kws...)


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
                evaluate_residuals!(du, shared.fp, pvalues)
            end
        end
    else
        funct = let prog=prog, shared=shared
            pvalues -> begin
                ProgressMeter.next!(prog; showvalues=() -> [(:fitstat, fitstat(shared.fp))])
                return evaluate_residuals!(shared.output, shared.fp, pvalues)
            end
        end
    end
    return prog, shared, funct
end


# --------------------------------------------------------------------
import LsqFit

struct lsqfit <: AbstractSolver end

function solve!(fitprob::FitProblem, wrap::WrapSolver{lsqfit})
    prog, shared, funct = eval_funct(fitprob)
    wrap.result = LsqFit.curve_fit((dummy, pvalues) -> funct(pvalues),
                                   1.:ndata(fitprob), fill(0., ndata(fitprob)),
                                   shared.guess, lower=shared.lowb, upper=shared.highb)
    ProgressMeter.finish!(prog)
    if !wrap.result.converged
        return SolverStatusError("Not converged")
    end

    set_bestfit!(fitprob, getfield.(Ref(wrap.result), :param), LsqFit.stderror(wrap.result))
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

function solve!(fitprob::FitProblem, wrap::WrapSolver{cmpfit})
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
        wrap.result = CMPFit.cmpfit(funct, guess, parinfo=parinfo, config=wrap.solver.config)
        if wrap.result.status <= 0
            return SolverStatusError("CMPFit status = $(wrap.result.status)")
        end

        if (wrap.result.status == 5)
            Δfitstat = (last_fitstat - wrap.result.bestnorm) / last_fitstat
            if Δfitstat > wrap.solver.ftol_after_maxiter
                println("Reached max. number of iteration but relative Δfitstat = $(Δfitstat) > $(wrap.ftol_after_maxiter), continue minimization...\n")
                last_fitstat = wrap.result.bestnorm
                guess = getfield.(Ref(wrap.result), :param)
                continue
            end
        end
        break
    end
    ProgressMeter.finish!(prog)

    set_bestfit!(fitprob,
                 getfield.(Ref(wrap.result), :param),
                 getfield.(Ref(wrap.result), :perror))

    if wrap.result.status == 2
        return SolverStatusWarn("CMPFit status = 2 may imply one (or more) guess values are too far from optimum")
    elseif wrap.result.status == 5
        return SolverStatusWarn("CMPFit status = 5, reached maximum allowed number of iteration.")
    end
    return SolverStatusOK()
end


# --------------------------------------------------------------------
function solve!(fitprob::FitProblem, wrap::WrapSolver{T}) where T <: NonlinearSolve.NonlinearSolveBase.AbstractNonlinearSolveAlgorithm
    prog, shared, funct = eval_funct(fitprob, nonlinearsolve=true)

    wrap.result = NonlinearSolve.solve(NonlinearSolve.NonlinearLeastSquaresProblem(
        NonlinearSolve.NonlinearFunction(funct, resid_prototype = zeros(ndata(fitprob))),
        shared.guess, shared), wrap.solver)
    ProgressMeter.finish!(prog)

    # TODO: AbstractSolverStatus
    set_bestfit!(fitprob, wrap.result.u, wrap.result.u .* 0)
    return SolverStatusOK()
end

end
