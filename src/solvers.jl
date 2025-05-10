module Solvers

using ProgressMeter

export AbstractSolverStatus, SolverStatusOK, SolverStatusWarn, SolverStatusError, AbstractSolver, WrapSolver, solve!, cmpfit

import ..GModelFit: FitProblem, free_params, nfree, ndata, fitstat, evaluate!, set_bestfit!
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

solve!(fitprob::FitProblem, solver::AbstractSolver) = solve!(fitprob, WrapSolver(solver))
solve!(fitprob::FitProblem, solver::NonlinearSolve.NonlinearSolveBase.AbstractNonlinearSolveAlgorithm) = solve!(fitprob, WrapSolver(solver))


# --------------------------------------------------------------------
import LsqFit

struct lsqfit <: AbstractSolver end

function solve!(fitprob::FitProblem, wrap::WrapSolver{lsqfit})
    params = free_params(fitprob)
    prog = ProgressUnknown(desc="Model (#free=$(nfree(fitprob))) evaluations:", dt=0.5, showspeed=true, color=:light_black)
    wrap.result = LsqFit.curve_fit((dummy, pvalues) -> begin
                                  ProgressMeter.next!(prog; showvalues=() -> [(:fitstat, fitstat(fitprob))])
                                  return evaluate!(fitprob, pvalues)
                              end,
                              1.:ndata(fitprob), fill(0., ndata(fitprob)),
                              getfield.(params, :val),
                              lower=getfield.(params, :low),
                              upper=getfield.(params, :high))
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

struct cmpfit <: AbstractSolver
    config::CMPFit.Config
    ftol_after_maxiter::Float64

    function cmpfit()
        out = new(CMPFit.Config(), 1e-4)
        out.config.maxiter = 1000
        return out
    end
end

function solve!(fitprob::FitProblem, wrap::WrapSolver{cmpfit})
    params = free_params(fitprob)
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

    evaluate!(fitprob, guess)
    last_fitstat = fitstat(fitprob)
    prog = ProgressUnknown(desc="Model (#free=$(nfree(fitprob))) evaluations:", dt=0.5, showspeed=true, color=:light_black)
    while true
        wrap.result = CMPFit.cmpfit((pvalues) -> begin
                                        ProgressMeter.next!(prog; showvalues=() -> [(:fitstat, fitstat(fitprob))])
                                        return evaluate!(fitprob, pvalues)
                                    end,
                                    guess, parinfo=parinfo, config=wrap.solver.config)
        if wrap.result.status <= 0
            return SolverStatusError("CMPFit status = $(wrap.result.status)")
        end

        if (wrap.result.status == 5)
            Δfitstat = (last_fitstat - wrap.result.bestnorm) / last_fitstat
            if Δfitstat > wrap.ftol_after_maxiter
                println("Reached max. number of iteration but relative Δfitstat = $(Δfitstat) > $(wrap.ftol_after_maxiter), continue minimization...\n")
                last_fitstat = wrap.result.bestnorm
                guess = getfield.(Ref(wrap.result), :param)
                continue
            end
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
end


# --------------------------------------------------------------------
function solve!(fitprob::FitProblem, wrap::WrapSolver{T}) where T <: NonlinearSolve.NonlinearSolveBase.AbstractNonlinearSolveAlgorithm
    params = free_params(fitprob)

    prog = ProgressUnknown(desc="Model (#free=$(nfree(fitprob))) evaluations:", dt=0.5, showspeed=true, color=:light_black)
    function local_evaluate!(du, u, fp)
        ProgressMeter.next!(prog; showvalues=() -> [(:fitstat, fitstat(fp))])
        du .= evaluate!(fp, u)
    end

    wrap.result = NonlinearSolve.solve(NonlinearSolve.NonlinearLeastSquaresProblem(
        NonlinearSolve.NonlinearFunction(local_evaluate!, resid_prototype = zeros(ndata(fitprob))),
        getfield.(params, :val), fitprob),
                                       wrap.solver)
    ProgressMeter.finish!(prog)

    # TODO: AbstractSolverStatus
    set_bestfit!(fitprob, wrap.result.u, wrap.result.u .* 0)
    return SolverStatusOK()
end

end
