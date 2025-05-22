# Solvers

The **GModelFit.jl** main purpose is to act as an high-level interface between the user and the underlying solver.
Currently supported solvers are:
- [LsqFit](https://julianlsolvers.github.io/LsqFit.jl/latest/) (default);
- [CMPFit](https://github.com/gcalderone/CMPFit.jl);
- [NonlinearSolve](https://docs.sciml.ai/NonlinearSolve/stable/).

More solvers may be added in the future.

To choose a specific solver add a third argument to the [`fit()`](@ref) or [`fit!()`](@ref) functions, e.g. 
```julia
fit(model, data)                 # use default solver: lsqfit
fit(model, data, cmpfit())       # use CMPFIT solver
fit(model, data, TrustRegion())  # use NonlinearSolve.TrustRegion solver
```


### Examples

```@example abc
using GModelFit
model = Model(:main => @fd (x, T=3.14) -> sin.(x ./ T) ./ (x ./ T))
data = GModelFit.mock(Measures, model, Domain(1:0.1:50), seed=1)
bestfit, fsumm = fit(model, data, lsqfit())
```
or
```@example abc
bestfit, fsumm = fit(model, data, cmpfit())
```
or
```@example abc
using NonlinearSolve
bestfit, fsumm = fit(model, data, NewtonRaphson())
```

The above solvers typically provide the same results, although in some complex case the results may differ.

!!! warning
    Unlike LsqFit and CMPFIT, the solvers from [NonlinearSolve](https://docs.sciml.ai/NonlinearSolve/stable/) do not provide best fit parameter uncertainties.


## The `cmpfit()` solver

The `cmpfit()` solver allows to specify several options to fine-tune the solver behaviour.  Specifically:
- the `CMPFit.Config` structure allows to specify the convergence criteria, the maximum number of iterations, etc. (see the "CONFIGURING MPFIT()" section [here](https://pages.physics.wisc.edu/~craigm/idl/cmpfit.html);
- the `ftol_after_maxiter` allows to specify a threshold on the relative difference in fit statistics before and after the `mpfit()` execution.  If the latter terminates because the maximum number of iterations has been reached, and the relative difference in fit statistics is still greater than `ftol_after_maxiter` the minimization process will continue.  E.g.:

  ```@example abc
  using GModelFit
  dom = Domain(1:0.1:50)
  model = Model(:main => @fd (x, T=3.14) -> sin.(x ./ T) ./ (x ./ T))
  data = GModelFit.mock(Measures, model, dom, seed=1)
  
  # Set solver options
  solver = GModelFit.cmpfit()
  solver.config.maxiter = 1
  solver.ftol_after_maxiter = 1e-8
  
  # Run the fit
  model[:main].T.val = 10  # guess value, purposely far from true one
  bestfit, fsumm = fit(model, data, solver)
  println(); # hide
  ```
