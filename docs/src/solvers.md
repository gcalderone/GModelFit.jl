# Minimizers

The **GModelFit.jl** main purpose is to act as an high-level interface between the user and the underlying minimizer.
Currently, two [non-linear least squares](https://en.wikipedia.org/wiki/Non-linear_least_squares) minimizers are supported:
- [LsqFit](https://github.com/JuliaNLSolvers/LsqFit.jl);
- [CMPFit](https://github.com/gcalderone/CMPFit.jl).

More minimizers may be added in the future.

To choose a specific minimizer use the `minimizer=` keyword when invoking the [`fit()`](@ref) function, e.g. `minimizer=GModelFit.lsqfit()` or `minimizer=GModelFit.cmpfit()`.  If the keyword is not provided the `lsqfit()` minimizer is used.


### Examples

```@example abc
using GModelFit
model = Model(:main => @fd (x, T=3.14) -> sin.(x ./ T) ./ (x ./ T))
data = GModelFit.mock(Measures, model, Domain(1:0.1:50), seed=1)
bestfit, stats = fit(model, data, lsqfit())
println(); # hide
```
or
```@example abc
bestfit, stats = fit(model, data, cmpfit())
println(); # hide
```

The above minimizers typically provide the same results, although in some complex case the [CMPFit](https://github.com/gcalderone/CMPFit.jl) may be more robust and less sensitive to initial guess parameters.


## The `cmpfit()` minimizer

The `cmpfit()` minimizer allows to specify several options to fine-tune the minimizer behaviour.  Specifically:
- the `CMPFit.Config` structure allows to specify the convergence criteria, the maximum number of iterations, etc. (see the "CONFIGURING MPFIT()" section [here](https://pages.physics.wisc.edu/~craigm/idl/cmpfit.html);
- the `ftol_after_maxiter` allows to specify a threshold on the relative difference in fit statistics before and after the `mpfit()` execution.  If the latter terminates because the maximum number of iterations has been reached, and the relative difference in fit statistics is still greater than `ftol_after_maxiter` the minimization process will continue.  E.g.:

  ```@example abc
  using GModelFit
  dom = Domain(1:0.1:50)
  model = Model(:main => @fd (x, T=3.14) -> sin.(x ./ T) ./ (x ./ T))
  data = GModelFit.mock(Measures, model, dom, seed=1)
  
  # Set minimizer options
  solver = GModelFit.cmpfit()
  solver.config.maxiter = 1
  solver.ftol_after_maxiter = 1e-8
  
  # Run the fit
  model[:main].T.val = 10  # guess value, purposely far from true one
  bestfit, stats = fit(model, data, solver)
  println(); # hide
  ```
