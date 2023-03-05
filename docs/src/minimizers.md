# Minimizers

The **GFit.jl** main purpose is to act as an high-level interface between the user and the underlying minimizer.
Currently, two [non-linear least squares](https://en.wikipedia.org/wiki/Non-linear_least_squares) minimizers are supported:
- [LsqFit](https://github.com/JuliaNLSolvers/LsqFit.jl);
- [CMPFit](https://github.com/gcalderone/CMPFit.jl);

and more may be added in the future.

To specify a minimizer simply use the `minimizer=` keyword when invoking the [`fit()`](@ref) function, e.g. `minimizer=GFit.lsqfit()` or `minimizer=GFit.cmpfit()`. If the keyword is not provided the `lsqfit()` minimizer is used.

There is also a dummy minimizer, `GFit.dry()`, whose purpose is to compare the model and the data, and to generate a `FitStats` object without modifying the model parameters.  The dry minimizer is used by the [`compare()`](@ref) function.


### Examples

```@example abc
using GFit
dom = Domain(1:0.1:50)
model = Model(dom, :main => @λ (x, T=3.14) -> sin.(x ./ T) ./ (x ./ T))
data = GFit.mock(Measures, model, seed=1)
best, fitstats = fit(model, data, minimizer=GFit.lsqfit())
println(); # hide
```
or
```@example abc
best, fitstats = fit(model, data, minimizer=GFit.cmpfit())
println(); # hide
```

The above minimizers typically provide the same results, although in some complex case the [CMPFit](https://github.com/gcalderone/CMPFit.jl) may be more robust and less sensitive to initial guess parameters.


## The `cmpfit()` minimizer

The `cmpfit()` minimizer allows to specify several options to fine-tune the minimizer behaviour.  Specifically:
- the `CMPFit.Config` structure allows to specify the convergence criteria, the maximum number of iterations, etc. (see the "*CONFIGURING MPFIT()*" section [here](https://pages.physics.wisc.edu/~craigm/idl/cmpfit.html);
- the `Δfitstat_threshold` allows to specify a threshold on the difference in fit statistics before and after the `mpfit()` execution to stop the fitting.  If such difference is greater than the `Δfitstat_threshold` the minimization process continue regardless of the number of iterations.

```@example abc
using GFit
dom = Domain(1:0.1:50)
model = Model(dom, :main => @λ (x, T=3.14) -> sin.(x ./ T) ./ (x ./ T))
data = GFit.mock(Measures, model, seed=1)

# Set minimizer options
mzer = GFit.cmpfit()
mzer.config.maxiter = 1000
mzer.Δfitstat_threshold = 1e-8

# Run the fit
best, fitstats = fit(model, data, minimizer=mzer)
println(); # hide
```
