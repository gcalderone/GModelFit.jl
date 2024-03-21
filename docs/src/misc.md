```@setup abc
include("setup.jl")
```

# Miscellaneous

## Generate mock datasets

In some case it is useful to test a model for robustness before the emprical data are available for fitting.  This can be achieved via the [`GModelFit.mock()`](@ref) function, whose purpose is to generate a mock dataset which simulates a measurement process by adding random noise to the foreseen ground-truth.

### Example

```@example abc
using GModelFit

model = Model(:main => @fd (x, T=3.14) -> sin.(x ./ T) ./ (x ./ T))

# Generate a mock dataset on a specific domain
dom = Domain(1:0.1:50)
data = GModelFit.mock(Measures, model, dom, seed=1)

# Fit model against the mock dataset
bestfit, stats = fit(model, data)
```



## Serialization

A few structures, namely [`GModelFit.ModelSnapshot`](@ref), [`GModelFit.FitStats`](@ref) and [`Measures{N}`](@ref), as well as `Vector`(s) of such structures can be *serialized*, i.e. stored in a file using a dedicated JSON format.  The structures can lated be *de-serialized* in a separata Julia session without the need to re-run the fitting process used to create them in the first place.

### Example

In the following we will generate a few **GModelFit.jl** objects and serialized them in a file.
```@example abc
using GModelFit

dom = Domain(1:0.1:50)
model = Model(:main => @fd (x, T=3.14) -> sin.(x ./ T) ./ (x ./ T))
data = GModelFit.mock(Measures, model, dom, seed=1)
bestfit, stats = fit(model, data)

# Serialize objects and save in a file
GModelFit.serialize("save_for_future_use.json", bestfit, stats, data)
println(); # hide
```

The same objects can be de-serialized in a different Julia session:
```@example abc
using GModelFit
bestfit, stats, data = GModelFit.deserialize("save_for_future_use.json")
```



## Quick plot (1D)

The **GModelFit.jl** package implements [**Gnuplot.jl**](https://github.com/gcalderone/Gnuplot.jl/) recipes to display plots of `Measures{1}` and `ModelSnapshot` objects., e.g.:

### Example

Create a model, a mock dataset and run a fit:
```@example abc
using GModelFit

dom = Domain(0:0.01:5)
model = Model(:bkg => GModelFit.OffsetSlope(1, 1, 0.1),
              :l1 => GModelFit.Gaussian(1, 2, 0.2),
              :l2 => GModelFit.Gaussian(1, 3, 0.4),
              :main => SumReducer(:bkg, :l1, :l2))
data = GModelFit.mock(Measures, model, dom)
bestfit, stats = fit(model, data)
println(); # hide
```

A plot of the dataset and of the best fit model can be simply obtained with
```@example abc
using Gnuplot
@gp data bestfit
saveas("gnuplot1") # hide
```
![](assets/gnuplot1.png)

You may also specify axis range, labels, title, etc. using the standard [**Gnuplot.jl**](https://github.com/gcalderone/Gnuplot.jl/) keyword syntax, e.g.:

```@example abc
using Gnuplot
@gp xr=[1, 4.5] xlabel="Wavelength" ylab="Flux" "set key outside" data bestfit
saveas("gnuplot2") # hide
```
![](assets/gnuplot2.png)


## GModelFit internals

(This section deals with **GModelFit.jl** internals, feel free to skip if not interested.)

During minimization a number of internal data structures are created to avoid reallocating heap memory at each iteration.  The most important of such structures are:

- [`GModelFit.CompEval`](@ref): a container to perform component evaluation on a specific domain.  This structure is relevant when defining [Custom components](@ref) as it is used to dispatch component evaluation to the proper `evaluate!` method;

- [`GModelFit.ModelEval`](@ref): a container for a [`Model`](@ref) evaluation on a specific domain.  This structure contains a dictionary of `CompEval` structures for all components in a model, as well as the values of *patched* parameters (see [Parameter constraints](@ref)), and is updated at each iteration of the minimizer to reflect the current model evaluaion.

  An important functionality of the `ModelEval` structure is that it detects the changes in the original `Model` even after it has been created, e.g.:
  ```@example abc
  # Create a Model
  model = Model(:comp1 => @fd (x, p1=1) -> p1 .* x)

  # Wrap the model into a ModelEval to perform evaluation on a specific domain
  dom = Domain(1:5)
  meval = GModelFit.ModelEval(model, dom)
  
  # Evaluate and print the maximum value
  GModelFit.update!(meval)
  println(maximum(GModelFit.last_evaluation(meval)))
  
  # Add a second component to the original model
  model[:comp2] = @fd (x, comp1, p2=1) -> comp1 .+ p2 .* x.^2
  
  # Re-evaluate the ModelEval (it will automatically detect the addition of :comp2)
  GModelFit.update!(meval)
  println(maximum(GModelFit.last_evaluation(meval)))
  ```

- [`GModelFit.Residuals`](@ref): container for a `ModelEval` object, a `Measures` object, a `Vector{Float64}` to store the normalized residuals, and a minimizer instance.  A `Residuals` object contains all the relevant informations to perform minimization, and is therefore the only argument required for the [`GModelFit.minimize!`](@ref) function.  Since `Residuals` wraps a `ModelEval` object it is also able to detect changes in the original model. 

  An example of its usage is as follows:
  ```@example abc
  data = GModelFit.mock(Measures, model, dom)
  resid = GModelFit.Residuals(meval, data, GModelFit.lsqfit())
  GModelFit.minimize!(resid)
  println() # hide
  ```

  The [`GModelFit.MultiResiduals`](@ref) has the same purpose in the [Multi-dataset fitting](@ref) case. 
