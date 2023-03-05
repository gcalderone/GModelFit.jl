```@setup abc
include("setup.jl")
```

# Miscellaneous

## Generate mock datasets

In some case it is useful to test a model for robustness before the emprical data are available for fitting.  This can be achieved via the [`GModelFit.mock()`](@ref) function, whose purpose is to generate a mock dataset which simulates a measurement process by adding random noise to the foreseen ground-truth.

### Example

```@example abc
using GModelFit

dom = Domain(1:0.1:50)
model = Model(dom, :main => @λ (x, T=3.14) -> sin.(x ./ T) ./ (x ./ T))

# Generate a mock dataset
data = GModelFit.mock(Measures, model, seed=1)

# Fit model against the mock dataset
best, fitstats = fit(model, data)
```



## Serialization

A few structures, namely [`GModelFit.ModelSnapshot`](@ref), [`GModelFit.FitStats`](@ref) and [`Measures{N}`](@ref), as well as `Vector`(s) of such structures can be *serialized*, i.e. stored in a file using a dedicated JSON format.  The structures can lated be *de-serialized* in a separata Julia session without the need to re-run the fitting process used to create them in the first place.

### Example

In the following we will generate a few **GModelFit.jl** objects and serialized them in a file.
```@example abc
using GModelFit

dom = Domain(1:0.1:50)
model = Model(dom, :main => @λ (x, T=3.14) -> sin.(x ./ T) ./ (x ./ T))
data = GModelFit.mock(Measures, model, seed=1)
best, fitstats = fit(model, data)

# Serialize objects and save in a file
GModelFit.serialize("save_for_future_use.json", best, fitstats, data)
println(); # hide
```

The same objects can be de-serialized in a different Julia session:
```@example abc
using GModelFit
best, fitstats, data = GModelFit.deserialize("save_for_future_use.json")
```



## Quick plot (1D)

The **GModelFit.jl** package implements [**Gnuplot.jl**](https://github.com/gcalderone/Gnuplot.jl/) recipes to display plots of `Measures{1}` and `ModelSnapshot` objects., e.g.:

### Example

Create a model, a mock dataset and run a fit:
```@example abc
using GModelFit

dom = Domain(0:0.01:5)
model = Model(dom, :bkg => GModelFit.OffsetSlope(1, 1, 0.1),
                   :l1 => GModelFit.Gaussian(1, 2, 0.2),
                   :l2 => GModelFit.Gaussian(1, 3, 0.4),
                   :main => SumReducer(:bkg, :l1, :l2))
data = GModelFit.mock(Measures, model)
best, res = fit(model, data)
println(); # hide
```

A plot of the dataset and of the best fit model can be simply obtained with
```@example abc
using Gnuplot
@gp data best
saveas("gnuplot1") # hide
```
![](assets/gnuplot1.png)

You may also specify axis range, labels, title, etc. using the standard [**Gnuplot.jl**](https://github.com/gcalderone/Gnuplot.jl/) keyword syntax, e.g.:

```@example abc
using Gnuplot
@gp xr=[1, 4.5] xlabel="Wavelength" ylab="Flux" "set key outside" data best
saveas("gnuplot2") # hide
```
![](assets/gnuplot2.png)
