```@setup abc
include("setup.jl")
```

# Miscellaneous

## Generate mock datasets

In some case it is useful to test a model for robustness before the emprical data are available for fitting.  This can be achieved via the [`GFit.mock()`](@ref) function, whose purpose is to generate a mock dataset which simulates a measurement process by adding random noise to the foreseen ground-truth.

### Example

```@example abc
using GFit

dom = Domain(1:0.1:50)
model = Model(dom, :main => @λ (x, T=3.14) -> sin.(x ./ T) ./ (x ./ T))

# Generate a mock dataset
data = GFit.mock(Measures, model, seed=1)

# Fit model against the mock dataset
best, fitstats = fit(model, data)
```



## Serialization

A few structures, namely [`GFit.ModelSnapshot`](@ref), [`GFit.FitStats`](@ref) and [`Measures{N}`](@ref), as well as `Vector`(s) of such structures can be *serialized*, i.e. stored in a file using a dedicated JSON format.  The structures can lated be *de-serialized* in a separata Julia session without the need to re-run the fitting process used to create them in the first place. 

### Example

In the following we will generate a few **GFit.jl** objects and serialized them in a file.
```@example abc
using GFit

dom = Domain(1:0.1:50)
model = Model(dom, :main => @λ (x, T=3.14) -> sin.(x ./ T) ./ (x ./ T))
data = GFit.mock(Measures, model, seed=1)
best, fitstats = fit(model, data)
  
# Serialize objects and save in a file
GFit.serialize("save_for_future_use.json", best, fitstats, data)
println(); # hide
```

The same objects can be de-serialized in a different Julia session:
```@example abc
using GFit
best, fitstats, data = GFit.deserialize("save_for_future_use.json")
```


## Viewers (1D case)

The [**GFitViewer.jl**](https://github.com/lnicastro/GFitViewer.jl) package provides a few tools to quickly inspect **GFit.jl** objects.  Specifically, it provides:
- [**Gnuplot.jl**](https://github.com/gcalderone/Gnuplot.jl/) recipes to display `Measures{1}`, and `ModelSnapshot` objects;
- The `viewer()` function to generate an HTML files showing a plot of the model and data, as well as the same text output you would see in a Julia REPL session.

### Example

In the following we will create a model, generate a mock dataset, and fit:
```@example abc
dom = Domain(0:0.01:5)
model = Model(dom, :bkg => GFit.OffsetSlope(1, 1, 0.1),
                   :l1 => GFit.Gaussian(1, 2, 0.2),
                   :l2 => GFit.Gaussian(1, 3, 0.4),
                   :main => SumReducer(:bkg, :l1, :l2))
model[:bkg].offset.val = 1
model[:bkg].offset.fixed = true
model[:bkg].slope.low  = 0
model[:bkg].slope.high = 0.2
model[:l2].norm.patch = :l1
model[:l2].sigma.patch = @λ m -> 2 * m[:l1].sigma
model[:l2].center.patch = @λ (m, v) -> v + m[:l1].center
model[:l2].center.val = 1   # guess value for the distance between the centers
model[:l2].center.low = 0   # ensure [l2].center > [l1].center
data = GFit.mock(Measures, model)
best, res = fit(model, data)
println(); #hide
```

The dataset and best fit model can be easily inspected with:
```@example abc
using Gnuplot, GFitViewer
@gp data best
saveas("example_viewers")
println(); # hide
```
![](assets/example_viewers.png)


Also, you can generate an HTML file (which will be automatically opened using your default browser) with:

```@example abc
using GFitViewer
viewer(best, fitstats, data);
```
