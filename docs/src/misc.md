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
