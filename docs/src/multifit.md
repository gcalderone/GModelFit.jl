```@setup abc
include("setup.jl")
```


# Multi-dataset fitting

**GModelFit.jl** is able to simultaneously fit several models against a corresponding number of datasets, while placing constraints among the models. Typical use cases are:
- a single phenomenon is observed with two (or more) instruments/detectors;
- a single phenomenon is observed at different times;

Fitting multiple datasets simultaneously may provide tighter constraints on the best fit parameters under the assumption that the models are somehow related, i.e. that their parameters are constrained (or *patched*).

To perform a multi-dataset fitting we should create one `Model` for each dataset in the usual way, collect them in a `Vector{Model}`, and define patch constraints among models.   The following example shows how to fit two Gaussian curves under the hypotesis that the center and normalization parameters are the same:
```@example abc
using GModelFit

# Create individual models and the Vector{Model} container
model1 = Model(GModelFit.Gaussian(1, 0., 1.))
model2 = Model(GModelFit.Gaussian(1, 0., 1.))
multi = [model1, model2]

# Patch parameters
multi[2][:main].norm.mpatch   = @fd m -> m[1][:main].norm
multi[2][:main].center.mpatch = @fd m -> m[1][:main].center

# Create datasets and fit
dom = Domain(-5.:5)
data1 = Measures(dom, [-0.006,  0.015,  0.001,  0.049,  0.198,  0.430,  0.226,  0.048,  0.017, -0.001, -0.006], 0.04)
data2 = Measures(dom, [-0.072, -0.033, -0.070,  0.108,  0.168,  0.765,  0.113, -0.054,  0.032,  0.013,  0.015], 0.04)
bestfit, fsumm = fit(multi, [data1, data2])
show((bestfit, fsumm)) # hide
```

The best fit models and values are returned as a `Vector{ModelSnapshot}` in `bestfit`, i.e.:
```@example abc
println("Width of Gaussian 1: ", bestfit[1][:main].sigma.val, " ± ", bestfit[1][:main].sigma.unc)
println("Width of Gaussian 2: ", bestfit[2][:main].sigma.val, " ± ", bestfit[2][:main].sigma.unc)
println("Reduced χ^2: ", fsumm.fitstat)
```
