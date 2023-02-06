```@setup abc
using Gnuplot
Gnuplot.quitall()
mkpath("assets")
Gnuplot.options.term = "unknown"
empty!(Gnuplot.options.init)
push!( Gnuplot.options.init, linetypes(:Set1_5, lw=2.5, ps=1.5))
saveas(file) = save(term="pngcairo size 550,350 fontscale 0.8", output="assets/$(file).png")
dumpjson(file, args...) = GFit.snapshot_json("assets/$(file).json", args...)
```


# Multi-dataset fitting

**GFit.jl** is able to simultaneously fit several models against a corresponding number of datasets, while placing constraints among the models. Typical use cases are:
- a single phenomenon is observed with two (or more) instruments/detectors;
- a single phenomenon is observed at different times;

Fitting multiple datasets simultaneously may provide tighter constraints on the best fit parameters under the assumption that the models are somehow related, i.e. that their parameters are constrained (or *patched*).

To perform a multi-dataset fitting we should create one `Model` for each dataset in the usual way, collect them in a `MultiModel` object, and define patch constraints among models.   The following example shows how to fit two Gaussian curves under the hypotesis that the center and normalization parameters are the same:
```@example abc
using GFit

# Create domain, individual models and the MultiModel container
dom = Domain(-5.:5)
model1 = Model(dom, GFit.Gaussian(1, 0., 1.))
model2 = Model(dom, GFit.Gaussian(1, 0., 1.))
multi = MultiModel(model1, model2)

# Patch parameters
multi[2][:main].norm.mpatch   = @λ m -> m[1][:main].norm
multi[2][:main].center.mpatch = @λ m -> m[1][:main].center

# Create datasets and fit
data1 = Measures(dom, [-0.006,  0.015,  0.001,  0.049,  0.198,  0.430,  0.226,  0.048,  0.017, -0.001, -0.006], 0.04)
data2 = Measures(dom, [-0.072, -0.033, -0.070,  0.108,  0.168,  0.765,  0.113, -0.054,  0.032,  0.013,  0.015], 0.04)
res = fit!(multi, [data1, data2])
dumpjson("ex_multifit", multi, data1, data2, res) # hide
show(res) # hide
```

The `MultiModel` object is very similar to a vector, hence the inidividual models can be retrieved by a numeric index, as in `multi[2][:main].norm`.  Also, the bestfit values can be retrieved by indexing the `bestfit` field of the results, e.g.:
```@example abc
println("Width of Gaussian 1: ", res.bestfit[1][:main].sigma.val, " ± ", res.bestfit[1][:main].sigma.unc, "\n")
println("Width of Gaussian 2: ", res.bestfit[2][:main].sigma.val, " ± ", res.bestfit[2][:main].sigma.unc, "\n")
println("Reduced χ^2: ", res.fitstat)
```
