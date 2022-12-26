# GFit.jl

**Note: This package is now mature enough to be released as v0.1, however it is not yet a registered package since it lacks the associated documentation.**

[![License](http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat)](LICENSE.md)
[![DocumentationStatus](https://img.shields.io/badge/docs-stable-blue.svg?style=flat)](https://gcalderone.github.io/GFit.jl/)

`GFit` is a general purpose, data-driven model fitting framework for Julia.
See below for a few examples, and [here](https://gcalderone.github.io/GFit.jl/v0.1.0/index.html) for extensive documentation.

## Installation

Install with:
```julia
]add GFit
```

## Examples

### Fit using an analytical formula
```julia
using GFit

# Prepare vectors with domain points, empirical measures and uncertainties, and
# use them to create the GFit input objects:
x    = [0.1, 1.1, 2.1, 3.1, 4.1]
meas = [6.29, 7.27, 10.41, 18.67, 25.3]
unc  = [1.1, 1.1, 1.1, 1.2, 1.2]
dom  = Domain(x)
data = Measures(dom, meas, unc)

# Create a model using an explicit mathematical expression, and provide the
# initial guess values:
model = Model(dom, @λ (x, a2=1, a1=1, a0=5) -> (a2 .* x.^2  .+  a1 .* x  .+  a0))

# Fit model to the data
res = fit!(model, data)
```

The output is as follows:
```julia
Best fit parameters:
╭───────────┬────────┬───────┬───────────┬───────────┬────────┬───────╮
│ Component │ Param. │ Range │ Value     │ Uncert.   │ Actual │ Patch │
├───────────┼────────┼───────┼───────────┼───────────┼────────┼───────┤
│ main      │ a2     │       │     1.201 │    0.3051 │        │       │
│           │ a1     │       │    -0.106 │     1.317 │        │       │
│           │ a0     │       │     6.087 │     1.142 │        │       │
╰───────────┴────────┴───────┴───────────┴───────────┴────────┴───────╯
Fit results:
    #Data :        5              #Free params  :          3
    DOF   :        2              Red. fit stat.:     1.0129
    Status:       OK              Elapsed time  :          0 s
```

You may plot the data and the best fit model with a plotting framework of your choice. E.g., with [Gnuplot.jl](https://github.com/gcalderone/Gnuplot.jl):
```julia
using Gnuplot
@gp coords(dom) values(data) uncerts(data) "w yerr t 'Data'" :-
@gp :- coords(dom) model() "w l t 'Best fit model'"
```
![Example plot](https://github.com/gcalderone/GFit.jl/blob/master/examples/ex0.png)



### Fit using built-in components

```julia
# Create a GFit model
x = 0:0.05:6
model = Model(Domain(x),
              :l1  => GFit.Gaussian(1, 2, 0.2),
              :l2  => GFit.Gaussian(1, 3, 0.4),
              :bkg => GFit.OffsetSlope(0.5, 1, 0.1),
              :main => SumReducer(:l1, :l2, :bkg))

# Constrain normalization of :l2 to be the same as :l1
model[:l2].norm.patch = :l1

# Constrain width of :l2 to be twice that of :l1
model[:l2].sigma.patch = @λ m -> 2 * m[:l1].sigma

# Generate a mock dataset with random noise
data = GFit.mock(Measures, model)

# Fit model to the data
res = fit!(model, data)
```
The output is as follows:
```julia
Best fit parameters:
╭───────────┬───────────────┬───────┬───────────┬───────────┬───────────┬───────────────────────╮
│ Component │ Param.        │ Range │ Value     │ Uncert.   │ Actual    │ Patch                 │
├───────────┼───────────────┼───────┼───────────┼───────────┼───────────┼───────────────────────┤
│ l1        │ norm          │ 0:Inf │     1.012 │   0.02514 │           │                       │
│           │ center        │       │     2.011 │  0.006576 │           │                       │
│           │ sigma         │ 0:Inf │    0.1974 │  0.004909 │           │                       │
├───────────┼───────────────┼───────┼───────────┼───────────┼───────────┼───────────────────────┤
│ l2        │ norm (FIXED)  │ 0:Inf │         1 │       NaN │     1.012 │ l1                    │
│           │ center        │       │     2.996 │   0.01845 │           │                       │
│           │ sigma (FIXED) │ 0:Inf │       0.4 │       NaN │    0.3949 │ m->2 * (m[:l1]).sigma │
├───────────┼───────────────┼───────┼───────────┼───────────┼───────────┼───────────────────────┤
│ bkg       │ offset        │       │    0.5081 │   0.01831 │           │                       │
│           │ x0 (FIXED)    │       │         1 │       NaN │           │                       │
│           │ slope         │       │   0.09736 │  0.006024 │           │                       │
╰───────────┴───────────────┴───────┴───────────┴───────────┴───────────┴───────────────────────╯
Fit results:
    #Data :      121              #Free params  :          6
    DOF   :      115              Red. fit stat.:    0.88978
    Status:       OK              Elapsed time  :      0.003 s
```

Use [GfitViewer](https://github.com/lnicastro/GFitViewer.jl) to display the results:
```julia
using GFitViewer
viewer(model, data, res)  # opens an HTML viewer in yuour browser

using GFitViewer, Gnuplot
@gp data model  # uses Gnuplot.jl recipes to display the plot
```
