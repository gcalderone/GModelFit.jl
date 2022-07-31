# GFit.jl

**Note: This package is now mature enough to be released as v0.1, however it is not yet a registered package since it lacks the associated documentation.**

[![License](http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat)](LICENSE.md)
[![DocumentationStatus](https://img.shields.io/badge/docs-stable-blue.svg?style=flat)](https://gcalderone.github.io/GFit.jl/v0.1.0/index.html)

`GFit` is a general purpose, data-driven model fitting framework for Julia.
See below for a simple example, and [here](https://gcalderone.github.io/GFit.jl/v0.1.0/index.html) for extensive documentation.

## Installation

Install with:
```julia
]add GFit
```

## Simple example

```julia
using GFit

# Prepare vectors with domain points, empirical measures an their uncertainties and
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

You may also plot the empirical data and the best fit model with a plotting framework of your choice, e.g. [Gnuplot.jl](https://github.com/gcalderone/Gnuplot.jl):
```julia
using Gnuplot
@gp coords(dom) values(data) uncerts(data) "w yerr t 'Data'"
@gp :- coords(domain(model)) model() "w l t 'Best fit model'"
```
![Example plot](https://github.com/gcalderone/GFit.jl/blob/master/examples/ex0.png)
