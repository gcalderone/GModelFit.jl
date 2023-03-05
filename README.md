# GFit.jl

[![License](http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat)](LICENSE.md)
[![DocumentationStatus](https://img.shields.io/badge/docs-stable-blue.svg?style=flat)](https://gcalderone.github.io/GFit.jl/)

`GFit` is a general purpose, data-driven model fitting framework for Julia.

## Installation

Install with:
```julia
]add GFit
```

## Example

```julia
using GFit

# Prepare vectors with domain points, empirical measures and uncertainties
x    = [0.1, 1.1, 2.1, 3.1, 4.1]
meas = [6.29, 7.27, 10.41, 18.67, 25.3]
unc  = [1.1, 1.1, 1.1, 1.2, 1.2]
dom  = Domain(x)
data = Measures(dom, meas, unc)

# Create a model using an explicit mathematical expression, and provide the
# initial guess values:
model = Model(dom, @λ (x, a2=1, a1=1, a0=5) -> (a2 .* x.^2  .+  a1 .* x  .+  a0))

# Fit model to the data
best, res = fit(model, data)
```

The output is as follows:
```julia
(Components:
╭───────────┬────────────┬─────────────┬───────────┬───────────┬───────────┬─────────╮
│ Component │ Type       │ Eval. count │ Min       │ Max       │ Mean      │ NaN/Inf │
├───────────┼────────────┼─────────────┼───────────┼───────────┼───────────┼─────────┤
│ main      │ GFit.FComp │ 76          │     6.088 │     25.84 │     13.56 │ 0       │
╰───────────┴────────────┴─────────────┴───────────┴───────────┴───────────┴─────────╯

Parameters:
╭───────────┬────────┬───────┬───────────┬───────────┬────────┬───────╮
│ Component │ Param. │ Range │ Value     │ Uncert.   │ Actual │ Patch │
├───────────┼────────┼───────┼───────────┼───────────┼────────┼───────┤
│ main      │ a2     │       │     1.201 │    0.3051 │        │       │
│           │ a1     │       │    -0.106 │     1.317 │        │       │
│           │ a0     │       │     6.087 │     1.142 │        │       │
╰───────────┴────────┴───────┴───────────┴───────────┴────────┴───────╯
, Fit results:
    #Data :        5         Red. fit stat.:     1.0129  (DOF: 2)
    #Free :        3         Elapsed time  :      0.004
    Status:       OK
```
