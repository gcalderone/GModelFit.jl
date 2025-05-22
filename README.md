# GModelFit.jl

[![License](http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat)](LICENSE.md)
[![DocumentationStatus](https://img.shields.io/badge/docs-stable-blue.svg?style=flat)](https://gcalderone.github.io/GModelFit.jl/)

`GModelFit` is a general purpose, data-driven model fitting framework for Julia.

> [!WARNING]
> The code in version 0.5.0 underwent a signficant refactoring, and a few details may break your code.
> **Please have a look at ChangeLog.md !!**

## Installation

Install with:
```julia
]add GModelFit
```

## Example

```julia
using GModelFit

# Prepare vectors with domain points, empirical measures and uncertainties
x    = [0.1, 1.1, 2.1, 3.1, 4.1]
meas = [6.29, 7.27, 10.41, 18.67, 25.3]
unc  = [1.1, 1.1, 1.1, 1.2, 1.2]
dom  = Domain(x)
data = Measures(dom, meas, unc)

# Create a model using an explicit mathematical expression, and provide the
# initial guess values:
model = Model(@fd (x, a2=1, a1=1, a0=5) -> (a2 .* x.^2  .+  a1 .* x  .+  a0))

# Fit model to the data
bestfit, fsumm = fit(model, data)
```

The output is as follows:
```julia
(Components:
╭───────────┬───────┬───────┬─────────────┬───────────┬───────────┬───────────┬─────────╮
│ Component │ Type  │ #Free │ Eval. count │ Min       │ Max       │ Mean      │ NaN/Inf │
├───────────┼───────┼───────┼─────────────┼───────────┼───────────┼───────────┼─────────┤
│ main      │ FComp │ 3     │ 66          │     6.088 │     25.84 │     13.56 │ 0       │
╰───────────┴───────┴───────┴─────────────┴───────────┴───────────┴───────────┴─────────╯

Parameters:
╭───────────┬───────┬────────┬──────────┬───────────┬───────────┬────────┬───────╮
│ Component │ Type  │ Param. │ Range    │ Value     │ Uncert.   │ Actual │ Patch │
├───────────┼───────┼────────┼──────────┼───────────┼───────────┼────────┼───────┤
│ main      │ FComp │ a2     │ -Inf:Inf │     1.201 │    0.3051 │        │       │
│           │       │ a1     │ -Inf:Inf │    -0.106 │     1.317 │        │       │
│           │       │ a0     │ -Inf:Inf │     6.087 │     1.142 │        │       │
╰───────────┴───────┴────────┴──────────┴───────────┴───────────┴────────┴───────╯
, Fit summary: #data: 5, #free pars: 3, red. fit stat.: 1.0129, status: OK
)
```
