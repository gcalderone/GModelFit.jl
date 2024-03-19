# API

## Index
```@index
```

## Exported symbols
The list of **GModelFit.jl** exported symbols is as follows:

```@docs
CartesianDomain{N}
Domain{N}
Measures{N}
Model
@fd
axis
compare
comptype
coords
domain
fit
fit!
freeze!
getindex
haskey
isfreezed
length
select_maincomp!
thaw!
uncerts
values
```


## Non-exported symbols
The following symbols are not exported by the **GModelFit.jl** package since they are typically not used in every day work, or aimed to debugging purposes.  Still, they can be useful in some case, hence they are documented here.

```@docs
GModelFit.evalcounter
GModelFit.evalcounters
GModelFit.CompEval
GModelFit.FitStats
GModelFit.FunctDesc
GModelFit.ModelEval
GModelFit.ModelSnapshot
GModelFit.MultiResiduals
GModelFit.Parameter
GModelFit.Residuals
GModelFit.comptypes
GModelFit.dependencies
GModelFit.evaluate!
GModelFit.last_evaluation
GModelFit.minimize!
GModelFit.mock
GModelFit.prepare!
GModelFit.serialize
GModelFit.update!
```
