```@setup abc
include("setup.jl")
```

# Basic concepts and data types
In order to exploit the **GModelFit.jl** model expressiveness a few concepts need to be introduced, along with their associated data types:

- *Domain*: an N-dimensional grid of points associated to empirical measures, and used to evaluate a model.  It is analogous to the independent varible $\vec{x}$ in the $f(\vec{x})$ notation. It is represented by either:
  - a [`Domain{N}`](@ref) object for linear domains, where the coordinates for each of the N dimensions are explicitly specified for all the points;
  - or a [`CartesianDomain{N}`](@ref) object where the coordinates are specified for each of the `N` axis and the coordinates for all points are obtained as the cartesian product of all the axes.  A cartesian domain is internally transformed into a linear one when needed;
  A domain object (either linear or cartesian) is required as first argument for the `Measures` constructor (see below).

- *Measures*: a container for the N-dimensional empirical data and their associated $1\sigma$ Gaussian uncertainties, represented by an object of type [`Measures{N}`](@ref) (further options may be available in the future, such as Poisson counts);

- *Model component*: the atomic building block of a (potentially complex) model, it is essentially a function used to map a `Domain` or `CartesianDomain` object into a `Vector{Float64}` representing the component evaluation.  A component is a structure inheriting from `GModelFit.AbstractComponent` and is typically characterized by one or more *parameters* (see below).   The **GModelFit.jl** package provides several [Built-in components](@ref), and new ones can be implemented by the user (see [Custom components](@ref)).  The memoization mechanism operates at the component level and aims to avoid unnecessary re-evaluation of the component if none of its parameter values has changed since last evaluation;

- *Parameter*: a single floating point number characterizing a specific aspect for the evaluation of a component (e.g. the slope of a power law or the width of a Gaussian profile). The parameter values are automatically varied during the fitting process until the residuals between the global model evaluation and the empirical data are minimized.  A parameter can be fixed to a specific value, limited in an interval, and/or be dynamically calculated (patched) according to the values of other parameters.  All parameters are represented by an object of type [`GModelFit.Parameter`](@ref);

- *Model*: is the overall model description, whose evaluation is supposed to be compared to a single `Measures` object and whose parameters are varied during fitting to reduce the residuals.  Internally, a model is implemented as a dictionary containing one or more *components*, each identified by a unique `Symbol` name (see [`Model`](@ref));

  - Component dependencies and *main component*: the evaluation of a component, say `A`, may use the outcome of another component, say `B`, to calculate its output, thus inducing a dependency between the two. In this case we say that `A` *depends* on `B`, and therefore `B` needs to be evaluated before `A` (circular dependencies are not allowed, and would raise an error if attempted).  The dependencies are automatically identified, and the last component being evaluated is dubbed *main component* since its output represent the overall model evaluation;

- *Multi-model*: a `Vector{Model}` containing two or more models, suitable to be compared to a corresponding `Vector{Measures}` to perform [Multi-dataset fitting](@ref);

- *Minimizer*: the **GModelFit.jl** package provides just the tools to define and manipulate a model, but the actual fitting (namely, the minimization of the residuals) is performed by an external *minimizer* library.  Two minimizers are currently available:
  - [LsqFit](https://github.com/JuliaNLSolvers/LsqFit.jl): a pure-Julia minimizer;
  - [CMPFit](https://github.com/gcalderone/CMPFit.jl): a C minimizer wrapped in a Julia package.
  Both are automatically installed with **GModelFit.jl**, and `LsqFit` is the default choice (unless otherwise specified in the [`fit()`](@ref) function call).

- *Model snapshot*: the best fit model, as well as the best fit parameter values and associated uncertainties, are returned by the [`fit()`](@ref) function as a [`GModelFit.ModelSnapshot`](@ref) structure, namely a *frozen snapshot* of the evaluation of a `Model` object on a given `Domain`.  Components, parameters and evaluations outcomes are accessed in exactly the same way on both `Model` and `ModelSnapshot` objects, the only difference being that the latter can't be re-evaluated on different domains or parameter values.

- *Fit statistics*: the purpose of fitting is to minimize the *distance* between the model and the data, as quantified by a proper fit statistic (typically a reduced $\chi^2$ for the Gaussian uncertainties case). Such statistic, as well as other information concerning the fit, are returned by the [`fit()`](@ref) function in a [`GModelFit.FitSummary`](@ref) structure;

- *function descriptor*: **GModelFit.jl** uses standard Julia function in two different contexts:
  - to calculate the value of a `Parameter` as a function of other `Parameter`'s values. In this case the parameters are said to be *patched*, or linked, since there is a constraint between their values.  Two (or more) parameters may be patched within the same model, or across models when performing [Multi-dataset fitting](@ref);
  - to define a model component using a standard Julia mathematical expression involving `Parameter`s values or other components;

  To use a standard function in this fashion it should be wrapped into a [`GModelFit.FunctDesc`](@ref) object which allows both to invoke the function itself, as well as to provide a string representation for display purposes.  In order to create a function descriptor object it typically is much easier to invoke the [`@fd`](@ref) macro rather than the `FunctDesc` constructor.
