# Basic concepts and data types

In order to exploit the **GFit.jl** model expressiveness we need to introduce a few concepts and the associated data types:

- *Domain*: an N-dimensional grid of points where the model is to be evaluated, it is analogous to the independent varible $\vec{x}$ in the $f(\vec{x})$ notation. It is represented by either:
  - a [`Domain{N}`](@ref) object for linear domains, i.e. where the coordinates for each of the N dimensions are explicitly specified for all the points;
  - or a [`CartesianDomain{N}`](@ref) object where the coordinates are specified for each of the `N` axis and the coordinates for all points are obtained as the cartesian product of the axis.  A carteisan domain is internally transformed into a linear one when needed;
  A domain object (either linear or cartesian) is required as first argument for the `Model` and `Measures` constructors (see below). However only the former is actually used during evaluation, while the latter is only used for visualization and consistency checks.

- *Measures*:  a container for the N-dimensional empirical data and their associated $1\sigma$ Gaussian uncertainties, represented by an object of type [`Measures{N}`](@ref) (further options may be available in the future, such as Poisson counts);

- *Model component*: the atomic building block of a (potentially very complex) model, it is essentially a function used to map a `Domain` or `CartesianDomain` object into a `Vector{Float64}` represeting the component evaluation.  All components are structures inheriting from `GFit.AbstractComponent` and are typically parametrized (see below).  The **GFit.jl** package provides several [Built-in components](@ref), and new ones can be implemented by the user.  The caching mechanism operates at the component level and aims to avoid unnecessary re-evaluation of the component if none has of its parameter values has changed since last evaluation;

- *Parameter*: a floating point value characterizing a specific aspect for the evaluation of a model component, e.g. the slope of a power law or the width of a Gaussian profile.  The parameter values are changed during the fitting process until the residuals between the global model evaluation and the empirical data are minimized.  A parameter can be fixed to a specific value, limited in an interval, and/or be dynamically calculated (patched) according to the values of other parameters.  All parameters are represented by an object of type [`GFit.Parameter`](@ref);

- *Model*: is the overall model description, whose evaluation is supposed to be compared to a single `Measures` objects and whose parameters are varied during fitting to reduce the residuals. All models are represented by an object of type [`Model`](@ref) containing a single `Domain` or `CartesianDomain` object representing the domain where the model will be evaluated, and one or more *components* characterizing the model itself.  Each component is identified by a unique name (actually a `Symbol`) within a model.
  - *Main component* of a model: the evaluation of a component, say `A`, may use the output of another component, say `B`, to calculate its output.  In this case we say tha `A` *depends* on `B`, and therefore `B` needs to be evaluated before `A` (circular dependencies are not allowed, and would raise an error if attempted).  The dependencies are automatically identified, and the last component being evaluated is dubbed *main component* since it is its output which represent the overall model evaluation;

- *Multi-model*: a container for two or more models, suitable to be compared to a corresponding number of `Measures` objects to perform multi-dataset fitting.  All models are identified by a unique integer identifier, starting from 1.  A multi-model is represented by an object of type `MultiModel`;


**λ-function**: TODO LComp, patch


- **Reducer**: a Julia function (encapsulated in a `Reducer` object) used to combine several component evaluations into a single evaluation output. If not explicitly mentioned when creating a `Model` object (see below), a default `Reducer` is created which simply performs an element-wise sum of all the components;



It is important to keep in mind the relationships among the above concepts:

```
Multi model (only if multiple models are involved)
 |
 + -- Model 1 (includes at least a reducer to combine the components into a model evaluation)
 |     |
 |     + -- Component 1
 |     |     |
 |     |     + -- Param. 1
 |     |     + -- Param. 2
 |     |     + -- ...
 |     + -- Component 2
 |     + -- etc.
 + -- Model 2
 + -- ...
```



The package most important function is `fit!`, whose purpose is to identify the best-fit parameter values which minimize the differences between the model evaluation and the empirical data.  The function arguments are a `Model` and a `Measures` objects, representing the fitting model and the empirical data respectively.  In the multi-dataset case the function accepts a `MultiModel` object and a vector of `Measures` objects.

The `fit!` function modifies its `Model` or `MultiModel` inputs (hence the exclamation mark in the name) by replacing the the initial parameter values with the best fit ones. The function returns an object of type `FitResult`, containing the fit statistics.
