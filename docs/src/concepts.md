# Basic concepts and data types

In order to exploit the **GFit.jl** model expressiveness we need to introduce a few concepts, and the associated data types:

- *Domain*: an N-dimensional grid of points where the model is to be evaluated, it is analogous to the independent varible $\vec{x}$ in the $f(\vec{x})$ notation. It is represented by either:
  - a [`Domain{N}`](@ref) object for linear domains, i.e. where the coordinates for each of the N dimensions are explicitly specified for all the points;
  - or a [`CartesianDomain{N}`](@ref) object where the coordinates are specified for each of the `N` axis and the coordinates for all points are obtained as the cartesian product of all the axes.  A cartesian domain is internally transformed into a linear one when needed;
  A domain object (either linear or cartesian) is required as first argument for the `Model` and `Measures` constructors (see below).

- *Measures*: a container for the N-dimensional empirical data and their associated $1\sigma$ Gaussian uncertainties, represented by an object of type [`Measures{N}`](@ref) (further options may be available in the future, such as Poisson counts);

- *Model component*: the atomic building block of a (potentially very complex) model, it is essentially a function used to map a `Domain` or `CartesianDomain` object into a `Vector{Float64}` representing the component evaluation.  A component is a structure inheriting from `GFit.AbstractComponent` and is typically characterized by one or more *parameters* (see below). One component's evaluation can also be used as input for another one's calculation, thus inducing a dependency between the two.  The **GFit.jl** package provides several [Built-in components](@ref), and new ones can be implemented by the user.  The memoization mechanism operates at the component level and aims to avoid unnecessary re-evaluation of the component if none of its parameter values has changed since last evaluation;

- *Parameter*: a single floating point number characterizing a specific aspect for the evaluation of a component (e.g. the slope of a power law or the width of a Gaussian profile). The parameter values are automatically varied during the fitting process until the residuals between the global model evaluation and the empirical data are minimized.  A parameter can be fixed to a specific value, limited in an interval, and/or be dynamically calculated (patched) according to the values of other parameters.  All parameters are represented by an object of type [`GFit.Parameter`](@ref);

- *Model*: is the overall model description, whose evaluation is supposed to be compared to a single `Measures` objects and whose parameters are varied during fitting to reduce the residuals. All models are represented by an object of type [`Model`](@ref) containing a single `Domain` or `CartesianDomain` object representing the domain where the model will be evaluated, and one or more *components* characterizing the model itself.  Each component is identified by a unique name (actually a `Symbol`) within a model.
  - Component dependencies and *main component*: the evaluation of a component, say `A`, may use the output of another component, say `B`, to calculate its output.  In this case we say that `A` *depends* on `B`, and therefore `B` needs to be evaluated before `A` (circular dependencies are not allowed, and would raise an error if attempted).  The dependencies are automatically identified, and the last component being evaluated is dubbed *main component* since its output represent the overall model evaluation;

- *Multi-model*: a container for two or more models, suitable to be compared to a corresponding number of `Measures` objects to perform [Multi-dataset fitting](@ref).  All models are identified by a unique integer identifier, starting from 1.  A multi-model is represented by an object of type `MultiModel`;

- *Fit results*: the purpose of fitting is to minimize the *distance* between the model and the data, as quantified by a proper fit statistic (typically a reduced $\chi^2$ for the Gaussian uncertainties case). Such statistic, as well as other information concerning the fit and the best fit parameter values and uncertainties, are returned by the [`fit!()`](@ref) function in a [`GFit.FitResult`](@ref) structure.

- *λ-function* ([`GFit.λFunct`](@ref)): is an anonymous function used in two different contexts within **GFit.jl**:
  - to calculate the value of a `Parameter` as a function of other `Parameter`'s values. In this case the parameters are said to be *patched*, or linked, since there is a constraint between their values.  Two (or more) parameters may be patched within the same model, or across models in a multi-model analysis;
  - to define a model component using a standard Julia mathematical expression involving `Parameter`s values or other components.
  In both cases the λ-function is generated using the [`@λ`](@ref) macro and the standard Julia syntax for anonymous functions.

- *Minimizer*: the **GFit.jl** package provides just the tools to define and manipulate a model, but the actual fitting (namely, the minimization of the residuals) is performed by an external *minimizer* library.  Two minimizers are currently available:
  - [LsqFit](https://github.com/JuliaNLSolvers/LsqFit.jl): a pure-Julia minimizer;
  - [CMPFit](https://github.com/gcalderone/CMPFit.jl): a C minimizer wrapped in a Julia package.
  Both are automatically installed with **GFit.jl**, and `LsqFit` is the default choice (unless otherwise specified in the [`fit!()`](@ref) function call).  However, for some of the most complex cases `CMPFit` seems to be more robust and less sensitive to initial guess parameters.

- *Mock data*: evaluating a a model may be useful even before actual data are available, e.g. to test its robustness and capabilities.  To this purpose **GFit.jl** provides the [`GFit.mock()`](@ref) function which is able to generate mock data set(s) using a (multi-)model as ground truth, and add a random noise to simulate the measurement process. This functionality is used in all the examples presented in the next sections.


## How to access the data structures

Many of the above mentioned data structures are accessible using indexing (as in dictionary or vectors) or a `struct`-like interface, starting from a single [`Model`](@ref) or a [`MultiModel`](@ref) object.  Hence it is important to keep in mind the relationships among these objects to be able to access them:
```
Multi-model (`multi`)
 |
 + -- Model 1 (`model`)
 |     |
 |     + -- Component1
 |     |     |
 |     |     + -- Param1
 |     |     + -- Param2
 |     |     + -- ...
 |     + -- Component2
 |     + -- etc.
 + -- Model 2
 + -- ...
```

E.g. the syntax to access the value of a parameter in a single model case is: `model[:Component1].Param1.val`.  In a multi-model case it is: `multi[1][:Component1].Param1.val`.

