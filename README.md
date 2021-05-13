# GFit.jl

**Note: This package is still in active development, and the interface may undergo significant breaking changes.**

[![Build Status](https://travis-ci.org/gcalderone/GFit.jl.svg?branch=master)](https://travis-ci.org/gcalderone/GFit.jl)

`GFit` is a general purpose, data-driven model fitting framework for Julia.

The typical use case for is as follows: you observe a physical phenomenon with one (or more) instrument(s) and wish to fit those empirical data against a (possibly very complex) theoretical model, in order to extract the characterizing physical quantities represented by the model parameters.

`GFit` provides the basic tools to build, evaluate, and fit complex models against empirical data.  The main functionalities are:

- it handles datasets of any dimensionality;
- the fitting model is built by combining one or more *components* (either built-in or implemented by the user) using a standard Julia mathematical expression, and evaluated on a user provided domain;
- all components results are cached so that repeated evaluations with the same parameter values do not involve further calculations;
- user provided components can pre-compute quantities based on the model domain, and store them in reserved areas for re-use;
- model parameters can be fixed to a specific value, limited in an interval, and/or be dynamically calculated (patched) according to the values of other parameters;
- multiple data sets can be fitted simultaneously;
- it support different minimizers ([LsqFit](https://github.com/JuliaNLSolvers/LsqFit.jl) and [CMPFit](https://github.com/gcalderone/CMPFit.jl) );
- it provides several facilities for interactive fitting and result displaying.

See below for a simple example, and the `examples` directory for more complex ones.

## Installation

`GFit` is not a registered package, hence the only way to install it is as follows:

```julia
] add https://github.com/gcalderone/GFit.jl#master
```

## Basic concepts

`GFit`  flexibility arise from a clear identification of the basic concepts involved in the fitting process, and a careful design of their corresponding data structures. The following list provides a brief description of such concepts:

- ***Domain***: an N-dimensional domain where the model is supposed to be evaluated.  It is analogous to the independent varible $\vec{x}$ in the $f(\vec{x})$ notation, and is represented by an object of type `Domain{N}`, where `N` (>=1) is the dimensionality.  All numbers in a domain are treated as `Float64`;

- ***Measures***: a container for the N-dimensional empirical data and their $1\sigma$ Gaussian uncertainties, represented by an object of type `Measures{N}` (in the future further options will be available, such as Poisson counts);

- ***Component***: a Julia function used to map a `Domain` object into a `Vector{Float64}` , represeting the component evaluation.  Components are typically parametrized, and the values for such parametrized are changed during the fitting process until an optimal match between the global model and the empirical data is achieved.  All component objects have a common supertype (`GFit.AbstractComponent`), and represent the atomic *building block* of a (potentially very complex) model. The `GFit` caching mechanism operates at the component level (i.e. a component is actually evaluated only if at least one of its parameter values has changed);

- ***Reducer***: a Julia function (encapsulated in a `Reducer` object) used to combine several component evaluations into a single evaluation output. If not explicitly mentioned when creating a `Prediction` object (see below), a default `Reducer` is created which simply performs an element-wise sum of all the components;

- ***Prediction***: a container for the following objects:
  
  - a `Domain` object;
  
  - one or more components, all evaluated on the same domain. Each component is identified by unique name (actually a `Symbol`) within a prediction;
  
  - a `Reducer` object, used to combine the component evaluations into a single evaluation output;
  
  The prediction evaluation is suitable to be compared to the empirical data in a `Measures` object. All predictions are represented by object of type `Prediction`;

- ***Model***: a container for one or more predictions, representing the global model, and suitable to be compared to several `Measures` objects to perform multi-dataset fitting.  All predictions in a model are identified by a unique integer ID, starting from 1. If only one prediction is present, the concept of *model* and *prediction* are almost equivalent, and the `GFit` API implicitly assumes `id=1` in all methods accepting a `Model` object. The mechanism to link a parameter value to an expression involving other parameters operates at the `Model` level.

The package most important function is `fit!`. The purpose of the function is to identify the parameter values which minimize the differences between the model evaluation and the empirical data.  The function arguments are:

- a`Model` object (containing one or more `Prediction` objects;

- either a scalar or a vector of `Measures` objects, representing the empirical data.

The `fit!` function modifies its `Model` input (hence the exclamation mark in the name) by replacing the the initial parameter values with the best fit ones. The function returns an object of type `BesftFitResult`, containing the best fit parameter values, along with their uncertainties.

## Examples

All the following example follow the same patterns:

1. Prepare a `Domain` object and define a model (combining one or more components);

2. Set the initial values for the model parameter to the *true* ones;

3. Evaluate the model in order to obtain the ground truth (i.e. the ideal output of a measurement process), and add some noise to generate a simulated dataset;

4. Fit the model against the dataset, and display the results.

In a real life example, where the ground truth is unknown and we wish to estimate the parameter values, we would need to provide sensible inital guess values in step 2, and use actual empirical dataset(s) in step 3.

### Simple example using a mathematical expression

The model to be compared with empirical data is a simple analytical formula:

```julia
f(x, p1, p2, p3, p4, p5) = @. (p1  +  p2 * x  +  p3 * x^2  +
                               p4 * sin(p5 * x))  *  cos(x)
```

where `p1`, ..., `p5` are the model parameters to be estimated.  We'll assume the *true* values for the parameters are:

```julia
true_vals = [1, 2, 1.e-3, 4, 5]
```

The code is as follows:

```julia
using GFit, Random

# Model formula and true parameter values
f(x, p1, p2, p3, p4, p5) = @. (p1  +  p2 * x  +  p3 * x^2  +
                               p4 * sin(p5 * x))  *  cos(x)
true_vals = [1, 2, 1.e-3, 4, 5]

# Prepare domain
dom = Domain(1:1e-2:10)

# Prepare a model with a single component of type GFit.FuncWrap
model = Model(dom, :f => GFit.FuncWrap(f, true_vals...))

# Create mock data set by evaluating the model and adding some
# random noise.  All uncertainties are equal to 1.
rng = MersenneTwister(0)
data = Measures(model() .+ randn(rng, length(dom)), 1.)

# Fit the model against the data
bestfit = fit!(model, data)
```

The `f` symbol is the arbitrary name we chose to identify the component in the model, while `GFit.FuncWrap` is the component type which (as the name suggests) is simply a wrapper to a user defined function, accepting the function itself and an initial guess estimate of the parameters.  Note that no `Prediction` object is involved since we only have one domain, and no `Reducer` is required since we have only one component.

The best fit value and uncertainty for the `p1` parameter can be printed with:

```julia
println(bestfit[:f].p[1].val, " ± ", bestfit[:f].p[1].unc)
```

A plot of the mock data set and the best fit model can be generated as follows (here I will use [Gnuplot.jl](https://github.com/gcalderone/Gnuplot.jl/), but any other plotting framework would do the job):

```julia
using Gnuplot
@gp    domain(model)[:] data.val data.unc "w yerr t 'Data'"
@gp :- domain(model)[:] model() "w l t 'Model' lw 5"
```

### Multiple components

A model can contain many components, combined with a standard mathematical expression.  The previous example can easily be re-implemented by splitting the analytical formula in 3 parts, each assigned to a different `FuncWrap` component, and combined using the `@reducer` macro:

```julia
using GFit, Random

# Model formula and true parameter values
f1(x, p1, p2, p3) = @.  p1  +  p2 * x  +  p3 * x^2
f2(x, p4, p5) = @. p4 * sin(p5 * x)
f3(x) = cos.(x)
true_vals = [1, 2, 1.e-3, 4, 5]

# Prepare domain
dom = Domain(1:1e-2:10)

# Prepare a model with three components of type GFit.FuncWrap
model = Model(dom,
              @reducer((f1, f2, f3) -> (f1 .+ f2) .* f3),
              :f1 => GFit.FuncWrap(f1, true_vals[1], true_vals[2], true_vals[3]),
              :f2 => GFit.FuncWrap(f2, true_vals[4], true_vals[5]),
              :f3 => GFit.FuncWrap(f3))

# Create mock data set by evaluating the model and adding some
# random noise.  All uncertainties are equal to 1.
rng = MersenneTwister(0)
data = Measures(model() .+ randn(rng, length(dom)), 1.)

# Fit the model against the data
bestfit = fit!(model, data)

# Plot data and best fit model
using Gnuplot
@gp    domain(model)[:] data.val data.unc "w yerr t 'Data'"
@gp :- domain(model)[:] model() "w l t 'Model' lw 5"
```

Note that the anonymous function used in the `@reducer` macro is just a standard Julia function, whose inputs are standard Julia `Vector`s.

Also note that the results are identical as those obtained with the previous example, but the execution time is now supposed to be much shorter (verify using a much finer domain, e.g. `Domain(1.:1e-5:10)`).  Such improvement is due to the caching mechanism which operates at the component level: a component is evaluated only if necessary, i.e. when one of the parameter value is modified by the minimizer. In this particular case the `f3` component, having no free parameter, is evaluated only once.

To check how many time each component and reducer are evaluated simply type the name of the `Model` object in the REPL (or invoke `show(model)`, and check the `Eval. count` column.

#### Add new components to a model

The model can be built iteratively, by adding further components as they become necessary.  In the following example we will start with by neglecting `f3`, and add it subsequently:

```julia
# Fit a model with just 2 component (neglect f3)
model = Model(dom,
              :f1 => GFit.FuncWrap(f1, true_vals[1], true_vals[2], true_vals[3]),
              :f2 => GFit.FuncWrap(f2, true_vals[4], true_vals[5]))
bestfit = fit!(model, data)

# The fit is clearly a bad one. Add the third component and
# adjust initial guess values
add!(model, @reducer((f1, f2, f3) -> (f1 .+ f2) .* f3),
           :f3 => GFit.FuncWrap(f3))
model[:f1].p[1].val = 1
model[:f1].p[2].val = 2
model[:f1].p[3].val = 0
model[:f2].p[1].val = 4
model[:f2].p[2].val = 5

bestfit = fit!(model, data)
@gp    domain(model)[:] data.val data.unc "w yerr t 'Data'"
@gp :- domain(model)[:] model() "w l t 'Model' lw 5"
```

### Multiple data sets

In the examples above we never used the concept of *prediction* since we always dealt with a single dataset.  Now suppose you observe the same phenomenon with two different instruments, or at two different times, and wish to use both data sets to constrain the model parameters.  The following example shows how to fit multiple data sets simultaneously:

```julia
using GFit, Random

# Model formula and true parameter values
f1(x, p1, p2, p3) = @.  p1  +  p2 * x  +  p3 * x^2
f2(x, p4, p5) = @. p4 * sin(p5 * x)
f3(x) = cos.(x)
true_vals = [1, 2, 1.e-3, 4, 5]

# Prepare first domain
dom1 = Domain(1:1e-1:10)
dom2 = Domain(2:5e-2:9)

# Prepare list of common components in a Dict:
comps = Dict(:f1 => GFit.FuncWrap(f1, true_vals[1], true_vals[2], true_vals[3]),
             :f2 => GFit.FuncWrap(f2, true_vals[4], true_vals[5]),
             :f3 => GFit.FuncWrap(f3))

# Prepare two predictions, one for each domain. Here we also assume
# that the second data set requires a multiplicative calibration
# factor as free parameter
pred1 = Prediction(dom1,  @reducer((f1, f2, f3) -> (f1 .+ f2) .* f3), comps);
pred2 = Prediction(dom2,
                   @reducer((calib, f1, f2, f3) -> (f1 .+ f2) .* f3 .* calib),
                   :calib => 1.5, comps...);

# Prepare global model
model = Model([pred1 , pred2])

# Prepare mock data sets, including an inter-calibration factor of 1.3.
# Note: we now need to append an index to `model` in order to refer to a
# specific prediction.
rng = MersenneTwister(0)
data1 = Measures(model[1]() .+ randn(rng, length(dom1)), 1.)
data2 = Measures(model[2]() .+ randn(rng, length(dom2)), 1.)

# Fit the model against the data
bestfit = fit!(model, [data1, data2])

# Plot data and best fit model
using Gnuplot
@gp    domain(model[1])[:] data1.val data1.unc "w yerr t 'Data 1'"
@gp :- domain(model[1])[:] model[1]() "w l t 'Model 1' lw 5"
@gp :- domain(model[2])[:] data2.val data2.unc "w yerr t 'Data 2'"
@gp :- domain(model[2])[:] model[2]() "w l t 'Model 2' lw 5"
```

```julia
model[2][:f1].p[1].fixed = 1
model[2][:f1].p[2].fixed = 1
model[2][:f1].p[3].fixed = 1
model[2][:f2].p[1].fixed = 1
model[2][:f2].p[2].fixed = 1
patch!(model) do m
    m[2][:f1].p[1] = m[1][:f1].p[1]
    m[2][:f1].p[2] = m[1][:f1].p[2]
    m[2][:f1].p[3] = m[1][:f1].p[3]
    m[2][:f2].p[1] = m[1][:f2].p[1]
    m[2][:f2].p[2] = m[1][:f2].p[2]
end
```

## Retrieve results

The best fit results are available as a `BestFitResult` structure, and returned by the `fit!` fuction.  From this structure the user can retrieve the parameter best fit values and uncertainties, the number of data samples, the number of degrees of freedom, the total chi-squared (`cost`) and the fitting elapsed time in seconds, e.g.:

```julia
println(bestfit[:comp1].p[1].val)
println(bestfit[:comp2].p[2].unc)
println(bestfit.ndata)
println(bestfit.dof)
println(bestfit.cost)
println(bestfit.elapsed)
```

In the following example I will use the [Gnuplot.jl](https://github.com/gcalderone/Gnuplot.jl) package to show how to plot the results:

```julia
using Gnuplot
@gp "set key left" :-
@gp :- domain(model) data.val "w p tit 'Data'" :-
@gp :- domain(model) model(:comp1) "w l tit 'comp1'" :-
@gp :- domain(model) model(:comp2) "w l tit 'comp2'" :-
@gp :- domain(model) model() "w lines tit 'Model' lw 3"
```

Similarly, for the second dataset:

```julia
@gp "set key left" :-
@gp :- domain(model, id=2) data2.val "w p tit 'Data'" :-
@gp :- domain(model, id=2) model(:comp1, id=2) "w l tit 'comp1'" :-
@gp :- domain(model, id=2) model(:comp2, id=2) "w l tit 'comp2'" :-
@gp :- domain(model, id=2) model(id=2) "w lines tit 'Model' lw 3"
```

To evaluate the model with a different parameter value:

```julia
model[:comp2].p[1].val = 0
evaluate(model)
@gp :- domain(model, id=2) model(id=2) "w lines tit 'comp2, p[1]=0' lw 3"
```

Alternatively, a plot can be obtained using the companion package [GFitViewer.jl](https://github.com/lnicastro/GFitViewer.jl):

```julia
using GFitViewer
viewer(model, [data, data2], bestfit)
```

## Built-in components

### `Funcwrap`

The `FuncWrap` is simply a wrapper to a user defined function of the form `f(x, [y], [z], [further dimensions...], p1, p2, [further parameters...])`.  The `x`, `y`, `z` arguments will be `Vector{Float64}` with the same number of elements, while `p1`, `p2`, etc. will be scalar floats.  The function must return a `Vector{Float64}` (regardless of thenumber of dimensions) with the same number of elements as `x`.  This components works with domains of any dimensionality.

The constructor is defined as follows:

```julia
FuncWrap(func::Function, args...)
```

where `args...` is a list of numbers.

The parameters are:

- `p::Vector{Parameter}`: vector of parameters for the user defined function.

### `SimplePar`

The `SimplePar` represent a scalar component in the model, whose value is given by the `val` parameter.  This components works with domains of any dimensionality.

The constructor is defined as follows:

```julia
SimplePar(val::Number)
```

The parameters are:

- `val::Parameter`: the scalar value.

## Using the `CMPFit` minimizer

The `GFit` package uses the [LsqFit](https://github.com/JuliaNLSolvers/LsqFit.jl) minimizer by default, but it allows to use the [CMPFit](https://github.com/gcalderone/CMPFit.jl) as well.  The latter typically provides better performances with respect to the former, but since `CMPFit` is a wrapper to a C library it is not a pure-Julia solution.   The better performances are due to a different minimization strategy, not to C vs. Julia differences.

To use the `CMPFit` minimzer (after having installed the package):

```julia
using CMPFit
GFit.@with_CMPFit
result = fit!(model, [data, data2], minimizer=cmpfit())
```

## Multidimensional domains

**IMPORTANT NOTE**: by default the `GFit` package defines structure only for 1D and 2D fitting.  To handle higher dimensional cases you should trigger definition of proper code with:

```julia
GFit.define_ndim(3)
```

`N`-dimensional `Domain` objects comes in two flavours: linear and cartesian ones:

- Linear domain: contains a `N`-dimensional tuple of coordinates, one for each data sample.  It is similar to `Vector{NTuple{N, Number}}`;
- Cartesian domain: contains `N` arrays of coordinates, one for each dimension. Optionally contains a 1D list of index to select a subset of all possible combinations of coordinates. It is similar to `Vector{Vector{Number}}`, whose length of outer vector is `N`;

Linear domains are created using the `Domain` function, providing as many arguments as the number of dimensions. Each argument can either be an integer (specifying how many samples are defined along each axis), or a vector of `float`s.  **The length of all dimensions must be exactly the same.**  Examples:

```julia
# 1D
dom = Domain(5)
dom = Domain(1.:5)
dom = Domain([1,2,3,4,5.])

# 2D
dom = Domain(5, 5)
dom = Domain(1.:5, [1,2,3,4,5.])
```

Note that the length of all the above domains is 5.

Cartesian domains are created using the `CartesianDomain` function, providing as many arguments as the number of dimensions.  There is no 1-dimensional cartesian domain, hence `CartesianDomain` requires at least two arguments.  Each argument can either be an integer (specifying how many samples are defined along each axis), or a vector of `float`s.  **The length of dimensions may be different.**  Examples:

```julia
# 2D
dom = CartesianDomain(5, 6)
dom = CartesianDomain(1.:5, [1,2,3,4,5,6.])
```

The length of all the above domains is 30, i.e. it is equal to the product of the lengths of all dimensions.

Typically, the model can be evaluated over the cartesian product of all dimensions.  However, there can be specific locations of the domain for which there is not empirical data to compare with, making the model evaluation useless.  In these cases it is possible to select a subset of the cartesian domain using a 1D linear index, e.g.:

```julia
dom = CartesianDomain(1.:5, [1,2,3,4,5,6.], index=collect(0:4) .* 6 .+ 1)
```

The length of such domain is 5, equal to the length of the vector passed as `index` keyword.

A cartesian domain can always be transformed into a linear domain, while the inverse is usually not possible.  To check how a "*flattened*" version of a cartesian domain looks like you can use the `GFit.flatten` function, i.e.:

```julia
dom = CartesianDomain(1.:5, [1,2,3,4,5,6.], index=collect(0:4) .* 6 .+ 1)
lin = GFit.flatten(dom)
```

Actually, all models are always evaluated on "*flattened*", i.e. linear, domains.

To get the vector of coordinates for dimensions 1, 2, etc. of the `dom` object use the `dom[1]`, `dom[2]`, etc. syntax.  For linear domains all such vectors have the same length, for cartesian domains the lengths may differ.

## Multidimensional fitting

As an example we will fit a 2D plane.  The analytic function will accept two vectors for coordinates x and y, and 2 parameters.

```julia
f(x, y, p1, p2) = @. p1 * x  +  p2 * y
```

This function will be called during model evaluation, where only linear domains are involved, hence we are sure that both the `x` and `y` vectors will have the same lengths.

To create a multidimensional `Measures` object we will populate a 2D array and use it as argument to the `Measures` function:

```julia
dom = CartesianDomain(30, 40)
d = fill(0., size(dom));
for i in 1:length(dom[1])
    for j in 1:length(dom[2])
        d[i,j] = f(dom[1][i], dom[2][j], 1.2, 2.4)
    end
end
data = Measures(d + randn(rng, size(d)), 1.)
```

To fit the model proceed in the usual way:

```julia
model = Model(dom, :comp1 => GFit.FuncWrap(f, 1, 2))
result = fit!(model, data)
```

## Interactive Use

`GFit.jl` provides several facilities for interactive use on the REPL:

- all the types (i.e. `Domain`, `Measures`, `Model` and `BestFitResult`) have a dedicated `show` method to quickly and easily inspect their contents.  Simply type the name of the object to run this method;
- To get the list of currently defined components in a model simply type `model.comps[:<TAB>`;
- To get a list of parameter for a component simply type `model[:<COMPONENT NAME>]`, e.g. `model[:comp1]`.  Remember that the component parameter can either be scalar `Parameter` or a `Vector{Parameter}`;
- To get the list of model parameter in a result simply type `result.comps[:<TAB>`;

## Parameter settings

Each model parameter has a few settings that can be tweaked by the user before running the actual fit:

- `.val` (float): guess initial value;
- `.low` (float): lower limit for the value (default: `-Inf`);
- `.high` (float): upper limit for the value (default: `+Inf`);
- `.fixed` (bool): false for free parameters, true for fixed ones (default: `false`);
- `.expr` (string): a mathematical expression to bind the parameter value to other parameters (default: `""`);

**Important note**: the default minimizer ([LsqFit](https://github.com/JuliaNLSolvers/LsqFit.jl)) do not supports bounded parameters, while  [CMPFit](https://github.com/gcalderone/CMPFit.jl) supports them.

Considering the previous example we can limit the interval for `p1`, and fix the value for `p2` as follows:

```julia
model[:comp1].p[1].val  = 1   # guess initial value
model[:comp1].p[1].low  = 0.5 # lower limit
model[:comp1].p[1].high = 1.5 # upper limit
model[:comp1].p[2].val  = 2.4
model[:comp1].p[2].fixed = true
result = fit!(model, data, minimizer=cmpfit())
```

To remove the limits on `p1` simply set their bounds to +/- Inf:

```julia
model[:comp1].p[1].low  = -Inf
model[:comp1].p[1].high = +Inf
```

## Parameter patching

A parameter may be "*patched*", i.e. its value being re-calculated before actually evaluating the model, through a common  Julia expression and set using the `@patch` macro.  A common case is to patch a parameter to link its value to another one:

```julia
@patch!(model, model[:comp1].p[2] .= model[:comp1].p[1])
model[:comp1].p[2].fixed = true
```

Note that in the above example we had to fix the `p[2]` parameter otherwise the minizer will try to find a best fit for a parameter which has no influence on the final model, since its value will always be overwritten by the expression.

Another possibility is to use one parameter value to calculate another one's:

```julia
@patch!(model, model[:comp1].p[2] .+= model[:comp1].p[1])
```

## Component settings

TODO: `thaw`, `freeze`

# Built-in components

The following components are available in the `GFit.Components` module:

- OffsetSlope (1D and 2D): an offset and slope component;
- Polynomial (only 1D): a n-th degree polynomial function (n > 1);
- Gaussian (1D and 2D): a Gaussian function;
- Lorentzian (1D and 2D): a Lorentzian function;

## OffsetSlope

An offset and slope component for 1D and 2D domains.  In 2D it represents a tilted plane.

The constructors are defined as follows:

- 1D: `GFit.Components.OffsetSlope(offset, x0, slope)`;
- 2D: `GFit.Components.OffsetSlope(offset, x0, y0, slopeX, slopeY)`;

The parameters are:

- 1D:
  - `offset::Parameter`: a global offset;
  - `x0::Parameter`: the X coordinate of the point where the component equals `offset`.  This parameter is fixed by default;
  - `slope::Parameter`: the slope of the linear function;
- 2D:
  - `offset::Parameter`: a global offset;
  - `x0::Parameter`: the X coordinate of the point where the component equals `offset`.  This parameter is fixed by default;
  - `y0::Parameter`: the Y coordinate of the point where the component equals `offset`.  This parameter is fixed by default;
  - `slopeX::Parameter` (only 2D): the slope of the plane along the X direction;
  - `slopeY::Parameter` (only 2D): the slope of the plane along the Y direction;

## Polynomial

A n-th degree polynomial function (n > 1) for 1D domains.

The constructor is defined as follows:

- `GFit.Components.Polynomial(args...)`;
  where `args...` is a list of numbers.

The parameters are:

- `coeff::Vector{Parameter}`: vector of polynomial coefficients.

## Gaussian

A normalized Gaussian component for 1D and 2D domains.

The constructors are defined as follows:

- 1D: `GFit.Components.Gaussian(norm, center, sigma)`;
- 2D: `GFit.Components.Gaussian(norm, centerX, centerY, sigma)` (implies `sigmaX=sigmaY`, `angle=0`);
- 2D: `GFit.Components.Gaussian(norm, centerX, centerY, sigmaX, sigmaY, angle)`;

The parameters are:

- 1D:
  
  - `norm::Parameter`: the area below the Gaussian function;
  - `center::Parameter`: the location of the center of the Gaussian;
  - `sigma::Parameter`: the width the Gaussian;

- 2D:
  
  - `norm::Parameter`: the volume below the Gaussian function;
  - `centerX::Parameter`: the X coordinate of the center of the Gaussian;
  - `centerY::Parameter`: the Y coordinate of the center of the Gaussian;
  - `sigmaX::Parameter`: the width the Gaussian along the X direction (when `angle=0`);
  - `sigmaY::Parameter`: the width the Gaussian along the Y direction (when `angle=0`);
  - `angle::Parameter`: the rotation angle of the whole Gaussian function.

## Lorentzian

A Lorentzian component for 1D and 2D domains.

The constructors are defined as follows:

- 1D: `GFit.Components.Lorentzian(norm, center, fwhm)`;
- 2D: `GFit.Components.Lorentzian(norm, centerX, centerY, fwhmX, fwhmY)`;

The parameters are:

- 1D:
  
  - `norm::Parameter`: the area below the Lorentzian function;
  - `center::Parameter`: the location of the center of the Lorentzian;
  - `fwhm::Parameter`: the full-width at half maximum of the Lorentzian;

- 2D:
  
  - `norm::Parameter`: the volume below the Lorentzian function;
  - `centerX::Parameter`: the X coordinate of the center of the Lorentzian;
  - `centerY::Parameter`: the Y coordinate of the center of the Lorentzian;
  - `fwhmX::Parameter`: the full-width at half maximum of the Lorentzian along the X direction (when `angle=0`);
  - `fwhmY::Parameter`: the full-width at half maximum of the Lorentzian along the Y direction (when `angle=0`).

## Examples:

### 1D: offset + two Gaussian profiles

```julia
x = Domain(1:0.05:10)
model = Model(x,
    :offset => 4,
    :line1  => GFit.Gaussian(1.1 , 4.4, 0.51),
    :line2  => GFit.Gaussian(0.52, 5.5, 1.2 ))

using Random
rng = MersenneTwister(0);
noise = maximum(model()) * 0.01
data = Measures(model() + noise * randn(rng, length(model())), noise);
ret1 = fit!(model, data)
```

To produce the plots I will use the [Gnuplot.jl](https://github.com/gcalderone/Gnuplot.jl) package, but the user can choose any other package:

```julia
using Gnuplot
@gp    "set multi layout 2,1" :-
@gp :- domain(model) data.val data.unc "w yerr tit 'Data'" :-
@gp :- domain(model) model(:line1) .+ model(:offset) "w l tit 'offset + line1'" :-
@gp :- domain(model) model(:line2) .+ model(:offset) "w l tit 'offset + line2'" :-
@gp :- domain(model) model() "w lines tit 'Model' lw 3" :-
@gp :- 2 x[1] (data.val - model()) ./ data.unc fill(1., length(data)) "w yerr tit 'Residuals'"
```

### 2D: tilted plane + 2D Gaussian profile

```julia
dom = CartesianDomain(-5:0.1:5, -4:0.1:4)
model = Model(dom,
              :background => GFit.OffsetSlope(0, 0, 0., 2., 3.),
              :psf => GFit.Gaussian(100., 0., 0., 1, 0.3, 15))


noise = maximum(model()) * 0.1
data = Measures(model() .+ 4 .+ noise .* randn(length(model())), noise);
ret1 = fit!(model, data)
```

To produce the plots I will use the [Gnuplot.jl](https://github.com/gcalderone/Gnuplot.jl) package, but the user can choose any other package:

```julia
using Gnuplot

# Plot the model...
@gsp dom[1] dom[2] reshape(model(), dom)

# ...and the residuals
@gsp dom[1] dom[2] reshape(data.val - model(), dom)

# Plot using pm3d style
@gsp "set pm3d" "set palette" dom[1] dom[2] reshape(model(), dom) "w dots"
```
