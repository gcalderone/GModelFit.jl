# GFit.jl

`GFit` is a general purpose, data-driven model fitting framework for Julia.

[![Build Status](https://travis-ci.org/gcalderone/GFit.jl.svg?branch=master)](https://travis-ci.org/gcalderone/GFit.jl)

It provides the basic tools to build, inspect and fit complex models against empirical data.

The typical use case for `GFit` is as follows: you observe a physical phenomenon with one (or more) instrument(s) and wish to fit those data against a (possibly very complex) theoretical model, in order to extract the characterizing quantities represented by the model parameters.  `GFit` provides the following functionalities:
- it handles data of any dimensionality;
- the fitting model is evaluated on a user provided domain;
- the fitting model is built up by combining one or more *components*, either builtin or implemented by the user, using a standard Julia mathematical expression;
- all components results are cached so that repeated evaluations with the same parameter values do not involve further calculations;
- User provided components can pre-compute quantities based on the model domain, and store them in a private structure;
- Model parameters can be fixed to a specific value, limited in an interval, or be dynamically calculated (patched) using a mathematical expression involving the other parameters;
- multiple data sets can be fitted simultaneously;
- it allows to use different minimizers, and compare their results and performances (currently only two minimizers are supported: [LsqFit](https://github.com/JuliaNLSolvers/LsqFit.jl) and [CMPFit](https://github.com/gcalderone/CMPFit.jl));
- it provides several facilities for interactive fitting and results displaying.

See below for a simple example, and the `examples` directory for more complex ones.

## Installation
```julia
] add GFit
```

## Simple example
Assume the model to be compared with empirical data has 5 parameters and the following analytical formula:
```julia
f(x, p1, p2, p3, p4, p5) = @. (p1  +  p2 * x  +  p3 * x^2  +  p4 * sin(p5 * x))  *  cos(x)
```

To simulate a measurement process we'll evaluate the model on a domain and add some random noise to a model realization:

```julia
# Actual model parameters:
params = [1, 1.e-3, 1.e-6, 4, 5]

# Domain for model evaluation
x = 1.:50:10000

# Evaluated model
y = f(x, params...);

# Random noise
using Random
rng = MersenneTwister(0);
noise = randn(rng, length(x));
```

In order to use the `GFit` framework we must create a `Domain` and a `Measures` objects to encapsulate the model domain and empirical data:
```julia
using GFit
data = Measures(y .+ noise, 1.)
```
The second argument to the `Measures` function is the (1 sigma Gaussian) uncertainty associated to each data sample.  It can either be an array with the same shape as the first argument or a scalar.  In the latter case all data samples are assumed to have the same uncertainty.

Also, we must create a `Model` object containing a reference to the analytical formula, and prepare it for evaluation on the domain `dom`
```julia
model = Model(Domain(x), :comp1 => GFit.FuncWrap(f, params...))
```
The `comp1` symbol is the name we chose to identify the component in the model.  Any other valid symbol could have been used.

The parameter initial values are those given in the component constructors.  Such values can be changed as follows:
```julia
model[:comp1].p[1].val = 1
model[:comp1].p[2].val = 1.e-3
 ```
etc.


Finally, we are ready to fit the model against empirical data:
```julia
bestfit = fit!(model, data)
```
Note that the `fit!` function modifies the `model1` objects as follows:
- it updates the model parameter with the best fit ones;
- it updates the internal cache of component evaluations with those resulting from best fit parameters;
- it updates the evaluation counter for each component (see below);

The procedure outlined above may seem cumbersome at first, however it is key to define very complex models and to improve performances, as shown below.



## Multiple components

A model is typically made up of many *components*, joined toghether with a standard mathematical expression.  The previous example can easily be re-implemented by splitting the analytical formula in 3 parts:
```julia
f1(x, p1, p2, p3) = @.  p1  +  p2 * x  +  p3 * x^2
f2(x, p4, p5) = @. p4 * sin(p5 * x)
f3(x) = cos.(x)

model = Model(Domain(x),
              @reducer((comp1, comp2, comp3) -> (comp1 .+ comp2) .* comp3),
              :comp1 => GFit.FuncWrap(f1, params[1], params[2], params[3]),
              :comp2 => GFit.FuncWrap(f2, params[4], params[5]),
              :comp3 => GFit.FuncWrap(f3))
bestfit = fit!(model, data)
```
Now we used 3 components (named `comp1`, `comp2` and `comp3`) and combined them with the mathematical expression in the second argument (the so-called *reducer*). **Any** valid mathematical expression is allowed.

Note that we obtained exactly the same result as before, but we **gained a factor ~3** in execution time.  Such perfromance improvement is due to the fact that the component evaluations are internally cached, and are re-evaluated only if necessary, i.e. when one of the parameter value is modified by the minimizer algorithm. In this particular case the `comp3` component, having no free parameter, is evaluated only once.

To check how many time each component and the model are evaluated simply type the name of the `Model` object in the REPL or call the `show` method, i.e.: `show(model)`, and check the `Counter` column.


## Fitting multiple data sets
`GFit` allows to fit multiple data sets simultaneously.

Suppose you observe the same phenomenon with two different instruments and wish to use both data sets to constrain the model parameters.  Here we will simulate a second data sets by adding random noise to the previously calculated model, and creating a second `Measures` object:
```julia
noise = randn(rng, length(x));
data2 = Measures(1.3 * (y + noise), 1.3)
```
Note that we multiplied all data by a factor 1.3, to simulate a different calibration between the instruments.  To take into account such calibration we add a second *prediction* into the model, as well as a further scalar component named `calib`:
```julia
add!(model, Prediction(Domain(x),
            @reducer((calib, comp1, comp2, comp3) -> calib .* (comp1 .+ comp2) .* comp3),
            :calib => 1, model.comps...))
```

Note that the new prediction uses a different expression to be evaluated.  If needed, also the `Domain` object may be different from the first one.  The latter posssibility allows, for instance, to fit multiple data sets each observed with different instruments.

Now we can fit both data sets as follows:
```julia
bestfit = fit!(model, [data, data2])
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
@patch!(model, model[:comp1].p[2] .= model[:comp1].p[1]
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
