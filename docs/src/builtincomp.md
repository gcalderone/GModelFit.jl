```@setup abc
include("setup.jl")
```


# Built-in components

The **GModelFit.jl** provides several built-in components which may be used to build arbitrarily complex models.

## OffsetSlope

An offset and slope component for 1D and 2D domains.

The constructors are defined as follows:

- 1D: `GModelFit.OffsetSlope(offset, x0, slope)`;
- 2D: `GModelFit.OffsetSlope(offset, x0, y0, slopeX, slopeY)`;

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


#### Example
```@example abc
using GModelFit

# Define a linear model using the OffsetSlope component
model = Model(:linear => GModelFit.OffsetSlope(2, 0, 0.5))

# Fit model against data
data = Measures([4.01, 7.58, 12.13, 19.78, 29.04], 0.4)
bestfit, fsumm = fit(model, data)
show((bestfit, fsumm)) # hide
```

The best fit parameter values can be retrieved with:
```@example abc
println("Best fit values:")
println("b:  ", bestfit[:linear].offset.val, " ± ", bestfit[:linear].offset.unc)
println("m:  ", bestfit[:linear].slope.val , " ± ", bestfit[:linear].slope.unc)
```

A similar example in 2D is as follows:
```@example abc
using GModelFit

# Define a linear model using the OffsetSlope component
model = Model(:plane => GModelFit.OffsetSlope(2, 0, 0, 0.5, 0.5))

# Fit model against data
dom = CartesianDomain(1:5, 1:5)
data = Measures(dom, [ 3.08403  3.46719  4.07612  4.25611  5.04716
                       3.18361  3.88546  4.52338  5.12838  5.7864
                       3.80219  4.90894  5.24232  5.06982  6.29545
                       4.34554  4.68698  5.51505  5.69245  6.35409
                       4.643    5.91825  6.18011  6.67073  7.01467], 0.25)
bestfit, fsumm = fit(model, data)
show((bestfit, fsumm)) # hide
```



## Polynomial

A *n*-th degree polynomial function (*n > 1*) for 1D domains.

The constructor is defined as follows:

- `GModelFit.Polynomial(p1, p2, ...)`;
  where `p1`, `p2`, etc. are the guess values for the coefficients of each degree of the polynomial.

The parameters are accessible as `p0`, `p1`, etc.

#### Example
```@example abc
using GModelFit

# Define domain and a linear model using the Polynomial component
model = Model(GModelFit.Polynomial(2, 0.5))

# Fit model against data
data = Measures([4.01, 7.58, 12.13, 19.78, 29.04], 0.4)
bestfit, fsumm = fit(model, data)
show((bestfit, fsumm)) # hide
```

Note that the numerical results are identical to the previous example involving the `OffsetSlope` component.  Also note that the default name for a component (if none is provided) is `:main`.  To use a 2nd degree polynomial we can simply replace the `:main` component with a new one:
```@example abc
model[:main] = GModelFit.Polynomial(2, 0.5, 1)
bestfit, fsumm = fit(model, data)
show((bestfit, fsumm)) # hide
```


## Gaussian

A normalized Gaussian component for 1D and 2D domains.

The constructors are defined as follows:

- 1D: `GModelFit.Gaussian(norm, center, sigma)`;
- 2D: `GModelFit.Gaussian(norm, centerX, centerY, sigma)` (implies `sigmaX=sigmaY`, `angle=0`);
- 2D: `GModelFit.Gaussian(norm, centerX, centerY, sigmaX, sigmaY, angle)`;

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
  - `angle::Parameter`: the rotation angle (in degrees) of the Gaussian.



#### Example
```@example abc
using GModelFit

# Define a model with a single Gaussian component
model = Model(GModelFit.Gaussian(1, 3, 0.5))

# Fit model against data
data = Measures([0, 0.3, 6.2, 25.4, 37.6, 23., 7.1, 0.4, 0], 0.6)
bestfit, fsumm = fit(model, data)
show((bestfit, fsumm)) # hide
```


A very common problem is to fit the histogram of a distribution with a Gaussian model.  The following example shows how to fit such Gaussian model to a distribution generated with `Random.randn`, and how to plot the results using [Gnuplot.jl](https://github.com/gcalderone/Gnuplot.jl/):
```@example abc
using Random, GModelFit, Gnuplot

# Calculate histogram of the distribution
hh = hist(randn(10000), bs=0.25)

# Define domain and data and fit a model
dom = Domain(hist_bins(hh, side=:center, pad=false))
data = Measures(dom, hist_weights(hh, pad=false), 1.)
model = Model(GModelFit.Gaussian(1e3, 0, 1))
bestfit, fsumm = fit(model, data)
show((bestfit, fsumm)) # hide
```

```@example abc
@gp hh coords(dom) bestfit() "w l t 'Model' lw 3"
saveas("gaussian") # hide
```
![](assets/gaussian.png)



A similar problem in 2D can be handled as follows:

```@example abc
using Random, GModelFit, Gnuplot

# Calculate histogram of the distribution
hh = hist(1 .+ randn(10000), 2 .* randn(10000))

# Define domain and data and fit a model
dom = CartesianDomain(hist_bins(hh, 1), hist_bins(hh, 2))
data = Measures(dom, hist_weights(hh) .* 1., 1.)
model = Model(GModelFit.Gaussian(1e3, 0, 0, 1, 1, 0))
bestfit, fsumm = fit(model, data)
show((bestfit, fsumm)) # hide
```


## FComp

As anticipated in [Basic concepts and data types](@ref) any Julia function can be used as a component to evaluate.  The corresponding component type is `FComp`, whose constructors are defined as follows:

```julia
FComp(funct::Function, deps=Symbol[]; par1=guess1, par2=guess2, ...)
FComp(funct::FunctDesc)
```
In the first constructor `funct` is the Julia function, `deps` is a vector of dependencies (either the domain dimensions or other component names) and `par1`, `par2` etc. are the named parameters with their corresponding initial guess values.

#### Example
```@example abc
using GModelFit

# Define a simple Julia function to evaluate a linear relationship
myfunc(x, b, m) = b .+ x .* m

# Define a model with a `FComp` wrapping the previously defined function.
# Also specify the initial guess parameters.
model = Model(:linear => GModelFit.FComp(myfunc, [:x], b=2, m=0.5))

# Fit model against a data set
data = Measures([4.01, 7.58, 12.13, 19.78, 29.04], 0.4)
bestfit, fsumm = fit(model, data)
show((bestfit, fsumm)) # hide
```

In the second constructor a [`GModelFit.FunctDesc`](@ref) object is accepted, as generated by the [`@fd`](@ref) macro).  The function is typically a mathematical expression combining any number of parameters and/or other component evaluations within the same model.  The expression should be given in the form:
```
@fd (x, [y, [further domain dimensions...],]
    [comp1, [comp2, [further components ...],]]
    [par1=guess1, [par2=guess2, [further parameters]]]) ->
    (mathematical expression)
```
where the mathematical expression returns a `Vector{Float64}` with the same length as the model domain.

The previous example can be rewritten as follows:

```@example abc
using GModelFit

# Define a linear model (with initial guess parameters)
model = Model(:linear => @fd (x, b=2, m=0.5) -> (b .+ x .* m))

# Fit model against data
data = Measures([4.01, 7.58, 12.13, 19.78, 29.04], 0.4)
bestfit, fsumm = fit(model, data)
show((bestfit, fsumm)) # hide
```
Note that a `FComp` component can be added to a model without explicitly invoking its constructor when the [`@fd`](@ref) macro is used.



The evaluation of a `FComp` component may also involve the outcomes from other components. Continuing from previous example, whose fit was clearly a poor one, we may add a quadratic term to the previously defined `linear` component:
```@example abc
model[:quadratic] = @fd (x, linear, p2=1) -> (linear .+ p2 .* x.^2)
bestfit, fsumm = fit(model, data)
show((bestfit, fsumm)) # hide
```

The keywords given when defining the function are interpreted as component parameters, hence their properties can be retrieved with:
```@example abc
println("Best fit values:")
println("b:  ", bestfit[:linear].b.val    , " ± ", bestfit[:linear].b.unc)
println("m:  ", bestfit[:linear].m.val    , " ± ", bestfit[:linear].m.unc)
println("p2: ", bestfit[:quadratic].p2.val, " ± ", bestfit[:quadratic].p2.unc)
```

## SumReducer

A component calculating the element-wise sum of a number of other components.

The `SumReducer` constructor is defined as follows:
```julia
SumReducer(args::AbstractSet{Symbol})
SumReducer(args::Vector{Symbol})
SumReducer(args::Vararg{Symbol})
```
where the `Symbol`s represent the component names

The `SumReducer` component has no parameter.

#### Example
```@example abc
using GModelFit

# Define domain and a linear model (with initial guess parameters)
model = Model(:linear => @fd (x, b=2, m=0.5) -> (b .+ x .* m))

# Add a quadratic component to the model
model[:quadratic] = @fd (x, p2=1) -> (p2 .* x.^2)

# The total model is the sum of `linear` and `quadratic`
model[:main] = SumReducer(:linear, :quadratic)

# Fit model against data
dom = Domain(1:5)
data = Measures(dom, [4.01, 7.58, 12.13, 19.78, 29.04], 0.4)
bestfit, fsumm = fit(model, data)
show((bestfit, fsumm)) # hide
```
