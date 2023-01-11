# Built-in components

The **GFit.jl** provides several built-in components which may be used to build arbitrarily complex models.


## λComp
The `λComp` component is a wrapper for a user defined lambda function (namely a [`GFit.λFunct`](@ref) object, as obtained by the [`@λ`](@ref) macro).  The function is typically a mathematical expression combining any number of parameters and/or other component evaluations within the same model.  The expression should be given in the form:
```
@λ (x, [y, [further domain dimensions...],]
    [comp1, [comp2, [further components ...],]]
    [p1=guess value, [p2=guess value, [further parameters]]]) ->
	(mathematical expression)
```
where the mathematical expression returns a `Vector{Float64}` with the same length as the model domain (regardless of the number of dimensions involved).


The `λComp` constructor is defined as follows:

```julia
λComp(f::λFunct)
```
however there is no need to explicitly invoke such constructor since a `λComp` object is automatically created whenever a `λFunct` is added to the model. This components works with domains of any dimensionality.


### Example
```@example abc
using GFit

# Prepare domain and a linear model (with initial guess parameters)
dom = Domain(1:5)
model = Model(dom, :linear => @λ (x, b=2, m=0.5) -> (b .+ x .* m))

# Fit model against data
data = Measures(dom, [4.01, 7.58, 12.13, 19.78, 29.04], 0.4)
res = fit!(model, data)
```

The evaluation of a `λComp` component may also involve the outcomes from other components. Continuing from previous example, whose fit was clearly a poor one, we may add a quadratic term to the previously defined `linear` component:
```@example abc
model[:quadratic] = @λ (x, linear, p2=1) -> (linear .+ p2 .* x.^2)
res = fit!(model, data)
```

Note that the keywords given when defining the function are interpreted as component parameters, hence their properties can be retrieved with:
```@example abc
println("Best fit values:")
println("b:  ", res.bestfit[:linear].b.val , " ± ", res.bestfit[:linear].b.unc)
println("m:  ", res.bestfit[:linear].m.val , " ± ", res.bestfit[:linear].m.unc)
println("p2: ", res.bestfit[:quadratic].p2.val, " ± ", res.bestfit[:quadratic].p2.unc)
```




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


### Example
```@example abc
using GFit

# Prepare domain and a linear model using the OffsetSlope component
dom = Domain(1:5)
model = Model(dom, :linear => GFit.OffsetSlope(2, 0, 0.5))

# Fit model against data
data = Measures(dom, [4.01, 7.58, 12.13, 19.78, 29.04], 0.4)
res = fit!(model, data)
```

Note that the numerical results are identical to the previous example (where an explicit mathematical expression for the `:linear` component was used).  Also, in this case, you may simply add a new quadratic term:
```@example abc
model[:quadratic] = @λ (x, linear, p2=1) -> (linear .+ p2 .* x.^2)
res = fit!(model, data)
```




```@example abc
using GFit

# Prepare domain and a linear model using the OffsetSlope component
dom = CartesianDomain(1:5, 1:5)
model = Model(dom, :plane => GFit.OffsetSlope(2, 0, 0, 0.5, 0.5))

# Fit model against data
data = Measures(dom, [ 3.08403  3.46719  4.07612  4.25611  5.04716
                       3.18361  3.88546  4.52338  5.12838  5.7864
                       3.80219  4.90894  5.24232  5.06982  6.29545
                       4.34554  4.68698  5.51505  5.69245  6.35409
                       4.643    5.91825  6.18011  6.67073  7.01467], 0.25)
res = fit!(model, data)
```






### Polynomial

A n-th degree polynomial function (n > 1) for 1D domains.

The constructor is defined as follows:

- `GFit.Components.Polynomial(args...)`;
  where `args...` is a list of numbers.

The parameters are:

- `coeff::Vector{Parameter}`: vector of polynomial coefficients.

### Gaussian

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

### Lorentzian

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
