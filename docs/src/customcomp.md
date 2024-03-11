```@setup abc
include("setup.jl")
```

# Custom components

Besides the [Built-in components](@ref), the user may define any number of custom components to be used in a model.

A user-defined component shall satisfy the following constraints:
- it shall be a structure inheriting from `GModelFit.AbstractComponent`;

- the component values shall be defined as fields with type `Parameter`, e.g.:
```julia
struct MyComponent <: AbstractComponent
	param1::Parameter
	param2::Parameter
	...
end
```
(see below for a complete example).
Alternatively, the parameters may be specified as a single field of type `OrderedDict{Symbol, Parameter}` (see the [`Polynomial`](https://github.com/gcalderone/GModelFit.jl/blob/master/src/components/Polynomial.jl) component for an example);

- the `GModelFit.evaluate!` function shall be extended to provide the component-specific code for evaluation.
Specifically, the `evaluate!` function should replace the content of the `buffer::Vector{Float64}` with the outcome of the new component evaluation, given the numerical values for the parameters, e.g.
```julia
function evaluate!(ceval::CompEval{MyComponent},
                   param1::Float64, param2::Float64...)
	ceval.buffer .= (component evaluation using param1 and param2 values)
end
```


Optionally, the user may chose to extend also the following functions:
- `prepare!`: to prepare the component for evaluations on a specific domain (e.g. to pre-compute quantities which depend only on the domain being used). The return value must be the buffer (of type `Vector{Float64}`) to accomodate the component evaluation.  The default implementation simply creates a buffer with the same length as the input domain;

- `dependencies`: return a `Vector{Symbol}` specifying the names of the component dependencies.  The evaluation of the latter will be made available as arguments to the `evaluate!` method. The default implementation returns an empty list `Symbol[]`.


### Life cycle of a component

All components "live" within `Model` object, which has a well defined domain associated (see constructor of `Model`).  The same domain will also be used for component evaluations.  The typical life cycle of a component is as follows:
1. the component is created invoking its constructor and providing an initial guess values for all parameters;
1. the component is added to the model. In this step the `prepare!` function is called to precompute component quantities (according to the domain associated to the `Model` object) and to allocate the buffer for evaluations.  Note that the `prepare!` function is invoked only once for each component;
1. the user may optionally modify the component parameter guess values, as well as their [Parameter constraints](@ref), before starting the fit;
1. before starting the fit process a dependency tree is generated to identify all component dependencies (by invoking the `dependencies()` function for all components). The tree is used to identify the proper order for component evaluations: the leaves will be evaluated first, and the remaining ones will be evaluated following the tree branches.  The evaluation of the last component (the root of the tree) will be compared to the data;
1. during the fitting process the component `evaluate!()` function is invoked whenever a change in its parameter values is detected, until a convergence criterion is satisfied.


### Example

A common case is to compare empirical data with a numerically evaluated theoretical model, possibly defined on a different grid with respect to the empirical one.  An interpolation is therefore required in order to compare the model to the data.

Let's assume the theoretical model is defined as follows:
```@example abc
theory_x = 0.:10
theory_y = [0, 0.841, 0.909, 0.141, -0.757, -0.959, -0.279, 0.657, 0.989, 0.412, -0.544]
println() # hide
```
while the empirical data are:
```@example abc
obs_x = [0.500, 2.071, 3.642, 5.212, 6.783, 8.354, 9.925]
obs_y = [2.048, 3.481, 1.060, 0.515, 3.220, 4.398, 1.808]
println() # hide
```

The following example shows how to implement a component aimed to interpolate the theoretical model onto a specific empirical domain, with the only parameter being a global scaling factor:
```@example abc
using GModelFit, Interpolations
import GModelFit.prepare!, GModelFit.evaluate!

# Define the component structure and constructor
struct Interpolator <: GModelFit.AbstractComponent
	theory_x::Vector{Float64}
	theory_y::Vector{Float64}
	interp_y::Vector{Float64}  # will contain the interpolated values
	scale::GModelFit.Parameter

	function Interpolator(theory_x, theory_y)
		scale = GModelFit.Parameter(1)
		scale.low = 0                  # ensure scale parameter is positive
		interp_y = Vector{Float64}()   # this will be populated in prepare!()
		return new(theory_x, theory_y, interp_y, scale)
	end
end

# Component preparation: invoked only once to precompute quantities
# and allocate evaluation buffer
function prepare!(comp::Interpolator, domain::AbstractDomain{1})
	# Pre-compute interpolation on the empirical domain
	itp = linear_interpolation(comp.theory_x, comp.theory_y)
	append!(comp.interp_y, itp(coords(domain)))
	return fill(NaN, length(comp.interp_y)) # buffer for evaluations
end

# Component evaluation
function evaluate!(ceval::GModelFit.CompEval{Interpolator, <: AbstractDomain{1}},
                   scale)
	ceval.buffer .= scale .* ceval.comp.interp_y
end
println() # hide
```

The following code shows how to prepare a `Model` including the interpolated theoretical model, and to take into account the possible background introduced by the detector used to obtain empirical data:
```@example abc
dom = Domain(obs_x)
model = Model(:theory => Interpolator(theory_x, theory_y),
              :background => GModelFit.OffsetSlope(1., 0., 0.2),
              :main => SumReducer(:theory, :background))
data = Measures(dom, obs_y, 0.2)
best, fitstats = fit(model, data)
dumpjson("ex_Customcomp", best, fitstats, data) # hide
show((best, fitstats)) # hide
```
