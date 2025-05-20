```@setup abc
include("setup.jl")
```

# Custom components

Besides the [Built-in components](@ref), the user may define any number of custom components.  The latter are structures satisfying the following constraints:

- They must be structures inheriting from `GModelFit.AbstractComponent`;

- The component parameters (if present) must be defined as fields with type `Parameter`, e.g.:
  ```julia
  struct MyComponent <: AbstractComponent
      param1::Parameter
      param2::Parameter
      ...
  end
  ```
  Alternatively, the parameters may be specified as a single field of type `OrderedDict{Symbol, Parameter}` (see the [`Polynomial`](https://github.com/gcalderone/GModelFit.jl/blob/master/src/components/Polynomial.jl) component for an example).  The structure may also contain further fields of any type;

- The [`GModelFit.evaluate!`](@ref) function should be extended with a dedicated method to evaluate the component, as shown below.


Optionally, the user may choose to extend also the following functions:
- [`GModelFit.prepare!`](@ref): to pre-compute quantities depending on the evaluation domain, and to allocate the evaluation buffer;

- [`GModelFit.dependencies`](@ref): to specify the list of the component dependencies.


Note that before being evaluated, all components need to be wrapped into a `CompEval` structure, and that the above mentioned `evaluate!` function requires a `CompEval` object as first argument.  The outcomes of the evaluations should be placed in `CompEval.buffer`.

The following example shows how to define a custom component:

```@example abc
using GModelFit
import GModelFit: prepare!, evaluate!

struct MyComponent <: GModelFit.AbstractComponent
    param1::GModelFit.Parameter

    function MyComponent(param1)
        println(" -> call to MyComponent constructor;")
        new(GModelFit.Parameter(param1))
     end
end

function prepare!(comp::MyComponent, domain::AbstractDomain)
    println(" -> call to prepare!()")
    return fill(NaN, length(domain)) # buffer for evaluations
end

function evaluate!(ceval::GModelFit.CompEval{MyComponent, <: AbstractDomain{1}},
                   param1)
    println(" -> call to evaluate!() with parameter value: ", param1)
    ceval.tpar.buffer .= param1
end

println() # hide
```



## Life cycle of a component

  The life cycle of a component is as follows:

1. The component is created by invoking its constructor, and is added to a [`Model`](@ref) object;

1. When the [`fit!`](@ref) function is invoked, all components in a `Model` are wrapped into `CompEval` objects;

    - During creation of the `CompEval` structure the [`GModelFit.prepare!`](@ref) function is invoked to allocate the proper buffer for evaluations.  Note that the `prepare!` function is called only once for each `fit!` invocation, hence it is the perfect place to pre-compute quantities which will be used during the component evaluation;

1. During the minimization process the [`GModelFit.evaluate!`](@ref) function is repeatedly invoked to evalute the component varying the parameter values until a convergence criterion is met.

The following example shows how to simulate the life cycle for the `MyComponent` structure defined above:
```@example abc
# Create a component and a domain for evaluation
comp = MyComponent(1)
dom = Domain(1:5)

# Create CompEval object (the `prepare!` function is invoked here):
ceval = GModelFit.CompEval(comp, dom)

# Repeated evaluations varying parameter value:
GModelFit.evaluate!(ceval, 1)
GModelFit.evaluate!(ceval, 2)
GModelFit.evaluate!(ceval, 3)

# Retrieve results
println(ceval.tpar.buffer)
```

The actual life cycle during minimization is slightly more complex since the `evaluate!` function is invoked only if a change in the parameter values with respect to previous evaluation has been detected.



## Complete example

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

The following example shows how to implement a component which interpolates a theoretical model onto a specific empirical domain, with the only parameter being a global scaling factor:
```@example abc
using GModelFit, Interpolations
import GModelFit: prepare!, result_length, evaluate!

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
end

result_length(comp::Interpolator, domain::AbstractDomain{1}) = length(domain)

# Component evaluation (apply scaling factor)
function evaluate!(comp::Interpolator, ::AbstractDomain{1}, output,
                   scale)
    output .= scale .* comp.interp_y
end
println() # hide
```

The following code shows how to prepare a `Model` including the interpolated theoretical model, and to take into account the possible background introduced by the detector used to obtain empirical data:
```@example abc
model = Model(:theory => Interpolator(theory_x, theory_y),
              :background => GModelFit.OffsetSlope(1., 0., 0.2),
              :main => SumReducer(:theory, :background))

data = Measures(Domain(obs_x), obs_y, 0.2)
bestfit, stats = fit(model, data)
show((bestfit, stats)) # hide
```
