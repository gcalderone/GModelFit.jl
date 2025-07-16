```@setup abc
include("setup.jl")
```

# Custom components

Besides the [Built-in components](@ref), the user may define any number of custom components.  The latter are structures satisfying the following constraints:

- They should inherit from `GModelFit.AbstractComponent`;

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
- [`GModelFit.dependencies`](@ref): to specify the list of the component dependencies.  If not defined, the fallback method returns `Symbol[]`;
- [`GModelFit.prepare!`](@ref): to pre-compute quantities depending on the evaluation domain.  If not defined, the fallback method does nothing;



The following example shows how to define a custom component and implement all of the above methods:

```@example abc
using GModelFit
import GModelFit: prepare!, dependencies, evaluate!

struct MyComponent <: GModelFit.AbstractComponent
    param1::GModelFit.Parameter

    function MyComponent(param1)
        println(" -> call to MyComponent constructor;")
        new(GModelFit.Parameter(param1))
     end
end

function dependencies(comp::MyComponent)
    println(" -> call to dependencies()")
    return Symbol[]
end

function prepare!(comp::MyComponent, domain::AbstractDomain)
    println(" -> call to prepare!()")
end

function evaluate!(comp::MyComponent, domain::AbstractDomain, output::Vector, param1)
    println(" -> call to evaluate!() with parameter value: ", param1)
    output .= param1
end

println() # hide
```


!!! note
    In the vast majority of cases there is no need to extend the `dependencies` and `prepare!` functions. The only mandatory implementaion is the one for `evaluate!`.


## Life cycle of a component

The life cycle of a component is as follows:

1. The component is created by invoking its constructor;

1. In order to prepare the data structures for component evaluation the following functions are invoked:
   - `GModelFit.dependencies`;
   - `GModelFit.prepare!`;

1. Finally the `GModelFit.evaluate!` method is invoked to actually evaluate the component.


A quick way to evaluate a component on a `Domain` is as follows:
```@example abc
# Create a domain and a component
dom = Domain(1:4)
comp = MyComponent(1)

# Evaluate component on the domain
comp(dom)
nothing  # hide
```

It is also possible to evaluate the component specifying a custom value for the component:
```@example abc
comp(dom, param1=4.5)
nothing  # hide
```

A similar life cycle is observed when the componens is evaluated from within a `Model`:
```@example abc
# Create a domain and a model
dom = Domain(1:4)
model = Model(:mycomp => MyComponent(1))

# Evaluate model on the domain
model(dom)
nothing  # hide
```

It is possible to modify the `param1` and re-evaluate with:
```@example abc
model[:mycomp].param1.val = 4.5
model(dom)
nothing  # hide
```

Finally, when fitting a model the `prepare!` is invoked only once.  The `dependencies` function is invoked a number of times to create the internal data structures, while `evaluate!` is invoked each time the solver needs to probe the model on a specific set of parameter values:
```@example abc
bestfit, fsumm = fit(model, Measures(dom, [9., 9., 9., 9.], 1.))
nothing  # hide
```


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
import GModelFit: prepare!, evaluate!

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
bestfit, fsumm = fit(model, data)
show((bestfit, fsumm)) # hide
```
