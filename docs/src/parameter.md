```@setup abc
using Gnuplot
Gnuplot.quitall()
mkpath("assets")
Gnuplot.options.term = "unknown"
empty!(Gnuplot.options.init)
push!( Gnuplot.options.init, linetypes(:Set1_5, lw=2.5, ps=1.5))
saveas(file) = save(term="pngcairo size 550,350 fontscale 0.8", output="assets/$(file).png")
```

# Parameter constraints

Models are characterized by *parameters* (see [Basic concepts and data types](@ref)) whose values are modified during fitting until a convergence criterion is met, and the *best fit* values are identified.  In many cases, however, the parameters can not vary arbitrarily but should satisfy some constraints for their values to be meaningful.  **GFit.jl** supports the definition of constraints by fixing the parameter to a specific value, limiting the value in a user defined range, or by dynamically calculating its value using a mathematical expression and the other parameter values. In the latter case the parameter is said to be *patched* through a mathematical relationship. Finally, **GFit.jl** also supports the definition of parametrized patch expressions, involving both the patched parameter as well as other ones.

An important concept to bear in mind is that the [`GFit.Parameter`](@ref) structure provides two field for the associated numerical value:
- `val`: is the parameter value which is being varied by the minimizer during fitting.  The value set before the fitting is the *guess* value.  The value after fitting is the *best fit* one;
- `actual`: is the result of the patch expression evaluation, and the actual value used when evaluating a component.  Note that this value will be overwitten at each model evaluation, hence setting this field has no effect. The `val` and `actual` values are identical if no patch constraint has been defined.

A parameter constraint is defined by explicitly modifiying the fields of the corresponding [`GFit.Parameter`](@ref) structure. More specifically:
1. to set a parameter to a specific value: set the `val` field to the numeric value and set the `fixed` field to `true`;
1. to set a parameter value range: set one or both the `low` and `high` fields (default values are `-Inf` and `+Inf` respectively);
1. to constraint a parameter to have the same numerical value as another one with the same name (but in another component): set the `patch` value to the component name (it must be a `Symbol`).  In this case the parameter is assumed to be fixed;
1. to dynamically calculate an `actual` value using a mathematical expression depending on other parameter values: set the `patch` field to a λ-function (generated with the [`@λ`](@ref) macro) which must accept a single argument (a dictionary of components) and return a scalar number.  In this case the parameter is assumed to be fixed;
1. to use a parametrized patch expression simply create a a λ-function with two arguments: the first has the same meaning as in the previous case, and the second will be the parameter value.  Note that in this case the patched parameter loses its original meaning, and becomes the parameter of the patch expression.
1. if the patch expression involves parameters from another model in a [Multi-dataset fitting](@ref) scenario simply use `mpatch` in place of `patch`, and the first argument to the λ-function will be a vector with as many elements as the number of models in the [`MultiModel`](@ref) object.

The following examples show how to define constraints for each of the afore-mentioned cases.

### Example

We will consider a model for a 1D domain consisting of the sum of a linear background component (named `bkg`) and two Gaussian-shaped features (`l1` and `l2`):
```@example abc
using GFit

dom = Domain(0:0.1:5)
model = Model(dom, :bkg => GFit.OffsetSlope(1, 1, 0.1),
                   :l1 => GFit.Gaussian(1, 2, 0.2),
                   :l2 => GFit.Gaussian(1, 3, 0.4),
                   :main => SumReducer(:bkg, :l1, :l2))
println() # hide
```

Assume that, for the model to be meaningful, the parameters should satisfy the following constraints:
- the `bkg` should have a fixed value of 1 at `x`=1, and a slope which is in the range [0:0.2]:
```@example abc
model[:bkg].offset.val = 1
model[:bkg].offset.fixed = true

model[:bkg].slope.low  = 0
model[:bkg].slope.high = 0.2
println() # hide
```
- the normalization of `l1` and `l2` must be the same:
```@example abc
model[:l2].norm.patch = :l1
println() # hide
```
- the width of `l2` must be twice that of `l1` (patched parameter):
```@example abc
model[:l2].sigma.patch = @λ m -> 2 * m[:l1].sigma
println() # hide
```
- the center of `l2` must be at a larger coordinate than the center of `l1`.  In this case we create a parametrized patch expression where the parameter is reinterpreted as the distance between the two centers:
```@example abc
model[:l2].center.patch = @λ (m, v) -> v + m[:l1].center
model[:l2].center.val = 1   # guess value for the distance between the centers
model[:l2].center.low = 0   # ensure [l2].center > [l1].center
println() # hide
```

We can fit the model against a mock dataset (see [Generate mock datasets](@ref)):
```@example abc
data = GFit.mock(Measures, model)
res = fit!(model, data)
```
and plot the results with [Gnuplot.jl](https://github.com/gcalderone/Gnuplot.jl):
```@example abc 
using Gnuplot
@gp    coords(dom) values(data) uncerts(data) "w yerr t 'Data'" :-
@gp :- coords(dom) model() "w l t 'Model'"
saveas("example_patch1") # hide
```
![](assets/example_patch1.png)
