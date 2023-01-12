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
              :psf => GFit.Gaussian(100., 0., 0., 1, 0.3, 15), 
			  :main => SumReducer(:background, :psf))
data = GFit.mock(Measures, model)
ret1 = fit!(model, data)
```

To produce the plots I will use the [Gnuplot.jl](https://github.com/gcalderone/Gnuplot.jl) package, but the user can choose any other package:

```julia
using Gnuplot

# Plot the model...
@gsp coords(dom, 1) coords(dom, 2) model()

# ...and the residuals
@gsp coords(dom, 1) coords(dom, 2) values(data) - model()

# Plot using pm3d style
@gsp "set pm3d" "set palette" coords(dom, 1) coords(dom, 2) model() "w dots"
```
