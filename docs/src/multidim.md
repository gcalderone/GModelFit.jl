```@setup abc
include("setup.jl")
```

# Multi-dimensional fitting

**GModelFit.jl** supports multi-dimensional data fitting, namely the case where both the model and the empirical measures are expressed as a function of an `N`-dimensional domain.

## Example in 2D

```@example abc
using GModelFit

model = Model(:background => GModelFit.OffsetSlope(0, 0, 0., 0.5, 1.),
              :psf => GModelFit.Gaussian(100., 0., 0., 1, 0.3, 15), 
			  :main => SumReducer(:background, :psf))

dom = CartesianDomain(-5:0.25:5, -4:0.25:4)
data = GModelFit.mock(Measures, model, dom)

bestfit, fsumm = fit(model, data)
show((bestfit, fsumm)) # hide
```

the results can be displayed with [Gnuplot.jl](https://github.com/gcalderone/Gnuplot.jl) as follows:
```@example abc
using Gnuplot

# Plot the model...
@gsp "set palette" axes(dom, 1) axes(dom, 2) model(dom) "w pm3d notit"

# ... and the data
@gsp :- axes(dom, 1) axes(dom, 2) values(data) "w dots notit lc rgb 'gray'"
saveas("ex2d") # hide
```
![](assets/ex2d.png)



## Example in 3D

```@example abc
using GModelFit

dom = CartesianDomain(-5:0.1:5, -5:0.1:5, -5:0.1:5)
model = Model(:source => @fd (x, y, z, cx=0., cy=0., cz=0.) -> @. sqrt((x-cx)^2 + (y-cy)^2 + (z-cz)^2))
data = GModelFit.mock(Measures, model, dom)
bestfit, fsumm = fit(model, data)
show((bestfit, fsumm)) # hide
```

and plot the best-fit model with:
```@example abc
using Gnuplot

m = bestfit()
@gp "set autoscale fix" "set key out" "set size ratio -1" xlab="X" ylab="Y" :-
for i in 1:size(m)[3]
	@gp :- i m[i, :, :] "w image t 'Z=$(string(axes(dom, 3)[i]))'" cbr=extrema(m) :-
end
@gp
Gnuplot.save("assets/animation.gif", term="gif animate size 600, 360 delay 5") # hide
```
![](assets/animation.gif)
