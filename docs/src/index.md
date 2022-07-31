# Gfit.jl
## A model fitting framework for Julia.
[![Stars](https://img.shields.io/github/stars/gcalderone/GFit.jl?style=social)](https://github.com/gcalderone/GFit.jl)


**Gfit.jl** is a general purpose, data-driven model fitting framework for Julia.

It provides the basic tools to define, interactively manipulate and efficiently evaluate a (possibly very complex) model, and to fit the latter to empirical data. 

The main functionalities are:

- it handles datasets of any dimensionality;
- the fitting model is evaluated on a user defined domain, and is the result of a combination of *components* (either built-in or implemented by the user), or  mathematical expressions (actually [lambda functions](https://en.wikipedia.org/wiki/Anonymous_function)), or any arbitrary mixture of the two;
- model parameters can be fixed to a specific value, limited in an interval, and/or be dynamically linked (patched) to the values of other parameters;
- multiple data sets can be fitted simultaneously against different models whose parameters can optionally be patched;
- all components results are cached so that repeated evaluations with the same parameter values do not involve further calculations;
- user provided components can pre-compute quantities based on the model domain, and store them in reserved areas for re-use;
- it support different minimizers ([LsqFit](https://github.com/JuliaNLSolvers/LsqFit.jl) and [CMPFit](https://github.com/gcalderone/CMPFit.jl));
- it provides several facilities for interactive fitting and result displaying.

The fitting process is actually carried out by one of the above mentioned minimizer packages, whose purpose is to automatically vary the parameter values (according to the user defined constraints) until the differences between the evaluated model and the empirical data are minimized.  The purpose of **Gfit.jl** is thus to act as an interface between the high-level model facing the user, and the low-level implementation details facing the minimizer.


# Installation

In the Julia REPL type:
```julia-repl
julia> ]add GFit
```
The `]` character starts the Julia [package manager](https://julialang.github.io/Pkg.jl/v1/getting-started.html#Basic-Usage-1). Hit backspace key to return to Julia prompt.


In order to easily visualize the outcomes of 1D analysis you may be interested in installing also [GFitViewer](https://github.com/lnicastro/GFitViewer.jl) and/or [Gnuplot.jl](https://github.com/gcalderone/Gnuplot.jl):
```julia-repl
julia> ]add GFitViewer
julia> ]add Gnuplot
```
