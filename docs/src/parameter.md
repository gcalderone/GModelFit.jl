# Parameter constraints

A parameter represents a quantity we wish to estimate from the comparison between a model and the empirical data.  During fitting the parameter value is modified until a convergence criterion is reached and a *best fit* value is identified.

Meaningful 

However, not all parameters are independent, or are free to vary on an unlimited range.  **GFit.jl** allows to specify constraint for each parameter.

All model parameters are represented by objects of type [`GFit.Parameter`](@ref).  Note that each parameter has two numerical values associated:
- `val`: is the parameter value which is being varied (in the range `low`:`high`) by the minimizer during fitting. The value set before the fitting is the *guess* value.  The value after fitting is the *best fit* one;
- `actual`: is the actual parameter being used to evaluate the model.  

If no constraint is defined for the parameter the 

The two values are identical, unless a patch prescription make 



containing both the parameter numerical value and the constraint prescriptions. Note that 

.  Such structure contains the parameter value (`val`) as well as other fields specifying the constraint prescriptions


The current value, in particular, is modified during the fitting process to minimize the differences between the model and the data.  The possible values for each parameter can be either unspecified (i.e. the value can be anything from `-Inf` to `+Inf`) or subject to a number of constraints. 

The possible constraints are:
- fixed parameter: the parameter is fixed to a specific value and can not be modified during the fitting process.  To fix a `Parameter` object simply set its `fixed` field to true;
- limited range: the parameter is allowed to vary only within a range limited by the `low` and `high` fields values;
- the current parameter value can be constrained to be equal to that of another parameter with the same name, but belonging to another component within the same model.  To establish such link simply set the component name in the `patch` field;
- the current parameter value can be dynamically calculated by a mathematical expression involving all other parameters within the same model. In this case the `patch` should be set to a `@λ` function with a single argument.  The parameter is au


- model parameters can be fixed to a specific value, limited in an interval, and/or be dynamically linked (patched) to the values of other parameters (see [Parameter constraints](@ref));

