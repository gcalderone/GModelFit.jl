```@setup abc
include("setup.jl")
```

# Instrument response


An instrument response allows to convert the *unfolded* model into a form suitable to be comapred with empirical data.  This may also involve a remapping from the unfolded model domain (namely, the domain where all model components are evaluated) to the folded model domain (namely, the domain where the measurements are evaluated).

**GModelFit.jl** provides only one built-in instrument response named `GModelFit.IdealInstrument` and representing an ideal instrument where the folded and unfolded models are identical.  No attempt is made to provide further instrument responses since these are strongly dependent on the instrument being used for the measurements.

It is however possible to implement custom instrument responses as follows:
- Define a new structure inheriting from `GModelFit.AbstractInstrumentResponse`;
- Implement a `model_domain` method whose purpose is to return the domain for the *unfolded* model;
- Implement an `apply_ir!` method whose purpose is to apply the instrument response on an unfolded model evaluation, and to populate the vector of the folded model;


As an example we copy here the implementation for the `IdealInstrument`:
```julia
# Define IdealInstrument structure (no fields are needed)
struct IdealInstrument <: AbstractInstrumentResponse
end

# Method to retrieve the unfolded model domain.  In this case it simply is the same as the folded domain.
model_domain(IR::IdealInstrument, folded_domain::AbstractDomain) = folded_domain

# Method to apply the instrument response.  In this case it simply copies all values from the "unfolded"
# vector to the "folded" one.
function apply_ir!(IR::IdealInstrument,
                   folded_domain::AbstractDomain  , folded::Vector,
                   unfolded_domain::AbstractDomain, unfolded::Vector)
    folded .= unfolded
end
```



A more complex example is as follows: suppose you need to evaluate the unfolded model on an domain which is oversampled with respect to the data domain, and then interpolate the unfolded model into the actual data grid when fitting.

The code to implement such `MyInstrument` instrument domain is as follows:
```@example abc
using GModelFit, Gnuplot
using Dierckx   # used for interpolations

# Import relevant methods
import GModelFit: model_domain, apply_ir!

# Create the new structure
struct MyInstrument <: GModelFit.AbstractInstrumentResponse
    unfolded_domain::GModelFit.AbstractDomain
    MyInstrument(r::StepRangeLen) = new(Domain(r))
end

# Return unfolded model domain 
model_domain(IR::MyInstrument, folded_domain::GModelFit.AbstractDomain) = IR.unfolded_domain

# Apply instrument response
function apply_ir!(IR::MyInstrument,
                   folded_domain::GModelFit.AbstractDomain, folded::Vector,
                   unfolded_domain::GModelFit.AbstractDomain, unfolded::Vector)
    folded .= Dierckx.Spline1D(coords(unfolded_domain), unfolded)(coords(folded_domain))
end
nothing # hide
```

The following example shows how to use the `MyInstrument` instrument response:

```@example abc
dom = Domain([0.5, 2, 3.5, 5.])
data = Measures(dom, [0.596, 0.8, -0.48, -0.98], 0.1)
model = Model(@fd (x, a=1, x0=1, ϕ=0)  -> (a .* sin.(x./x0 .+ ϕ)))
set_IR!(model, MyInstrument(0:0.1:2pi))
bestfit, fsumm = fit(model, data)
```

```@example abc
@gp data bestfit
saveas("ir_example"); # hide
```
![](assets/ir_example.png)
