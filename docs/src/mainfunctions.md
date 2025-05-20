```@setup abc
include("setup.jl")
```

# Main functionalities

- *Preparation of empirical data*: both the data domain and empirical values (with associated uncertainties) should be wrapped into `Domain` (or `CartesianDomain`) and `Measures` objects respectively.  Such objects are created by simply passing `AbstractVector{<: Real}` to their respective constructors, e.g.:
  ```@example abc
  using GModelFit
  dom  = Domain([0.1, 1.1, 2.1, 3.1, 4.1])
  data = Measures(dom, [6.29, 7.27, 10.41, 18.67, 25.3],
                        [1.1,  1.1,   1.1,   1.2, 1.2])
  println() # hide
  ```

- *Component creation*: a component is a structure inheriting `GModelFit.AbstractComponent` and hosting one or more fields with type [`GModelFit.Parameter`](@ref).  It is created by simply invoking its constructor, e.g.:
  ```@example abc
  using GModelFit
  
  # Create a stand-alone component (i.e. a component used outisde a model)
  comp = GModelFit.Gaussian(1, 0, 1) # numbers represent the parameter values
  
  # A stand-alone component can be evaluated on a user provided domain as follows:
  comp(Domain(-4:0.1:4))
  
  # Evaluate the component providing custom parameter values:
  comp(Domain(-4:0.1:4), center=0.1, sigma=1.3)
  println() # hide
  ```
  The list of available components is available in [Built-in components](@ref).

  Note: a component with dependencies can't be evaluated as a stand-alone since it requires the corresponding dependencies to be available in a model.


- *Model definition* and *manipulation*: a [`Model`](@ref) object is essentially a dictionary of components with `Symbol` keys.  The `keys()`, `haskey()` and `iterate()` methods defined for the `Model` object provide the usual functionalities as for any dictionary.  .  A model object can be created and manipulated as follows:
  ```@example abc
  using GModelFit
  
  # Create an empty model
  model = Model()
  
  # Add a two Gaussian components, and a third one representing their sum
  model[:comp1] = GModelFit.Gaussian(1, 3, 1)
  model[:comp2] = GModelFit.Gaussian(0.5, 4, 0.3)
  model[:sum] = @fd (comp1, comp2) -> comp1 .+ comp2
  
  # Modify a parameter value:
  model[:comp1].center.val = 5
  
  # Evaluate the model on a user defined domain
  dom = Domain(0:0.1:10)
  model(dom)
  
  # Evaluate the model, but retrieve the outcome of the :comp2 component
  model(dom, :comp2)
  println() # hide
  ```

- *Mock data*: the [`GModelFit.mock()`](@ref) function allows to generate mock data set(s) using a (multi-)model as ground truth, and add a random noise to simulate the measurement process.  An example using the previously defined model and domain is as follows:
  ```@example abc
  data = GModelFit.mock(Measures, model, dom)
  println() # hide
  ```
  This functionality is used in the examples of the next sections to generate the mock datasets.

- *Fitting*: the main functions to fit a model (represented by a [`Model`](@ref) object) to an empirical dataset (represented by a [`Measures`](@ref) object) are [`fit`](@ref) and [`fit!`](@ref).  The latter provide the same functionality as the former with the only difference that upon return the `Model` object will have their parameters set to the best fit values.  In both cases the `Model` object will be evaluated on the same domain associated with the `Measures` object.  An overview of the fit workflow is as follows:

  ![](assets/schema.svg)

  The following code shows how to fit the previously generated mock data set to the above model:
  ```@example abc
  bestfit, stats = fit(model, data)
  ```

  The [`fit`](@ref) function returns a tuple with:
  - a [`GModelFit.ModelSnapshot`](@ref) structure containing a snapshot of the best fit model;
  - a [`GModelFit.FitSummary`](@ref) structure containing statistics on the fit.

  To perform a [Multi-dataset fitting](@ref) simply pass a `Vector{Model}` and a `Vector{Measures` to the `fit` function.

- *Serialization*: a few structures (such as  [`GModelFit.ModelSnapshot`](@ref), [`GModelFit.FitSummary`](@ref) and [`Measures{N}`](@ref)) can be *serialized*, i.e. stored in a file, and later *de-serialized* in a separata Julia session.  This is useful when the best fit model and associated informations must be saved for a later use, without the need to re-run the fitting.  The best fit model, fit statistics and mock dataset used above can be serialized with:
  ```@example abc
  GModelFit.serialize("my_snapshot.json", bestfit, stats, data)
  println() # hide
  ```
  In a separate Julia session, you can obtain a copy of exactly the same data with
  ```@example abc
  using GModelFit
  (bestit, stats, data) = GModelFit.deserialize("my_snapshot.json")
  println() # hide
  ```
