# Version 0.3.0
- Breaking Changes:
	* `Model` constructor no longer accepts a `Domain` argument;
	
	* `mock` now accepts the `Domain` argument;

	* fix -> fix!

    * getindex(d::Union{Domain, CartesianDomain}, dim::Integer) is no longer available

    * @λ -> @fd
	
    * Minimizer status
	
# Version 0.2.1

- New features:
	* Components can now be directly evaluated on a domain by invoking them as function, e.g.
	```
	comp = GModelFit.Gaussian(1, 0, 1);
	comp(Domain(-5:5))
	```

- Performance improvements:
	* During a fit the `Model.maincomp` is temporarily set to the main component name.  This allow avoiding unnecessay invocations of `find_maincomp()` during model evaluation;

	* Refactored code in GModelFit.PV


- Bugfix:
	* Fixed a bug in `show()` when a component evaluates to NaN;


# Version 0.2.0

- New features:
	* using PrecompileTools to reduce time-to-first-run in Julia v1.9;

    * Implemented `GModelFit.evalcounter()` to retrieve the number of times a component has been evaluated;

	* Implemented `show()` method for `ModelSnapshot` objects;

	* Refactored serialization code;

- Bugfix:
	* Fixed accessibility issue for parameter of `FCompv`;


# Version 0.1.0
- First release.
