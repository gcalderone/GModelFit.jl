# Version 0.3.2
    * Updated docstrings and documentation;

	* Internals: ModelEval is now able to track the changes in the original Model after it has been created;


# Version 0.3.1
    * Updated docstrings

    * Removed unused dependency


# Version 0.3.0
- Breaking Changes:
    * The `Model` constructor no longer accepts a `Domain` argument;

    * The `mock()` function now requires a `Domain` argument;

    * The `fit()` function no longer modifies the provided `Model` argument. The old behaviour is available with the newly added `fit!()` function;

    * The `@λ` macro has been renamed to `@fd` to remind that the resulting value is a `FunctDesc` structure;

    * Domain objects can no longer be indexed as if they were vectors.  The same functionality is available via the `coords()` or `axis()` functions;

    * The minimizer status is no longer an `Enum`, the same information is now indicated by the corresponding subtypes of `AbstractMinimizerStatus`;


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
