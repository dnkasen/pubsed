=======================================
Adding Runtime Parameters
=======================================

All runtime parameters can be viewed in the defaults/sedona_defaults.lua file

Runtime parameters are blocked into different structures that group together
related parameters (e.g., opacity_ , hydro_, ...)

To avoid an over-abundance of runtime parameters, please think if the functionality
can be accomplished in a way that does not require additional parameters.

The guidelines for adding a new runtime parameter are

* Choose a clear name that begins with its structure name
* Add the parameter into the **sedona_defaults.lua** file with a sensible default value
* Add documentation for the parameter into the appropriate parameter tables in the user's guide
* Check that adding the parameter has not broken any of the test suite problems
