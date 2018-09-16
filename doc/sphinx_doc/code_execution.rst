=================================
Basic Code Execution
=================================




Several example setups for different sorts of science runs are given in the ``examples/`` directory. Simpler tests are in the ``tests/`` directory.

|sedona| everywhere uses cgs units.


.. _model_files:

----------------------------------
Input Model Files
----------------------------------

Model file set the density, composition, temperature at each 


.. _parameter_files:

----------------------------------
Setting Runtime Parameters
----------------------------------

Runtime parameters are currently set in the param file, by default is assumed to be
called **param.lua** 

The parameter files uses the Lua scripting language, which allows for 
using math expressions and function calls (may be replaced with a python option).
One can also grab environment variables using, e.g.,::

	sedona_home   = os.getenv('SEDONA_HOME')

which sets the local variable **sedona_home** based on the environment variable **$SEDONA_HOME**.

Default values of all parameters are set in a **defaults** file, which is also in
the Lua language. The **param.lua** file must point to the defaults file using, e.g., ::

	defaults_file    = sedona_home.."/defaults/sedona_defaults.lua"

This points to the standard defaults file provided with |sedona|. (N.B. the .. notation means string concatination in Lua)

Scalar parameters are set as e.g.,::

	tstep_time_stop  = 100.0

String parameters (such as filenames) are set using quotes, e.g.,::

	transport_module = "monte_carlo"

Vectors parameters are currently set as e.g.,:

	spectrum_nu_grid   = {start,stop,delta}

where gives a uniform vector between values **start** and **stop** with spacing **delta**.
To use a logarithmically spaced grid, we add an extra entry of 1::

	spectrum_nu_grid   = {start,stop,delta,1}

where now the spacing between points is dx = x*delta



----------------------------------
Controlling Time Evolution
----------------------------------

Calculations in |sedona| can be run either as time evolving or steady-state models.
This is controlled by time-stepping parameters



.. list-table:: Time Stepping Parameters
	:widths: 15,10,40
	:header-rows: 1

	* - parameter
	  - values
	  - definition
   	* - tstep_max_steps
   	  - <integer>
   	  - Maximum number of time steps to take before exiting
   	* - tstep_time_start
   	  - <real>
   	  - Start time (in seconds)
	* - tstep_time_stop
	  - <real>
	  - Stop time (in seconds)
	* - tstep_max_dt
	  - <real>
	  - Maximum value of a time step (in seconds)
	* - tstep_min_dt
	  - <real>
	  - Minimum value of a time step (in seconds)
	* - tstep_max_delta
	  - <real>
	  - Maximum fractional size of a timestep -- restricts dt to the specified value multiplied by the current time

Times are always in seconds.


----------------------------------
Output Files
----------------------------------
