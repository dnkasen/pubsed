===========================
Input  Files
===========================

Several example setups for different sorts of science runs are given in the ``examples/`` directory. Simpler tests are in the ``tests/`` directory.

|sedona| everywhere uses cgs units.



.. _parameter_files:

----------------------------------
Parameter Files
----------------------------------

Runtime parameters are set in a parameter file that by default is assumed to be
called **param.lua**

All runtime parameters are collected in the appendix, as well as described throughout
this documentation.

The parameter files uses the Lua scripting language, which allows for
using math expressions and function calls (may be replaced with a python option).
One can also grab environment variables using, e.g.,::

  sedona_home   = os.getenv('SEDONA_HOME')

which sets the local variable **sedona_home** based on the environment variable **$SEDONA_HOME**.

Default values for all parameters are set in a **defaults** file, which is also in
the Lua language. The **param.lua** file must point to the defaults file using, e.g., ::

  defaults_file    = sedona_home.."/defaults/sedona_defaults.lua"

This points to the standard defaults file provided with |sedona|. (N.B. the .. notation means string concatenation in Lua)

Scalar parameters are set as e.g.,::

  tstep_time_stop  = 100.0

String parameters (such as filenames) are set using quotes, e.g.,::

  transport_module = "monte_carlo"

Vectors parameters are currently set as e.g.,::

  transport_nu_grid   = {start,stop,delta}

where gives a uniform vector between values **start** and **stop** with spacing **delta**.
To use a logarithmically spaced grid, we add an extra entry of 1::

  transport_nu_grid   = {start,stop,delta,1}

where now the spacing between points is dx = x*delta




.. _model_files:

----------------------------------
Model Files
----------------------------------

Model file set the density, composition, temperature, etc...  of the initial model.

.. list-table:: Model File Parameters
        :widths: 15,10,40
        :header-rows: 1

        * - parameter
          - values
          - definition
        * - grid_type
          - "grid_1D_sphere", "grid_2D_cyln", "grid_3D_cart"
          - grid geometry; must match input model
        * - model_file
          - <string>
          - Name of model file




---------------------------
Atomic Data Files
---------------------------

Atomic datafiles hold detailed information about atoms, and
can be found in the **data/** folder.

Essential atomic data is compiled into a single hdf5 file. Additional atomic line data
(e.g., Kurucz line lists) can optionally be accessed using "fuzzline'' files

.. list-table:: Atomic Data Files
        :widths: 15,10,40
        :header-rows: 1

        * - parameter
          - values
          - definition
        * - data_atomic_file
          - <string>
          - name of the atomic data file
        * - data_fuzzline_file
          - <string>
          - name of fuzzline file to include extra "fuzz" lines
