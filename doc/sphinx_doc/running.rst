=================
Running Sedona
=================

Before running [sedona] it is helpful to set the environment variable  ``SEDONA_HOME`` to point to the base directory of sedona. In the bash shell, for example, add this line to your ``.bash_profile`` or ``.bashrc``::

	export SEDONA_HOME=/Users/kasen/sedona6/

replacing ``/Users/kasen/sedona6/`` with the base directory of your sedona


To run [sedona] on a single processor, copy the compiled ``sedona6.ex`` executable into the directory where you plan to run.  Make sure the required input files exist (see ... ) and run the code with the command::

  ./sedona6.ex my_parameter_file.lua

where ``my_parameter_file.lua`` is the name of the input parameter file containing runtime parameters.  If you run ``./sedona6.ex`` with no subsequent argument, the code will assume the parameter file is name *param.lua*

For  calculations on multiple processing cores, [sedona] uses a mixture of MPI and openMP parallelism. You should refer to the instructions on the system you are using to determine the appropriate commands to execute a parallel calculation.   If using mpich, for example, the command to run with MPI parallelism is::

  mpirun -n 2 ./sedona6.ex my_parameter_file.lua

where the number after ``-n`` is the number of mpi ranks to use (in this case 2).

To run using openMP  threading, often you set the number of threads as an environment variable.  In the bash shell, for example::

  export OMP_NUM_THREADS=4

where the value of ``OMP_NUM_THREADS`` is the number of threads to use (in this case 4).

A key distinction between the two parallelism models is that openMP threads share the same memory whereas MPI ranks do not. Thus, each MPI rank used in [sedona] allocates a replica of the entire simulation grid and independently transports a subset of the Monte Carlo photon particles, with the results being communicated among all processors at the end of each transport step. In contrast, all openMP threads (on a given MPI rank) share and operate on the same simulation grid.

The ``examples/`` directory if the [sedona] installation provides example setups for a range of different science problems. This is a good place to start becoming familiar with calculations; below we describe a few instructive examples.  If you are acustomed to runnning Jupyter notebooks, the directory ``examples/jupyter_notebooks`` provides some example notebooks.

--------------------------------------
Example #1: Spherical Lightbulb Test
--------------------------------------

The *spherical lightbulb* is a quick test problem consisting of a spherical surface that
uniformly emits blackbody radiation into an empty medium. The radial dependence of the radiation field in steady state can be calculated analytically

.. math::

  T_{\rm rad}(r) = T_{\rm ph} \left[ \frac{1}{2} \left(1 - \sqrt{1 - R_{\rm ph}^2/r^2} \right) \right]^{1/4}

where :math:`R_{\rm ph}` is the radius and :math:`T_{\rm ph}` the blackbody temperature of the emitting spherical surface.  The observed spectrum should be a blackbody of temperature :math:`T_{\rm ph}`.

To run the problem, change to the ``examples/spherical_lightbulb/1D`` directory and copy the executable *sedona6.ex* there. You will see two input files: a parameter file and model file.  See section for description of the files.

Sets core by using core_radius, core_temperature.  Sets the number of photons emitted. Sets a vacuum by using grey opacity with a very small value....

Produces output...

Can change parameters, opacity to see the effects.  

-------------------------------------------
Example #2: Type Ia Supernova Spectrum
-------------------------------------------


----------------------------------------------
Example #3: Type Ia Supernova Light Curve
----------------------------------------------
