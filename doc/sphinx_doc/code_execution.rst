=================================
Code Execution
=================================

|sedona| can be run in one of two modes:


* A :ref:`steady_state` calculation solves the radiation transport for a time-independent "snapshot" in time, whith the matter density and velocity profiles assumed to be fixed. If the radiation itself modifies the matter temperature or ionization/excitation profile, then multiple iterations can be used to converge the solution.

* A :ref:`time_dependent` calculation evolves the transport in time according to user defined time-stepping criteria.

In both time dependent and steady state calculations, the velocity of the gas influences the transport calculation by including Doppler shifts, abberation, adiabatic losses.


.. _steady_state:
----------------------------------
Steady State Calculations
----------------------------------

To run a time-independent steady-state calculation set with (for example) 4 iterations, set::

  transport_steady_iterate = 4

The code will calculate radiation transport assuming the matter distribution (i.e., the density and velocity) is fixed.   This should be a reaonable approximation when the timescale for photons to diffuse through the ejecta is much faster than the time it takes for the  atter distribution to change significantly.
All Monte Carlo photon packets will be propagated through the matter until they either escape or are destroyed.

A single steady-state *iteration* consists of calculating the matter opacities/emissivites, solving the radiation transport, and then updating the matter temperature structure. Multiple iterations are generally required to converge to the  correct temperature structure. To caclulate the temperature under the commonly applied assumption of *radiative equilibrium* (i.e., when the heating by absorption of photons is balanced by the cooling by the emission of photons) set::

  transport_radiative_equilibrium  = 1

Radiative equilibrium is a reasonable assumption when the timescale for matter heating and cooling is short enough that the gas has time to reach equilbrium. When the radiative equilibrium parameter is set to be non-zero, the temperature at each iterate will be set to be.

When radiative equilibrium is set, photon packets that are absorbed are not removed from the calculation.  Insteady the are "effectively scattered", i.e., they are immediately re-emitted with a different direction and a different frequency (sampled from the local emissivity).  This enforces the radiative energy balance (heating = cooling) in the transport propagation.

See Example Problem for a specific example

.. _time_depdendent:
----------------------------------
Time Dependent Calculations
----------------------------------

A time-dependent |sedona| calculation evolves a system from a user-defined start time to stop time through a series of discrete time steps. The time duration of the run is set in the :ref:`parameter_files`; for example, to run a calculation that starts at 100 seconds and evolves to 1000 seconds, set::

  tstep_time_start = 100
  tstep_time_stop  = 1000

The size of each time step is set relative to the current time. For example, to set the time step duration to be 10% of the current simulation time, set::

  tstep_max_delta = 0.1

To limit the maximum of minumum time step size, set::

  tstep_min_dt  = 1
  tstep_max_dt  = 20

Other constraints on the time step may be applied from other physics in the code. For example, if hydrodynamics is being set, the time step size will be limited by a courant condition.

Unless otherwise specified, a |sedona| calculation will calculate the time-dependent radiation transport under the assumption that the matter distribution is fixed with time. If desired, the matter distribution can also be evolved (in limited ways) by setting the ``hydro_module`` parameter.  For example, setting::

  hydro_module = "homologous"

Will assume the gas velocity profile is homologous (i.e., at all points, velocity is proportional to radius, :math:`r = v t` where :math:`t` is the time since homology set it). Such a velocity profile is often applied modeling supernova ejecta.  In a homologous calculation, the radius of the matter scales out linearly with time and the density drops as :math:`\rho \propto t^{-3}`.

The model file has a start time.

For general velocity profiles, evolving the matter distribution requires solving the equations of hydrodynamics. The current |sedona| distribution has limited support for radiation-hydrodynamics. In 1D spherical geometries, one can set::

  hydro_module = 1D_lagrangian

which calculates hydrodynamics using 1D langragian method with an artificlal viscosity. Parameters for controlling hydrodynamics are given below.


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
          - Maximum fractional size of a timestep, restricts dt to the specified value multiplied by the current time

Times are always in seconds.
