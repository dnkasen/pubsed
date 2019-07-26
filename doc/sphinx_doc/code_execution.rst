=================================
Basic Code Execution
=================================




----------------------------------
Time Dependent Calculations
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
          - Maximum fractional size of a timestep, restricts dt to the specified value multiplied by the current time

Times are always in seconds.


----------------------------------
Steady State Calculations
----------------------------------
