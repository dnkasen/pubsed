======================
Output Files
======================


----------------------------------
Spectrum Files
----------------------------------

Files with the name **spectrum_?.h5** carry the emergent radiation of
the model. Output is always in frequency space, with the output luminosity
in units of ergs/sec/Hz.

Note that the time spacing and frequency spacing of the output spectrum
need not be identical to that used in the transport calculation. However,
it is sensible to use time and frequency resolutions that are comparable to
or coarser than that used in the transport run.



To set the frequency grid of the output spectrum, set the runtime parameter::

  spectrum_nu_grid   = {start,stop,delta}

where gives a uniform frequency grid between values **start** and **stop** with spacing **delta**.
To use a logarithmically spaced frequency grid,  add an extra entry of 1::

  spectrum_nu_grid   = {start,stop,delta,1}

where now the spacing between points is dnu = nu*delta.  To perform
a single frequency (i.e., grey) calculation, simply set start >= stop, e.g.,::

  spectrum_nu_grid   = {1,1,1}

The spacing of the time bins is set in a similar way, e.g.,::

  spectrum_time_grid   = {start,stop,delta}

which gives a uniformly spaced time gridbetween values **start** and **stop** with spacing **delta**.



Parameters controlling spectrum output

.. list-table:: Output Spectrum Parameters
        :widths: 15,10,40
        :header-rows: 1

        * - parameter
          - values
          - definition
        * - spectrum_name
          - <string>
          - name of the output spectrum files, if output_write_radiation is enabled
        * - spectrum_time_grid
          - <float vector>
          - time grid for the spectrum file
        * - spectrum_nu_grid
          - <float vector>
          - frequency grid for the spectrum file
        * - spectrum_n_mu
          - <integer>
          - number of evenly spaced mu (viewing angles in theta (polar coord.) direction mu = cos theta)
        * - spectrum_n_phi
          - <integer>
          - number of evenly spaced phi (viewing angles in phi (polar coord.) direction)
        * - gamma_name
          - <string>
          - name of the output gamma-ray spectrum, if radioactivity is being used
        * - gamma_nu_grid
          - <float vector>
          - grid for output gamma-rays; dimensions here are MeV

----------------------------------
Plt Files
----------------------------------

Plt files contain data describing the physical properties (e.g.,
temperature, density) of the model
and different time steps or iterations. Exhaustive information
is contained in the **plt_?????.h5** files, in hdf5 form.  For 1D
calculations, some information is also output (for convenience) in
ascii **plt_?????.dat** files.

Several runtime parameters control when and what is written to
plt files.  To control when plt files are output, set e.g., ::

  output_write_plt_file_time = 1000.0

which will write a plt file every 1000.0 seconds of simulation time.
To write plt files out with logarithmic spacing use e,g.,::

  output_write_plt_log_space = 0.5

In which case, if a plt file is written at time t0, the next
plt file will be written at time t1 = t0*(1 + 0.5).  This parameter
will override the **output_write_plt_file_time** parameter


.. list-table:: plt File Output Parameters
        :widths: 15,10,50
        :header-rows: 1

        * - parameter
          - values
          - definition
        * - output_write_plt_file_time
          - <float>
          - interval of simulation time (in seconds) before writing next plt file
        * - output_write_plt_log_space
          - <float>
          - using logarithmic spacing for plt file output. If equal to 0, use
            equal spacing set by output_write_plt_file_time. If > 0 will
            override write_plt_file_time
        * - output_write_radiation
          - 0 = no | 1 = yes
          - Write out frequency dependent radiation properites (e.g., opacity, emissivity, Jnu)
            for every zone
        * - output_write_atomic_levels
          - 0 = no | 1 = yes
          - Write out detailed level populations for every zone
        * - output_write_mass_fractions
          - 0 = no | 1 = yes
          - Write out the composition (mass fractions) for every zone




----------------------------------
Checkpoint Files
----------------------------------

Checkpoint files provide a complete dump of the state of the
program, and are used to restart a calculation.  By default, checkpointing
will not happen.

Several parameters are available to control checkpointing.  For example, setting::

  run_do_checkpoint = 1
  run_chk_timestep_interval = 10

will write a checkpoint file every 10 time steps.  Meanwhile setting::

  run_do_checkpoint = 1
  run_chk_simtime_interval = 1000.0

will write a checkpoint file every 1000.0 seconds of time elapsed in the simulation.  To
make sure that a checkpoint is written before a job is scheduled to die, use::

  run_do_checkpoint = 1
  run_chk_walltime_max = 12*60.0*60.0
  run_chk_walltime_max_buffer = 1.1

will make sure that a checkpoint is written when the code thinks the next
time step will not complete within the set max walltime of 12 hours. This estimated
by determining if the time left before 12 hours is less than the time it took to compute
the last time step, multiplied by a buffer (here = 1.1) for safety.



.. list-table:: Checkpoint Parameters
        :widths: 15,10,80
        :header-rows: 1

        * - parameter
          - values
          - definition
        * - run_do_restart
          - 0 = no | 1 = yes
          - Whether or not to restart from a checkpoint file. If 0, starts a fresh run. Otherwise, restarts from run_restart_file.
        * - run_restart_file
          - <string>
          - Name of file to restart from (e.g., chk.h5)
        * - run_do_checkpoint
          - 0 = no | 1 = yes
          - Whether or not to writeout checkpoint files. Note, that one of the interval
            parameters below must also be specified to write checkpoints
        * - run_checkpoint_name_base
          - <string>
          - Filename prefix for checkpoint files
        * - run_chk_timestep_interval
          - <int>
          - If 0, don't checkpoint based on simulation iteration number. Otherwise, checkpoint every $run_chk_timestep_interval timesteps.
        * - run_chk_walltime_interval
          - <float>
          - If 0, don't checkpoint based on wallclock time. Otherwise, checkpoint $run_chk_walltime_interval after the last checkpoint in wallclock time. Measured in seconds
        * - run_chk_simtime_interval
          - <float>
          -  If 0, don't checkpoint based on simulation time. Otherwise, checkpoint $run_chk_simtime_interval after the last checkpoint in simulation time. Measured in seconds,
        * - run_chk_walltime_max
          - <float>
          - If 0, don't checkpoint based on when the simulation will end. Otherwise, checkpoint when the simulation thinks it might not finish before $run_chk_walltime_max
            of wallclock time has elapsed since the start of the run. Checkpoints based on this condition happen when ${run_chk_walltime_max_buffer} *
            (walltime duration of last timestep) + (current walltime) >= ${run_chk_walltime_max}. Measured in seconds, default is 0. This time should probably be the wallclock limit on your run.
        * - run_chk_walltime_max_buffer
          - <float>
          - See above. Default is 1.1. Setting this to 0 will also turn off checkpointing based on run_chk_walltime_max
        * - run_chk_number_start
          - <int>
          -  Number with which to start checkpoint file numbering.
        * - run_do_checkpoint_test
          - 0 = no | 1 = yes
          - Whether to save out a checkpoint file immediately after reading in a restart file. If you choose to run this test, running h5diff on the
              restart file and this initial checkpoint file (named {$run_checkpoint_name_base}_init.h5) should return empty.
