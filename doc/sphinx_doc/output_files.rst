======================
Output Files
======================


----------------------------------
Spectrum Files
----------------------------------


Spectrum files are output to hdf5 files named **spectrum_?.h5**, where ? is an integer that indexes the order spectrum files were written. For viewing angle independent calculations (i.e., n_mu = n_phi = 1) files are also output in ascii format for convience, and named **spectrum_?.dat**.


Controlling Spectrum File Output
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Files with the name **spectrum_?.h5** and **spectrum_?.dat** carry the emergent radiation of the model. To set the frequency binning of the output spectrum, set the runtime parameter::

  spectrum_nu_grid   = {start,stop,delta}

where gives a uniform frequency grid between values **start** and **stop** with spacing **delta**.
To use a logarithmically spaced frequency grid,  add an extra entry of 1::

  spectrum_nu_grid   = {start,stop,delta,1}


where now the spacing between points is :math:`\Delta \nu = \nu \times {\rm delta}`.  To perform a single frequency (i.e., grey) calculation, simply set start >= stop, e.g.,::

  spectrum_nu_grid   = {1,1,1}

The spacing of the time bins is set in a similar way, e.g.,::

  spectrum_time_grid   = {start,stop,delta}

which gives a uniformly spaced time grid between values **start** and **stop** with spacing **delta**.

The time and frequency binning of the output spectrum need *not* be identical to that used in the transport calculation. However, it is recommended to use time and frequency resolutions that are comparable to or coarser than that used in the transport run.

In multi-dimensional calculations, the spectrum can be output for different viewing angles.  To set the number of viewing angle bins in the polar angles :math:`(\theta, \phi)` use, for example::

    spectrum_n_phi = 10
    spectrum_n_mu  = 10

where :math:`{\rm mu} = \cos(\theta)`.  The output bins are equally spaced in mu and phi, such that each viewing angle bin is equally likely to be observed.


Spectrum File Format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Spectrum files are written in hdf5 format and (for lower dimensional problmes) as ascii
format. The ascii files begin with a one-line header that gives the dimensions of the output::

  #  n_time    n_freq      n_mu      n_phi

where :math:`n_{\rm time}, n_{\rm freq}, n_{\rm mu}, n_{\rm phi}` are the number
of bins used in each of the time, frequency, mu and phi dimensions.

After the header, the ascii files have columns of data which provide::

  (time)   (frequency)   luminosity     particle_count

The *(time)* column  will be absent in a time-indepdendent (n_times = 1) calculation, while the *(frequency)* column will be absent in a frequncy-independent (n_frequency = 1).
The time and frequency values are given at the center of the bin, and the units are cgs (time in seconds, frequency in Hz)

In a frequency independent calculation, the *luminosity* column gives the bolometric
luminosity (units ergs/sec).  In a frequency dependent calculation, the *luminosity* column gives the observed specific luminosity in the spectral bin, which is
in frequency units  (:math:`L_\nu`, units ergs/sec/Hz). To convert :math:`L_\nu` to the commonly used wavelength units (:math:`L_\lambda`, units ergs/sec/Angstrom) use the conversion formulae:

.. math::

  \lambda = \frac{c}{\nu} \times 10^8 ~~~~~~~~~~~~~~~~~~~~~ L_\lambda = L_\nu  \left( \frac{\nu^2}{c}  \right) \times 10^{-8}

where the factor of :math:`10^8` appears to convert between centimeters to angstroms.

In calculations with multiple viewing angles (n_mu and/or n_phi > 1) the luminosity given is an *isotropic equivalent* luminosity.  That is, it is the total luminosity that would be emitted by an isotropic source with the same spectrum.

The *particle_count* column gives the number of photon packets that were recorded in the given bin. It can be used to roughly estimate the signal to noise in the bin,
:math:`S/N \sim 1/\sqrt{N}`.  However, because different photon packets can have different weights (i.e., energies) this is not a rigorous determination of the S/N.

The hdf5 format psectrum files have a series of datasets that provide the same information. Reading datasets from the file into numpy arrays is simple using the python h5py package::

  import h5py, numpy
  import sedonalib as sed

  fin = h5py.File("spectrum_1.h5","r")
  nu  = numpy.array(fin["nu"])
  Lnu = numpy.array(fin["Lnu"])

The output spectrum is contained in the multi-dimensional Lnu array.

.. list-table:: plt File hdf5 datasets
        :widths: 15,25,40
        :header-rows: 1

        * - **dataset**
          - **type**
          - **description**
        * - Lnu
          - real :math:`[n_{\rm time}][n_{\rm freq}][n_{\rm mu}][n_{\rm phi}]`
          - Output luminosity (erg/sec/Hz or erg/sec)
        * - nu
          - real :math:`[n_{\rm freq}]`
          - frequncy bin center values (units Hz)
        * - nu_edges
          - real :math:`[n_{\rm freq} + 1]`
          - frequency bin edge values

If any of these dimensions is equal to 1, then it will be omitted from the shape of the Lnu array. For example, if
:math:`n_{\rm mu} = n_{\rm phi} = 1` the shape of Lnu will
be ':math:`[n_{\rm time}][n_{\rm freq}]`, or if

----------------------------------
Plt Files
----------------------------------

Plt files contain data describing the physical properties (e.g
.,
temperature, density) of the model
and different time steps or iterations. Exhaustive information
is contained in the **plt_?.h5** files, in hdf5 form.  For 1D
calculations, some information is also output (for convenience) in
ascii **plt_?.dat** files.


Controlling Plt File Output
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Several runtime parameters control when and what is written to
plt files.  To control when plt files are output, set e.g., ::

  output_write_plt_file_time = 1000.0

which will write a plt file every 1000.0 seconds of simulation time.
To write plt files out with logarithmic spacing use e,g.,::

  output_write_plt_log_space = 0.5

In which case, if a plt file is written at time t0, the next
plt file will be written at time t1 = t0*(1 + 0.5).  This parameter
will override the **output_write_plt_file_time** parameter


To write out the detailed properties of the radiation and opacity in each zone set::

  output_write_radiation = 1

This will create a

To write more out detailed properties of the gas excitation/ionization state, set::

  output_write_atomic_levels = 1

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





------------------------------------------------
Table of Runtime Parameters Controlling Output
------------------------------------------------


.. list-table::  Spectrum File Output Parameters
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



.. list-table:: plt File Output Parameters
        :widths: 25,10,40
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



.. list-table:: Checkpoint Parameters
        :widths: 30,10,65
        :header-rows: 1

        * - **parameter**
          - **values**
          - **definition**
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
          - Whether to save out a checkpoint file immediately after reading in a restart file. If you choose to run this test, running h5diff on the restart file and this initial checkpoint file (named {$run_checkpoint_name_base}_init.h5) should return empty.
