===================================
All Runtime Parameters
===================================

|sedona| uses cgs units everywhere

|




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

|




.. list-table:: Hydrodynamics Parameters
        :widths: 15,10,40
        :header-rows: 1

        * - parameter
          - values
          - definition
        * - hydro_module
          - <string>
          - choose hydro module to evolve rho and T, options are homologous, none, 1D_lagrangian
        * - hydro_gamma_index
          - <float>
          - gamma used in the eos
        * - hydro_mean_particle_mass
          - <float>
          - mean particle mass mu, mass = mu * proton mass
        * - hydro_cfl
          - <float>
          - cfl parameter used in hydro simulation to control time steps
        * - hydro_v_piston
          - <float>
          -
        * - hydro_viscosity_parameter
          -
          -
        * - hydro_central_point_mass
          -
          -
        * - hydro_use_gravity
          -
          - Whether to include gravity for the hydro simulation
        * - hydro_use_transport
          - 0 = no | 1 = yes
          - Whether to use radiation transport (?)
        * - hydro_accrete_radius
          -
          -
        * - hydro_bomb_radius
          - <float>
          - Radius of the bomb in cm
        * - hydro_bomb_energy
          - <float>
          - Energy of the bomb in erg
        * - hydro_boundary_outflow
          - 0 = no | 1 = yes
          -
        * - hydro_boundary_rigid_outer_wall
          - 0 = no | 1 = yes
          -


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

|
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



.. list-table:: Opacity parameters
        :header-rows: 1
        :widths: 20,10,40

        * - parameter
          - values
          - definition
        * - opacity_grey_opacity
          - <real>
          - value of grey opacity to use (in cm^2/g). Will override all other opacity settings
        * - opacity_zone_specific_grey_opacity
          - 0 = no | 1 = yes
          - Use a zone dependent grey opacity dataset set in an hdf5 input model file
        * - opacity_user_defined
          - 0 = no | 1 = yes
          - Calculate opacities by calling the function
        * - opacity_epsilon
          - <float>
          - The fraction of
        * - opacity_atom_zero_epsilon
          - <int>
          -
        * - opacity_electron_scattering
          - 0 = no | 1 = yes
          - include electron scattering opacity
        * - opacity_line_expansion
          - 0 = no | 1 = yes
          - include binned line expansion opacity
        * - opacity_fuzz_expansion
          - 0 = no | 1 = yes
          - include binned line expansion opacity, taken from a fuzz file
        * - opacity_bound_free
          - 0 = no | 1 = yes
          - include bound-free (photoionization) opacity
        * - opacity_free_free
          - 0 = no | 1 = yes
          - include free-free opacity
        * - opacity_bound_bound
          - 0 = no | 1 = yes
          - include bound-bound (resolved line) opacity
        * - opacity_use_nlte
          - 0 = no | 1 = yes
          - include nlte opacity
        * - opacity_atoms_in_nlte
          - <int vector>
          - A vector of atomic numbers of the species to be treated in NLTE
        * - opacity_use_collisions_nlte
          - 0 = no | 1 = yes
          - only matters if use_nlte == 1, include collisions for nlte calculations
        * - opacity_no_ground_recomb
          - 0 = no | 1 = yes
          - Suppress all recombination transitions to the ground state in the NLTE level population solve
        * - opacity_minimum_extinction
          - <float>
          - Minimum value of the extinction coefficient (units 1/cm) in any zone
        * - opacity_maximum_opacity
          - <float>
          - Minimum value of the extinction coefficient (units 1/cm) in any zone
        * - opacity_no_scattering
          - 0 = no | 1 = yes
          - if = 1, will not include any kind of scattering opacity
        * - dont_decay_composition
          -
          -
        * - opacity_compton_scatter_photons
          -
          -
        * - line_velocity_width
          - <float>
          - velocity in cm/s used to doppler broaden the (bound bound?) lines
        * - line_x_extent
          -
          -

|

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

|
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

|
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



.. list-table:: Radiation Transport Parameters
        :widths: 15,10,40
        :header-rows: 1

        * - parameter
          - values
          - definition
        * - transport_module
          - "monte_carlo"
          - What method to use for transport. Currently only monte carlo is implemented.
        * - transport_nu_grid
          - <float vector>
          - Define frequency grid used for transport and opacities
        * - transport_radiative_equilibrium
          - 0 = no | 1 = yes
          - Whether to solve for radiative equilibrium
        * - transport_steady_iterate
          - <integer>
          - Do a steady-state calculation with this number of iterations
        * - transport_boundary_in_reflect
          - 0 = no | 1 = yes
          -
        * - transport_boundary_out_reflect
          - 0 = no | 1 = yes
          -
        * - transport_store_Jnu
          - 0 = no | 1 = yes
          -
        * - transport_use_ddmc
          - 0 = no | 1 = yes
          - Whether to use discrete diffusion monte carlo
        * - transport_ddmc_tau_threshold
          - <float>
          - At what optical depth ddmc takes over
        * - transport_fleck_alpha
          - <float>
          - fleck alpha parameter (needs to be between 0.5 and 1 for ddmc)
        * - transport_solve_Tgas_with_updated_opacities
          - 0 = no | 1 = yes
          - whether to solve for Tgas after updating opacities
        * - transport_fix_Tgas_during_transport
          - 0 = no | 1 = yes
          - whether to fix Tgas
        * - transport_set_Tgas_to_Trad
          - 0 = no | 1 = yes
          - whether to set Tgas to Trad instead of solving for it

|


.. list-table:: Radiating Core Parameters
        :widths: 15,10,40
        :header-rows: 1
        
        * - parameter
          - values
          - definition
        * - core_n_emit
          - <integer>
          - Number of particles to emit from core per time step (or iteration)
        * - core_radius
          - <real>
          - Radius (in cm) of emitting spherical core
        * - core_luminosity
          - <real> or <function>
          - Luminosity (in erg/s) emitted from core
        * - core_temperature
          - <real>
          - Blackbody spectrum of core emission, if using blackbody emission
        * - core_photon_frequency
          - <real>
          - Frequency of photons emitted from core, if using monochromatic emission
        * - core_timescale
          - <real>
          - ?
        * - core_spectrum_file
          - <string>
          - filename of file to read to set spectrum of core emission
        * - core_fix_luminosity
          - 0 = no | 1 = yes
          - In steady state calculations, will rescale to fix output luminosity
        * - particles_max_total
          - <float>
          - maximum number of particles (photons) allowed on the grid at the same time
        * - particles_n_emit_radioactive
          - <integer>
          - number of particles emitted through radioactivity per timestep
        * - particles_n_emit_thermal
          - <integer>
          - number of thermal particles emitted per timestep
        * - particles_n_initialize
          - <integer>
          - number of particles used to initialize the simulation
        * - particles_n_emit_pointsources
          - <integer>
          -
        * - particles_pointsource_file
          - <string>
          -
        * - particles_last_iter_pump
          -
          -
        * - multiply_particles_n_emit_by_dt_over_dtmax
          - 0 = no | 1 = yes
          -
        * - force_rprocess_heating
          - 0 = no | 1 = yes
          -


|


