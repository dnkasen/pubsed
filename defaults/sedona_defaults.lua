-----------------------------------------------------
-- This file sets the defaults for all parameters
-- used in the sedona radiation transport code
-----------------------------------------------------

-- atomic data files
data_atomic_file   = "../../data/cmfgen_atomdata.hdf5"
data_fuzzline_file = ""

-- grid
grid_type      = "grid_1D_sphere"  -- grid geometry; must match input model

-- job run (i.e., checkpoint/restart) parameters
run_do_restart              = 0
run_do_checkpoint           = 0
run_do_checkpoint_test      = 0
run_restart_file            = "chk.h5"
run_checkpoint_name_base    = "chk"
run_chk_timestep_interval   = 0
run_chk_walltime_interval   = 0. 
run_chk_simtime_interval    = 0.
run_chk_walltime_max_buffer = 0.
run_chk_walltime_max       = 0.
run_chk_number_start        = 0 -- number of first checkpoint file

-- default hydro module is none
hydro_module      = "none"
hydro_gamma_index = 5.0/3.0
hydro_mean_particle_mass = 0.7
hydro_cfl         = 0.1
hydro_v_piston    = 0
hydro_viscosity_parameter = 5
hydro_central_point_mass  = 0
hydro_use_gravity    = 0
hydro_use_transport  = 0
hydro_accrete_radius = 0
hydro_bomb_radius    = 0
hydro_bomb_energy    = 0
hydro_boundary_outflow = 1
hydro_boundary_rigid_outer_wall = 0


-- default nu_grid is single frequency
transport_module = "monte_carlo"
transport_nu_grid  = {1,1,1}
transport_radiative_equilibrium  = 0
transport_steady_iterate         = 0
transport_boundary_in_reflect    = 0
transport_boundary_out_reflect   = 0
transport_store_Jnu              = 1
transport_use_ddmc               = 0
transport_ddmc_tau_threshold     = 100
transport_fleck_alpha            = 0

-- whether or not to fix RNG seed
transport_fix_rng_seed           = 0
-- what a fixed RNG seed should be set to
transport_rng_seed               = 0

-- options for handling temperature solve
transport_solve_Tgas_with_updated_opacities = 0
transport_fix_Tgas_during_transport         = 0
transport_set_Tgas_to_Trad                  = 0

-- inner source emission = none
core_n_emit           = 0
core_radius           = 0
core_luminosity       = 0
core_temperature      = 0
core_photon_frequency = 0
core_timescale        = 0
core_spectrum_file    = ""
core_fix_luminosity   = 0

-- default particle params
particles_max_total = 1e7
particles_n_emit_radioactive   = 0
particles_n_emit_thermal       = 0
particles_n_initialize         = 0
particles_n_emit_pointsources  = 0
particles_pointsource_file     = ""
particles_last_iter_pump        = 1
multiply_particles_n_emit_by_dt_over_dtmax = 0
force_rprocess_heating         = 0

-- time stepping
tstep_max_steps  = 1000
tstep_time_start = 0.0
tstep_time_stop  = 100.0
tstep_max_dt     = 1000.0
tstep_min_dt     = 0.0
tstep_max_delta  = 0.1

-- output parameters
output_write_atomic_levels        = 0
output_write_radiation            = 0
output_write_plt_file_time        = 1
output_write_plt_log_space        = 0  -- use logarithimic spacing of write times
output_write_mass_fractions       = 0

-- limiting values for calculation
limits_temp_max = 1e8
limits_temp_min = 1000




-- opacity calculation defaults
opacity_grey_opacity				= 0
opacity_zone_dependent_grey_opacity	= 0
opacity_user_defined		        = 0
opacity_epsilon             		= 1.0
opacity_atom_zero_epsilon  			= {}
opacity_electron_scattering 		= 0
opacity_line_expansion      		= 0
opacity_fuzz_expansion      		= 0
opacity_bound_free          		= 0
opacity_bound_bound         		= 0
opacity_free_free           		= 0
opacity_use_nlte            		= 0
opacity_atoms_in_nlte       		= {}
opacity_use_collisions_nlte           	= 1 -- only matters if use_nlte == 1
opacity_no_ground_recomb                = 0
opacity_minimum_extinction  		= 0
opacity_maximum_opacity     		= 1e40
opacity_no_scattering       		= 0
dont_decay_composition      		= 0

opacity_compton_scatter_photons = 0;

-- line treatment parameters
line_velocity_width         = 0
line_profile                = "voigt"
line_x_extent               = 100


-- output spectrum information
spectrum_name      = "spectrum";
spectrum_time_grid = {0,1,1}
spectrum_nu_grid   = {1,1,1}
spectrum_n_mu      = 1
spectrum_n_phi     = 1
spectrum_suppress_txt = 0

-- output gamma-ray spectrum
gamma_name     = ""
gamma_nu_grid  = {1,1,1}

spectrum_particle_list_name = "escaped_particles"
spectrum_calc_particle_files = {}
spectrum_calc_file_to_rank = {}
spectrum_calc_chk_file = "chk.h5"
spectrum_calc_out_file = "spectrum_out"
spectrum_calc_n_ranks_old = 1
