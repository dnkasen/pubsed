-----------------------------------------------------
-- This file sets the defaults for all parameters
-- used in the sedona radiation transport code
-----------------------------------------------------

-- atomic data files
data_atomic_file   = "../../data/cmfgen_atomdata.hdf5"
data_fuzzline_file = ""

-- grid
grid_type      = "grid_1D_sphere"  -- grid geometry; must match input model

-- default hydro module is none
hydro_module      = "none"
hydro_gamma_index = 4.0/3.0
hydro_cfl         = 0.1
hydro_v_piston    = 0
hydro_viscosity_parameter = 5
hydro_central_point_mass  = 0
hydro_use_gravity    = 0
hydro_accrete_radius = 0
hydro_bomb_radius    = 0
hydro_bomb_energy    = 0

-- default nu_grid is nothing
transport_module = "monte_carlo"
transport_nu_grid  = {1,1,1}
transport_radiative_equilibrium  = 0
transport_steady_iterate         = 0
transport_boundary_in_reflect    = 0
transport_boundary_out_reflect   = 0
transport_store_Jnu              = 1

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
output_write_mass_fractions       = 0
output_write_plt_file_time        = 1
output_write_plt_log_space        = 0  -- use logarithimic spacing of write times
output_write_mass_fractions       = 0

-- limiting values for calculation
limits_temp_max = 1e8
limits_temp_min = 1000

-- opacity calculation defaults
opacity_grey_opacity        = 0
opacity_user_defined        = 0
opacity_epsilon             = 1.0
opacity_atom_zero_epsilon   = {}
opacity_electron_scattering = 0
opacity_line_expansion      = 0
opacity_fuzz_expansion      = 0
opacity_bound_free          = 0
opacity_bound_bound         = 0
opacity_free_free           = 0
opacity_use_nlte            = 0
opacity_minimum_extinction  = 0
opacity_maximum_opacity     = 1e40
opacity_no_scattering       = 0

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

-- output gamma-ray spectrum
gamma_name     = ""
gamma_nu_grid  = {1,1,1}



