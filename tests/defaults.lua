
-- atomic data files
data_atomic_file   = "../../data/cmfgen_atomdata.hdf5"
data_fuzzline_file = ""

-- grid
grid_type      = "grid_1D_sphere"  -- grid geometry; must match input model
grid_write_out = 1 -- how often to write out grid data

-- default hydro module is none
hydro_module = "none"

-- default nu_grid is nothing
transport_nu_grid  = {1,1,1}
transport_radiative_equilibrium  = 1
transport_steady_iterate         = 1

-- inner source emission = none
core_n_emit           = 0
core_radius           = 0
core_luminosity       = 0
core_temperature      = 0
core_photon_frequency = 0
core_spectrum_file    = ""

-- default particle params
particles_max_total = 1e6
particles_step_size = 0.3
particles_n_emit_radioactive   = 0
particles_n_initialize         = 0

-- time stepping
tstep_max_steps  = 1000
tstep_time_start = 0.0
tstep_time_stop  = 100.0
tstep_max_dt     = 1000.0
tstep_min_dt     = 0.0
tstep_max_delta  = 0.1

-- limiting values for calculation
limits_temp_max = 1e8
limits_temp_min = 1000

-- opacity calculation defaults
opacity_grey_opacity        = 0
opacity_epsilon             = 1.0
opacity_electron_scattering = 0
opacity_line_expansion      = 0
opacity_fuzz_expansion      = 0
opacity_bound_free          = 0
opacity_bound_bound         = 0
opacity_free_free           = 0
opacity_lines               = 0
opacity_use_nlte            = 0

-- line treatment parameters
line_velocity_width         = 0
line_profile                = "voigt"

-- output spectrum information
spectrum_name      = "spectrum";
spectrum_time_grid = {0,1,1}
spectrum_nu_grid   = {1,1,1}
spectrum_n_mu      = 1
spectrum_n_phi     = 1

-- output gamma-ray spectrum
gamma_name     = ""
gamma_nu_grid  = {1,1,1}



