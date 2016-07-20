defaults_file = "../defaults.lua"
model_file    = "onezone.mod"       
data_atomic_file     = "../../data/chianti_atomdata.hdf5" -- H_3lev_plusC_atomdata.hdf5"

transport_nu_grid  = {3.5e14,30.0e15,0.0001,1}  -- frequency grid
transport_steady_iterate         = 0
transport_radiative_equilibrium  = 1
transport_boundary_in_reflect    = 1
transport_boundary_out_reflect   = 1

t_lc = 1e12/3e10
tstep_max_steps  = 1
tstep_time_start = 0.0
tstep_time_stop  = 2*t_lc
tstep_max_dt     = 0.1*t_lc
tstep_min_dt     = 0.1*t_lc


-- inner source emission
particles_max_total  = 1e8
core_n_emit          = 0
core_radius          = 1.0e14
core_temperature     = 1.0e4
core_photon_frequency   =  5e15
core_luminosity         = 7.12115e40
particles_n_initialize  = 5e8

-- output spectrum
spectrum_nu_grid    = transport_nu_grid
output_write_levels = 1

-- opacity information
opacity_use_nlte      = 1
opacity_bound_bound   = 0
opacity_bound_free    = 1
opacity_free_free     = 0

opacity_grey_opacity  = 0
opacity_epsilon       = 1
line_velocity_width   = 1e8
opacity_electron_scattering = 0





