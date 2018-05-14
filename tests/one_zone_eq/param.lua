defaults_file = "../../defaults/sedona_defaults.lua"
model_file    = "onezone.mod"       
data_atomic_file     = "../../data/chianti_atomdata.hdf5" -- H_3lev_plusC_atomdata.hdf5"

transport_nu_grid  = {3.5e13,3.e17,0.01,1}  -- frequency grid
transport_steady_iterate         = 0
transport_radiative_equilibrium  = 0
transport_boundary_in_reflect    = 1
transport_boundary_out_reflect   = 1

t0 = 1e4
tstep_max_steps  = 6
tstep_time_start = 0.0
tstep_time_stop = 3 * t0
tstep_max_dt     = t0
tstep_min_dt     = t0
output_write_times  = t0

-- inner source emission
particles_max_total  = 1e7
core_n_emit          = 0
core_radius          = 1.0e14
particles_n_emit_thermal  = 1e6
particles_n_initialize  = 0

-- output spectrum
spectrum_nu_grid    = transport_nu_grid
output_write_levels = 1
output_write_radiation = 1

-- opacity information
opacity_grey_opacity  = 0

opacity_use_nlte      = 1
opacity_bound_bound   = 1
opacity_bound_free    = 1
opacity_free_free     = 1
opacity_electron_scattering = 1
line_velocity_width   = 3e8









