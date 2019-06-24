sedona_home   = os.getenv('SEDONA_HOME')
defaults_file    = sedona_home.."/defaults/sedona_defaults.lua"
data_atomic_file = sedona_home.."/data/atom_data.hdf5"
model_file    = "onezone.h5"       

transport_nu_grid  = {3.5e13,3.5e13,0.1,1}  -- frequency grid
transport_steady_iterate         = 0
--transport_radiative_equilibrium  = 0
skip_gas_temp_update_during_transport   = 1
transport_boundary_in_reflect    = 1
transport_boundary_out_reflect   = 1
transport_fleck_alpha            = 0.5

tstep_max_steps    = 1e5
tstep_time_stop    = 2.e-7
tstep_max_dt       = 1.e-11
tstep_min_dt       = 1.e-11
tstel_max_delta    = 1.0

hydro_module = "1D_lagrangian"
hydro_use_transport = 1
hydro_boundary_outflow = 0
hydro_boundary_rigid_outer_wall = 1

particles_n_emit_thermal  = 0
particles_n_initialize    = 1e4

-- output spectrum
spectrum_nu_grid    = transport_nu_grid
output_write_atomic_levels = 0
output_write_radiation = 1
output_write_plt_file_time = 1.e-11
output_write_plt_log_space = 0.5

-- opacity information
opacity_grey_opacity  = 0.4

opacity_use_nlte      = 0
opacity_bound_bound   = 0
opacity_bound_free    = 0
opacity_free_free     = 0
opacity_electron_scattering = 0









