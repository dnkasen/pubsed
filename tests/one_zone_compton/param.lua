sedona_home   = os.getenv('SEDONA_HOME')
defaults_file    = sedona_home.."/defaults/sedona_defaults.lua"
data_atomic_file = sedona_home.."/data/atom_data.hdf5"
model_file    = "compton.h5"       

transport_nu_grid  = {2.e13,6.e17,0.01,1}  -- frequency grid
transport_radiative_equilibrium  = 0
transport_steady_iterate         = 0
transport_boundary_in_reflect    = 1
transport_boundary_out_reflect   = 1

t_c =  74361.
tstep_max_steps    = 30
tstep_time_stop    = 30 * t_c
tstep_max_dt       = t_c
tstep_min_dt       = t_c
tstel_max_delta    = 1.0


-- inner source emission
particles_max_total  = 1e4
core_n_emit          = 0
particles_n_emit_thermal  = 0
core_radius          = 0.
particles_n_initialize  = 1e4
particles_init_photon_frequency = 2.08366179e14;

-- output spectrum
spectrum_nu_grid   		    = transport_nu_grid
spectrum_time_grid = {0.,30. * t_c,t_c}
output_write_atomic_levels	= 0
output_write_radiation		= 1
output_write_plt_file_time  = t_c

-- opacity information
opacity_grey_opacity  = 0
opacity_use_nlte      = 0
opacity_bound_bound   = 0
opacity_bound_free    = 0
opacity_free_free     = 0
opacity_electron_scattering = 1
opacity_compton_scatter_photons = 1
