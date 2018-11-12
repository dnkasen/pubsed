
grid_type    = "grid_1D_sphere"   
model_file   = "data_ejecta_paul_faily_radmetz_t9.0e2_i150_t9.0e2_vmin4.71e-3_vout6.0e-1.txt"
hydro_module = "homologous"

sedona_home        = os.getenv('SEDONA_HOME')
sedona_scratch     = os.getenv('SEDONA_SCRATCH')
defaults_file      = sedona_scratch.."/data/sedona_defaults.lua"
data_atomic_file   = sedona_scratch.."/data/ASD_atomdata.hdf5"
data_fuzzline_file = sedona_scratch.."/data/kurucz_cd23_fuzz.hdf5"

-- helper variable
days = 3600.0*24
hours = 3600.

tstep_max_steps              = 1e8 -- meant to be very large
tstep_time_start             = 900.0
tstep_time_stop              = 0.5*hours
tstep_max_dt                 = 0.25*hours
tstep_max_delta              = 0.1
particles_n_initialize = 100.0
particles_n_emit_radioactive = 0.0
force_rprocess_heating         = 1
dont_decay_composition         = 1

-- grid to calculate and store opacities
nu1 = 1e13
nu2 = 2e16
transport_nu_grid   = {nu1,nu2,0.0003,1}

-- grid to calculate output spectrum
spectrum_name = "spectrum_600h";
nu1_spec = nu1*1.1
spectrum_nu_grid   = {nu1_spec,nu2,0.005,1}
spectrum_time_grid = {-1*hours,600*hours,0.25*hours}
output_time_write  = 1*days
output_write_plt_file_time = 48 * 0.25*hours

-- opacity settings
opacity_grey_opacity         = 0.01
opacity_epsilon              = 1.0
opacity_electron_scattering  = 0 -- 1
opacity_free_free            = 0
opacity_bound_bound          = 0
opacity_line_expansion       = 0
opacity_fuzz_expansion       = 0 -- 1

-- transport settings
transport_steady_iterate        = 0
transport_radiative_equilibrium = 1
