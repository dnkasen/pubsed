
grid_type    = "grid_1D_sphere"   
model_file   = "w7.mod"
hydro_module = "homologous"

sedona_home        = os.getenv('SEDONA_HOME')
defaults_file      = sedona_home.."/defaults/sedona_defaults.lua"
data_atomic_file   = sedona_home.."/data/ASD_atomdata.hdf5"
data_fuzzline_file = sedona_home.."/data/kurucz_cd23_fuzz.hdf5"

-- helper variable
days = 3600.0*24

particles_n_emit_radioactive = 1e4
tstep_time_start             = 2*days
tstep_time_stop              = 50*days
tstep_max_dt                 = 1.0*days
tstep_max_delta              = 0.1

-- grid to calculate and store opacities
nu1 = 1e14
nu2 = 2e16
transport_nu_grid   = {nu1,nu2,0.0003,1}

-- grid to calculate output spectrum
nu1_spec = nu1*1.1
spectrum_nu_grid   = {nu1_spec,nu2,0.005,1}
spectrum_time_grid = {-5*days,100*days,1.0*days}
output_time_write  = 1*days
output_write_radiation = 1

-- opacity settings
opacity_grey_opacity         = 0.0
opacity_epsilon              = 1.0
opacity_electron_scattering  = 1
opacity_free_free            = 0
opacity_bound_bound          = 0
opacity_line_expansion       = 0
opacity_fuzz_expansion       = 1

-- transport settings
transport_steady_iterate        = 0
transport_radiative_equilibrium = 1
