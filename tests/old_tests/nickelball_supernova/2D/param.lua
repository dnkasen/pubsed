grid_type    = "grid_2D_cyln"    -- grid geometry; must match input model
model_file   = "../models/nickelball_2D.h5"   -- input model file
hydro_module = "homologous"

sedona_home        = os.getenv("SEDONA_HOME")
defaults_file      = sedona_home.."/tests/defaults.lua"
data_atomic_file   = sedona_home.."/data/ASD_atomdata.hdf5"
data_fuzzline_file = sedona_home.."/data/kurucz_cd23_fuzz.hdf5"

-- helper variable
days = 3600.0*24

particles_n_emit_radioactive = 1e4
tstep_time_start = 1*days
tstep_time_stop  = 50.*days
tstep_max_dt     = 1.0*days
tstep_max_delta  = 0.02

-- grid to calculate and store opacities
nu1 = 1e14
nu2 = 2e16
transport_nu_grid   = {nu1,nu2,0.003,1}

-- grid to calculate output spectrum
nu1_spec = nu1*1.1
spectrum_nu_grid   = {nu1_spec,nu2,0.01,1}
spectrum_time_grid = {-5*days,100*days,0.5*days}
spectrum_n_mu      = 10
output_write_times = 2*days
output_write_grid   = 0

-- opacity settings
opacity_grey_opacity         = 0.0
opacity_epsilon              = 1.0
opacity_electron_scattering  = 1
opacity_free_free            = 1
opacity_bound_bound          = 0
opacity_line_expansion       = 0
opacity_fuzz_expansion       = 1

line_velocity_width          = 5e7

-- transport settings
transport_steady_iterate        = 0
transport_radiative_equilibrium = 1



