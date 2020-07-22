grid_type    = "grid_1D_sphere"   
model_file   = "toy_SNIa_1D_d20_500.mod"
hydro_module = "homologous"

sedona_home        = os.getenv('SEDONA_HOME')
defaults_file      = sedona_home.."/defaults/sedona_defaults.lua"
data_atomic_file   = sedona_home.."/data/cmfgen_levelcap100.hdf5"
data_fuzzline_file = sedona_home.."/data/kurucz_cd23_fuzz.hdf5"


-- helper variable
days = 3600.0*24

particles_n_emit_radioactive = 1e4
tstep_time_start             = 20.0*days

-- grid to calculate and store opacities
nu1 = 1e14
nu2 = 2e16
transport_nu_grid   = {nu1,nu2,0.0003,1}

-- grid to calculate output spectrum
nu1_spec = nu1*1.1
spectrum_nu_grid   = {nu1_spec,nu2,0.005,1}
spectrum_time_grid = {20.0*days,20.0*days,1}
output_time_write  = 1*days
output_write_radiation = 1
output_write_atomic_levels = 1

-- opacity settings
opacity_grey_opacity         = 0.0
opacity_use_nlte             = 1
opacity_atoms_in_nlte = { 20 }
opacity_epsilon              = 1.0
opacity_electron_scattering  = 1
opacity_free_free            = 1
opacity_bound_free           = 1
opacity_bound_bound          = 1
opacity_line_expansion       = 0
opacity_fuzz_expansion       = 0
line_velocity_width   = 5.e7

-- transport settings
transport_steady_iterate        = 5
transport_radiative_equilibrium = 1
particles_last_iter_pump        = 10
