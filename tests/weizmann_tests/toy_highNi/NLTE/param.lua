grid_type    = "grid_1D_sphere"   
model_file   = "snia_toy_m1e1nip6.mod"
hydro_module = "homologous"

defaults_file      = "/homes/nroth1/sedona6/defaults/sedona_defaults.lua"
data_atomic_file   = "atom_data_levelcap100.hdf5"

-- helper variable
days = 3600.0*24

particles_n_initialize       = 5e5
particles_n_emit_radioactive = 1e5
--particles_n_emit_radioactive_per_zone = 1000
tstep_time_start             = 2*days
tstep_time_stop              = 50*days
tstep_max_dt                 = 1.0*days
tstep_max_delta              = 0.1

-- grid to calculate and store opacities
nu1 = 2.9979e13
nu2 = 2.9979e17
transport_nu_grid   = {nu1,nu2,0.0003,1}

-- grid to calculate output spectrum
nu1_spec = nu1*1.1
spectrum_nu_grid   = {nu1_spec,nu2,0.005,1}
spectrum_time_grid = {-5*days,100*days,1.0*days}
output_time_write  = 1*days
output_write_radiation = 1
output_write_mass_fractions = 1

-- opacity settings
opacity_grey_opacity         = 0.0
opacity_use_nlte             = 1
opacity_electron_scattering  = 1
opacity_free_free            = 1
opacity_bound_free           = 1
opacity_bound_bound          = 1
opacity_line_expansion       = 0
opacity_fuzz_expansion       = 0
line_velocity_width   = 1.e8


-- transport settings
transport_steady_iterate        = 0
transport_radiative_equilibrium = 1
