grid_type    = "grid_1D_sphere"   
model_file   = "tde_static_1D_atmosphere.mod"
defaults_file      = "/homes/nroth1/sedona6/defaults/sedona_defaults.lua"
data_atomic_file   = "/homes/nroth1/sedona6/data/atom_data.hdf5"

-- inner source emission
core_n_emit          = 2e8
core_radius          = 1.0e14
core_luminosity  = 8.45e46
core_temperature = 3.29e5

particles_n_emit_thermal  = 0
particles_n_initialize  = 1e7
particles_n_emit_radioactive = 0

tstep_max_steps    = 40
tstep_time_stop    = 1.6e7
tstep_max_dt       = 4.e5
tstep_min_dt       = 4.e5
tstel_max_delta    = 1.0

-- grid to calculate and store opacities
nu1 = 2.9979e13
nu2 = 2.9979e17
transport_nu_grid   = {nu1,nu2,0.0025,1}

-- grid to calculate output spectrum
spectrum_nu_grid   = transport_nu_grid
spectrum_time_grid = {0.,1.6e7,4.e5}
output_write_radiation = 1

-- opacity settings
opacity_grey_opacity         = 0
opacity_use_nlte             = 1
opacity_atoms_in_nlte = {1,2,6,7,8}
opacity_electron_scattering  = 1
opacity_free_free            = 1
opacity_bound_free           = 1
opacity_bound_bound          = 1
opacity_line_expansion       = 0
opacity_fuzz_expansion       = 0
line_velocity_width          = 1.e9

-- transport settings
transport_steady_iterate        = 0
transport_radiative_equilibrium = 2
