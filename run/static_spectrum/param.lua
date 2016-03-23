grid_type  = "grid_1D_sphere"    -- grid geometry; must match input model
model_file = "input_model.mod"   -- input model file

defaults_file      = "./defaults.lua"
data_atomic_file   = "../../data/cmfgen_atomdata.hdf5"

-- helper variable
days = 3600.0*24

-- inner source "lightbulb"
core_n_emit      = 1e5
core_radius      = 1e9*50*days
core_luminosity  = 1e43
core_temperature = 0.6e4

particles_step_size = 0.2

-- grid to calculate and store opacities
nu2 = 3e10/(3000*1e-8)
nu1 = 3e10/(10000*1e-8)
dnu = (nu2-nu1)/500.0
transport_nu_grid   = {nu1,nu2,dnu}

-- grid to calculate output spectrum
nu1_spec = nu1*1.1
spectrum_nu_grid   = {nu1_spec,nu2,dnu}

-- opacity settings
opacity_epsilon              = 0.0
opacity_electron_scattering  = 1
opacity_line_expansion       = 0
opacity_lines                = 1
line_velocity_width          = 5e7

-- transport settings
transport_steady_iterate        = 3
transport_radiative_equilibrium = 1



