grid_type  = "grid_1D_sphere"    -- grid geometry; must match input model
model_file = "cc1.mod"   -- input model file

defaults_file      = "defaults.lua"
data_atomic_file   = "../../data/cmfgen_atomdata.hdf5"

-- inner source "lightbulb"
core_n_emit      = 0
-- radioacxtive emission parameters
particles_n_emit_radioactive = 1e4

particles_step_size = 0.2

-- time stepping
days = 3600.0*24
tstep_max_steps  = 1
tstep_time_start = 500.0*days
tstep_time_stop  = 500.0*days
tstep_max_dt     = 0.5*days
tstep_min_dt     = 0.0
tstep_max_delta  = 0.05

-- grid to calculate and store opacities
-- don't need to resolve wavelength for gamma-ray transport
nu2 = 3e10/(10000*1e-8)
nu1 = 3e10/(10000*1e-8)
dnu = (nu2-nu1)/100
transport_nu_grid   = {nu1,nu2,dnu}

-- grid to calculate output spectrum
spectrum_nu_grid   = transport_nu_grid

-- opacity settings
opacity_grey_opacity         = 0.00000000001
opacity_epsilon              = 1

-- transport settings
transport_steady_iterate        = 1
transport_radiative_equilibrium = 1



