defaults_file = "../defaults.lua"
model_file    = "vacuum.mod"       

transport_nu_grid  = {0.2e14,5.0e15,0.2e14}  -- frequency grid
transport_radiative_equilibrium  = 1
transport_steady_iterate         = 1

-- inner source emission
core_n_emit      = 1e5
core_radius      = 5.0e14
core_luminosity  = 1.0e43
core_temperature = 1.0e4

-- output spectrum
spectrum_nu_grid   = transport_nu_grid

-- opacity information
opacity_grey_opacity = 1e-20







