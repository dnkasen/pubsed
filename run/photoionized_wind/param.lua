defaults_file = "./defaults.lua"
model_file    = "wind.mod"       
data_atomic_file     = "../../data/chianti_atomdata.hdf5"


transport_radiative_equilibrium  = 0
transport_steady_iterate         = 4

output_write_levels = 1

transport_nu_grid  = {1e14,8e15,0.00002,1}  -- frequency grid

core_luminosity       = 1e37
core_temperature      = 40000.0
core_radius           = 1e15
core_n_emit           = 2e6
core_spectrum_file    = ""


-- output spectrum
spectrum_nu_grid   = {1e14,5e15,0.00005,1}  -- frequency grid

-- opacity information
opacity_electron_scattering = 0
opacity_bound_bound  = 1
opacity_free_free    = 1
opacity_bound_free   = 1
opacity_use_nlte     = 1
line_velocity_width  = 0






