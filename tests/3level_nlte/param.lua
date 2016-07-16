defaults_file = "../defaults.lua"
model_file    = "vacuum.mod"       
data_atomic_file     = "../../data/H_3lev_plusC_atomdata.hdf5"

transport_nu_grid  = {0.5e14,2.0e16,0.0005,1}  -- frequency grid
transport_radiative_equilibrium  = 0
transport_steady_iterate         = 1

output_write_levels = 1
-- inner source emission
core_n_emit      = 2e6
core_radius      = 3.0e13
core_temperature = 5.0e4
core_luminosity  = 0

-- output spectrum
spectrum_nu_grid   = transport_nu_grid

-- opacity information
opacity_use_nlte    = 1
opacity_bound_free  = 1
opacity_bound_bound = 1
line_velocity_width = 1e8





