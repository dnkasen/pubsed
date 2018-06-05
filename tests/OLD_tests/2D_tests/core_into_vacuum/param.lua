defaults_file = "../../defaults.lua"
model_file    = "vacuum_2d.h5"       
data_atomic_file     = "../../../data/2level_atomdata.hdf5"
grid_type  = "grid_2D_cyln"    -- grid geometry; must match input model

transport_nu_grid  = {1e13,5.0e15,0.01,1}  -- frequency grid
transport_radiative_equilibrium  = 1
transport_steady_iterate         = 1

-- inner source emission
core_n_emit      = 1e6
core_radius      = 5.0e14
core_luminosity  = 1.0e43
core_temperature = 1.0e4
-- output spectrum
spectrum_nu_grid   = transport_nu_grid
spectrum_n_mu      = 10

-- opacity information
opacity_grey_opacity = 1e-10
opacity_use_nlte = 0






