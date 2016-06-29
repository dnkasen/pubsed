defaults_file = "../defaults.lua"
model_file    = "vacuum.mod"       
data_atomic_file     = "../../data/H_2lev_noC_atomdata.hdf5"

transport_nu_grid  = {0.25e15,1.0e16,1e13}  -- frequency grid
transport_radiative_equilibrium  = 1
transport_steady_iterate         = 1

-- inner source emission
pi = 3.14159
sb  = 5.6704e-5 
core_n_emit      = 5e5
core_radius      = 3.0e13
core_temperature = 5.0e4
core_luminosity  = 4.0*pi*sb*core_radius^2*core_temperature^4

-- output spectrum
spectrum_nu_grid   = transport_nu_grid

-- opacity information
opacity_use_nlte = 1






