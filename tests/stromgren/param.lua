defaults_file = "../defaults.lua"
model_file    = "constant.mod"       
--data_atomic_file     = "../../data/H_3lev_plusC_atomdata.hdf5" 
data_atomic_file     = "../../data/cmfgen_atomdata.hdf5" 


transport_radiative_equilibrium  = 0
transport_steady_iterate         = 30

-- inner source emission
h       = 6.6260755e-27      -- planck's constant (ergs-s)
rydberg = 2.1798741e-11      --  rydberg constant in ergs
nu_t    = rydberg/h          -- threshold frequency
nu      = 1.1*nu_t           -- emit just above threshold
Q   = 5e48                   -- photons/s emitted

core_luminosity       = Q*h*nu
core_spectrum_file    = "inner_spectrum.dat"
core_temperature      = 5e5
core_radius           = 0
core_n_emit           = 0.5e5

-- frequency grid
transport_nu_grid  = {0.1*nu_t,5*nu_t,0.0001,1}  

-- output spectrum
spectrum_nu_grid   = transport_nu_grid

output_write_levels = 1

-- opacity information
opacity_grey_opacity = 0
opacity_free_free    = 0
opacity_bound_bound  = 1
opacity_bound_free   = 1
opacity_use_nlte     = 1
line_velocity_width  = 5e7






