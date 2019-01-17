sedona_home   = os.getenv('SEDONA_HOME')

defaults_file    = sedona_home.."/defaults/sedona_defaults.lua"
data_atomic_file = sedona_home.."/data/H_2lev_noC_atomdata.hdf5"
model_file           = "vacuum.mod"       


transport_nu_grid  = {2e14,4e15,0.0001,1}  -- frequency grid
transport_radiative_equilibrium  = 0
transport_steady_iterate         = 2

-- inner source emission
core_n_emit      = 1e6
core_radius      = 3.0e13
core_temperature = 5.0e4
core_luminosity  = 0 -- set automatically from blackbody T,R

-- output spectrum
spectrum_nu_grid   = transport_nu_grid
output_write_radiation            = 1
output_write_atomic_levels        = 1

-- opacity information
opacity_use_nlte    = 1
opacity_atoms_in_nlte       = {1}
opacity_bound_bound = 1
line_velocity_width = 5e7





