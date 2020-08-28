
sedona_home   = os.getenv('SEDONA_HOME')

defaults_file    = sedona_home.."/defaults/sedona_defaults.lua"
data_atomic_file = sedona_home.."/data/2level_atomdata.hdf5"

model_file    = "../models/vacuum_1D.mod"       

-- transport properites
transport_nu_grid  = {0.2e14,5.0e15,0.01,1}  -- frequency grid
transport_radiative_equilibrium  = 1
transport_steady_iterate         = 1

-- inner source emission
core_n_emit      = 2e5
core_radius      = 5.0e14
core_luminosity  = 1.0e43
core_temperature = 1.0e4

-- output spectrum
spectrum_nu_grid   = transport_nu_grid

-- opacity information
opacity_grey_opacity = 1e-10

-- output files
output_write_radiation = 1

limits_temp_min = 0.
limits_temp_max = 1e6