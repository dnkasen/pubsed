
-- get environment variable specifying location of sedona
sedona_home   = os.getenv('SEDONA_HOME')

-- name of file which defines the default parameter values
defaults_file    = sedona_home.."/defaults/sedona_defaults.lua"

-- name of file holding atomic data
data_atomic_file = sedona_home.."/data/2level_atomdata.hdf5"

-- name of model file specifying initial conditions
model_file    = "../models/vacuum_1D.mod"

-- -------------------------------------------------

-- Do steady state run, with one iteration
transport_steady_iterate         = 1

-- Assume radiative heating = cooling to calculate gas temp
transport_radiative_equilibrium  = 1

-- -------------------------------------------------
-- Frequency grid to use for the calculation
-- Format is:
--     {nu_start, nu_stop, nu_delta, do_log}
-- since do_log != 0, this creates a grid from nu_start
-- to nu_stop with logarithmically spacing nu_delta
-- -------------------------------------------------
transport_nu_grid  = {1.0e14,1.0e16,0.01,1}  -- frequency grid

-- Frequency grid of the output spectrum, set to be the same
-- as the transport grid (this is recommended but not required)
spectrum_nu_grid   = transport_nu_grid

-- Radius and temperature of the spherical inner core
-- code will calculate luminosity as L = 4*pi*R^2*sigma*T^4
core_radius      = 6.0e14
core_temperature = 1.0e4

-- number of photon packets to emit from the core
core_n_emit      = 2e5

-- assume use grey opacity with the set value
opacity_grey_opacity = 1e-10

-- Flag to write out the full radiation properties in each zone
output_write_radiation = 1
