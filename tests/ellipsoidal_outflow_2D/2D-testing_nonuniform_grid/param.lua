sedona_home   = os.getenv('SEDONA_HOME')

defaults_file    = sedona_home.."/defaults/sedona_defaults.lua"
data_atomic_file = sedona_home.."/data/ASD_atomdata.hdf5"

-- model file
model_file   = "../models/ellipsoidal_outflow_2D-testing_nonuniform_grid.h5"

-- grid parameters
grid_type    = "grid_2D_cyln"

-- hydro parameters
hydro_module = "homologous"

-- helper variable
hours = 3600.0
days = 3600*24.0

-- transport parameters
transport_module = "monte_carlo"
transport_radiative_equilibrium = 1
transport_use_ddmc = 0
-- transport_ddmc_tau_threshold = 3

-- particle parameters
particles_n_initialize       = 0   -- total number of particles used initially
particles_n_emit_radioactive = 5e4 -- number of particles emitted per time step from radioactivity
particles_n_emit_thermal     = 0   -- number of particles emitted per time step from thermal emission

-- time stepping
tstep_max_steps  = 1e6
tstep_time_start = 0.1*days
tstep_time_stop  = 10.1*days
-- tstep_min_dt     = 0.1*days
tstep_max_dt     = 0.1*days
tstep_max_delta  = 0.1

-- output parameters
output_write_radiation     = 1
output_write_plt_file_time = 10.0*days

-- opacity information
opacity_grey_opacity = 1.e-1

-- output spectrum parameters
spectrum_name	   = "optical_spectrum"
spectrum_time_grid = {0.1*days,10.0*days,0.1*days}
spectrum_n_mu      = 10