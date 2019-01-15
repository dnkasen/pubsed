sedona_home   = os.getenv('SEDONA_HOME')
defaults_file    = sedona_home.."/defaults/sedona_defaults.lua"
data_atomic_file = sedona_home.."/data/2level_atomdata.hdf5"

grid_type  = "grid_1D_sphere"
model_file   = "../models/constant.mod"
hydro_module = "homologous"

-- time stepping
days = 3600.0*24
tstep_max_steps    = 10000
tstep_time_stop    = 80.0*days
tstep_max_dt       = 0.2*days
tstep_min_dt       = 0.2*days
tstep_max_delta    = 0.2

-- particles and opacity
particles_n_initialize  = 1e5
opacity_grey_opacity             = 50000
transport_radiative_equilibrium  = 1
output_write_plt_file_time       = 1*days

transport_use_ddmc = 2
transport_ddmc_tau_threshold = 0.1
