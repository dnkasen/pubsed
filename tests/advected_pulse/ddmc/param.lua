sedona_home   = os.getenv('SEDONA_HOME')

defaults_file    = sedona_home.."/defaults/sedona_defaults.lua"
data_atomic_file = sedona_home.."/data/2level_atomdata.hdf5"

grid_type    = "grid_1D_sphere"    -- grid geometry; must match input model
model_file   = "../models/constant.mod"      -- input model file

-- time stepping
tstep_max_steps    = 1000
tstep_time_stop    = 10e-10
tstep_max_dt       = 2e-12
tstep_min_dt       = 2e-12
tstel_max_delta    = 1.0

-- transport
particles_n_initialize       = 1e6
transport_use_ddmc           = 1
transport_imd_ddmc_switch    = 2
transport_ddmc_tau_threshold = 10

-- opacity
opacity_grey_opacity             = 30.0
opacity_epsilon                  = 0.0

output_write_plt_file_time   = 5e-11
output_write_radiation       = 1
