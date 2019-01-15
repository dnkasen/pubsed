sedona_home   = os.getenv('SEDONA_HOME')

defaults_file    = sedona_home.."/defaults/sedona_defaults.lua"
data_atomic_file = sedona_home.."/data/2level_atomdata.hdf5"

grid_type    = "grid_1D_sphere"    -- grid geometry; must match input model
model_file   = "../models/constant.mod"      -- input model file

-- time stepping
tstep_max_steps    = 1000
tstep_time_stop    = 1e-8
tstep_max_dt       = 5e-10
tstep_min_dt       = 0.0
tstel_max_delta    = 1.0

-- inner source
particles_n_initialize       = 1e5
transport_use_ddmc           = 0
transport_ddmc_tau_threshold = 0.1

-- opacity
opacity_grey_opacity             = 1.0
opacity_epsilon                  = 0.0

output_write_plt_file_time        = 5e-10
spectrum_time_grid = {0, 1e-9, 1e-7}
