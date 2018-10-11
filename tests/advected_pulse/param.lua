sedona_home   = os.getenv('SEDONA_HOME')

defaults_file    = sedona_home.."/defaults/sedona_defaults.lua"
data_atomic_file = sedona_home.."/data/2level_atomdata.hdf5"

grid_type    = "grid_1D_sphere"    -- grid geometry; must match input model
model_file   = "constant.mod"      -- input model file

-- time stepping
tstep_max_steps    = 1000
tstep_time_stop    = 2e-10
tstep_max_dt       = 2e-12
tstep_min_dt       = 2e-12
tstel_max_delta    = 1.0

-- inner source
particles_n_initialize  = 1e6

-- opacity
opacity_grey_opacity             = 1.0
opacity_epsilon                  = 0.0

output_write_plt_file_time        = 1e-11
