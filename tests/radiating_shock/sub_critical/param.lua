sedona_home   = os.getenv('SEDONA_HOME')

defaults_file    = sedona_home.."/defaults/sedona_defaults.lua"
data_atomic_file = sedona_home.."/data/ASD_atomdata.hdf5"

grid_type  = "grid_1D_sphere"    -- grid geometry; must match input model
model_file = "../models/ensman_model.mod"        -- input model file

hydro_module     = "1D_lagrangian"
hydro_gamma_index = 5.0/3.0
hydro_v_piston    = 6e5
transport_module = "monte_carlo"

particles_n_initialize    = 1000
particles_n_emit_thermal  = 10000
opacity_grey_opacity      = 0.4
opacity_epsilon           = 1

-- time stepping
tstep_max_steps    = 200000
tstep_time_stop    = 0.9e5

output_write_plt_file_time        = 1e4




