sedona_home   = os.getenv('SEDONA_HOME')

defaults_file    = sedona_home.."/defaults/sedona_defaults.lua"
data_atomic_file = sedona_home.."/data/ASD_atomdata.hdf5"

grid_type  = "grid_1D_sphere"    -- grid geometry; must match input model
model_file = "sod_model.mod"        -- input model file

hydro_module     = "1D_lagrangian"
hydro_gamma_index = 1.4
hydro_cfl         = 0.1
hydro_viscosity_parameter = 3
hydro_mean_particle_mass = 1.0

-- use no transport
transport_module = ""

-- time stepping
tstep_max_steps    = 100000
tstep_time_stop    = 1
tstep_max_dt       = 1

output_write_plt_file_time  = 0.1
