defaults_file = "../defaults.lua"
grid_type  = "grid_1D_sphere"    -- grid geometry; must match input model
model_file = "sod_model.mod"        -- input model file

hydro_module     = "1D_lagrangian"
hydro_gamma_index = 1.4
transport_module = ""


-- time stepping
tstep_max_steps    = 10000
tstep_time_stop    = 1
tstep_max_dt       = 1

output_write_levels = 0
output_write_grid   = 0
output_write_times  = 0.1



