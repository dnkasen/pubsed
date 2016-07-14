defaults_file = "../defaults.lua"
grid_type  = "grid_1D_sphere"    -- grid geometry; must match input model
model_file   = "constant.mod"        -- input model file



-- time stepping
days = 3600.0*24
tstep_max_steps    = 200
tstep_time_stop    = 2e-10
tstep_max_dt       = 2e-12
tstep_min_dt       = 2e-12
tstel_max_delta    = 1.0
grid_write_out     = 1e-11


transport_nu_grid  = {1,1,1}

-- inner source
particles_n_initialize  = 200000
core_n_emit        =   0
core_radius        =   0
core_luminosity    =   0

-- output spectrum
spectrum_name  = "out";
spectrum_time_grid = {1e-10,1e-8,1e-10}
spectrum_nu_grid   = nu_grid
spectrum_n_mu      = 1
spectrum_n_phi     = 1



-- opacity infor
opacity_grey_opacity             = 1
opacity_epsilon                  = 0

transport_steady_iterate         = 0
transport_radiative_equilibrium  = 0



