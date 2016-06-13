defaults_file = "../defaults.lua"
grid_type  = "grid_1D_sphere"    -- grid geometry; must match input model
model_file   = "constant.mod"        -- input model file
hydro_module = "homologous"



-- time stepping
days = 3600.0*24
tstep_max_steps    = 100
tstep_time_stop    = 80.0*days
tstep_max_dt       = 1*days
tstep_min_dt       = 1*days
tstel_max_delta    = 1.0


transport_nu_grid  = {0.2e14,5.0e15,0.2e14}

-- inner source
particles_n_initialize  = 10
core_n_emit        =   0
core_radius        =   0
core_luminosity    = 1e43

-- output spectrum
spectrum_name  = "out";
spectrum_time_grid = {-10*days,10*days,100*days}
spectrum_nu_grid   = nu_grid
spectrum_n_mu      = 1
spectrum_n_phi     = 1



-- opacity infor
opacity_grey_opacity             = 5
transport_steady_iterate         = 0
transport_radiative_equilibrium  = 1
transport_write_out              = 1*days



