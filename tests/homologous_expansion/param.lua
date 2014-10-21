grid_type  = "grid_1D_sphere"    -- grid geometry; must match input model
model_file   = "constant.mod"        -- input model file
hydro_module = "homologous"


--
max_total_particles = 1e6
step_size = 0.3

-- time stepping
days = 3600.0*24
max_tsteps = 100
t_stop     = 80.0*days
dt_max     = 0.5*days
dt_min     = 0.0
dt_del     = 1.0


nu_grid  = {0.2e14,5.0e15,0.2e14}

-- inner source
init_particles = 10
n_emit_radioactive = 0;
n_emit_core    =  0
r_core         = 0

-- output spectrum
spectrum_name  = "out";
spec_time_grid = {-10*days,10*days,100*days}
spec_nu_grid   = nu_grid
gamma_nu_grid  = {1,1,1}
spec_n_mu      = 1
spec_n_phi     = 1


-- opacity infor
grey_opacity = 5

steady_iterate = 0
radiative_eq   = 1
write_out      = 1



