grid_type    = "grid_1D_sphere"        -- grid geometry; match input model
model_file   = "models/lucy_1D.mod"    -- input model file
atomic_data  = "/Users/kasen/codes/sedona_base/opacity_module/data/atoms"


hydro_module = "homologous"

-- time stepping
days = 3600.0*24
max_tsteps = 1000
t_stop     = 70.0*days
dt_max     = 0.5*days
dt_min     = 0.0
dt_del     = 0.05

--
max_total_particles = 1e6
step_size = 0.3
nu_grid  = {1e15,1.0e15,1e14}

-- emission parameters
init_particles = 0
n_emit_radioactive = 1e3
n_emit_core  = 0
r_core       = 0

-- output spectrum
spectrum_name  = "output";
spec_time_grid = {-5*days,100*days,1*days}
spec_nu_grid   = nu_grid
spec_n_mu      = 1
spec_n_phi     = 1
gamma_nu_grid  = {1,1,0.01}

-- opacity parameters
grey_opacity = 0.1
epsilon = 1

steady_iterate = 0
radiative_eq   = 1
write_out      = 100

