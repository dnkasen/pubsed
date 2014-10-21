grid_type  = "grid_1D_sphere"    -- grid geometry; must match input model
model_file = "vacuum.mod"        -- input model file
hydro_module = "none"

-- helper variables
days = 3600.0*24

--
max_total_particles = 1e6
step_size = 0.3

nu_grid  = {0.2e14,5.0e15,0.2e14}

-- inner source
init_particles = 0
n_emit_core  = 1e5
r_core       = 0.5e15
L_core       = 1e43
T_core       = 1.0e4
n_emit_radioactive = 0;

-- output spectrum
spectrum_name  = "out";
spec_time_grid = {-10*days,10*days,100*days}
spec_nu_grid   = nu_grid
gamma_nu_grid  = {1,1,1}
spec_n_mu      = 1
spec_n_phi     = 1


-- opacity infor
grey_opacity = 0.00001
radiative_eq   = 1
steady_iterate = 1
write_out      = 1
epsilon        = 1

atomic_data = ""
