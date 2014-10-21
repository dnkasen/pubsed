grid_type  = "grid_1D_sphere"    -- grid geometry; must match input model
model_file = "constant.mod"        -- input model fil
atomic_data  = "/Users/kasen/codes/sedona_base/opacity_module/data/atoms"

hydro_module = "none"

-- helper variables
days = 3600.0*24

--
max_total_particles = 1e6
step_size = 0.3

nu_grid  = {0.2e14,5.0e15,0.1e14}

-- inner source
init_particles = 0
n_emit_core  = 2e4
r_core       = 0.5e15
L_core       = 1e42
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
grey_opacity   = 0
epsilon        = 1
radiative_eq   = 1
steady_iterate = 3
write_out      = 1



