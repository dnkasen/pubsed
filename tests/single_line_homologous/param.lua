defaults_file = "../defaults.lua"
grid_type  = "grid_1D_sphere"    -- grid geometry; must match input model
model_file = "constant_standard.mod"        -- input model file

-- helper variables
days = 3600.0*24

-- inner source
core_n_emit      = 1e5
core_radius      = 1e9*10*days
core_luminosity  = 1e42
core_temperature = 0.0e4

particles_step_size = 0.2

nu2 = 3e10/(1050*1e-8)
nu1 = 3e10/(1400*1e-8)
dnu = 0.0001 --(nu2-nu1)/500.0

transport_nu_grid   = {nu1,nu2,dnu,1}
spectrum_nu_grid   = {nu1, nu2, dnu,1}

opacity_epsilon             = 0.0
opacity_electron_scattering = 0
opacity_line_expansion      = 0
opacity_bound_bound         = 1
line_velocity_width         = 5e7

transport_steady_iterate = 1
transport_radiative_equilibrium = 1



