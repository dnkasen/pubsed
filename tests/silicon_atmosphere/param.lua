defaults_file = "/Users/kasen/codes/sedona6/tests/defaults.lua"
grid_type  = "grid_1D_sphere"    -- grid geometry; must match input model
model_file = "/Users/kasen/codes/sedona6/tests/silicon_atmosphere/constant.mod"        -- input model fil

-- helper variables
days = 3600.0*24

-- inner source
core_n_emit      = 1e5
core_radius      = 1e9*10*days
core_luminosity  = 1e42
core_temperature = 0.6e4

particles_step_size = 0.1

transport_nu_grid   = {0.5e14,1.0e15,1e12}
spectrum_nu_grid   = transport_nu_grid

opacity_epsilon  = 0.0
opacity_electron_scattering = 0
opacity_line_expansion      = 0
opacity_lines               = 1
line_velocity_width         = 1e8

transport_steady_iterate = 1
transport_radiative_equilibrium = 1



