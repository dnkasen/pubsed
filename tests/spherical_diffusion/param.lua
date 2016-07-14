defaults_file = "../defaults.lua"
model_file    = "powerlaw.mod"       
data_atomic_file     = "../../data/ground_only_atomdata.hdf5"


transport_radiative_equilibrium  = 1
transport_steady_iterate         = 0

-- time stepping

days = 3600.0*24.0
tstep_max_steps  = 1000
tstep_time_stop  = 200*days
tstep_max_dt     = 1*days
tstep_min_dt     = 1*days
grid_write_out   = 10*days

transport_nu_grid  = {1e13,1e17,0.01,1}  -- frequency grid

core_luminosity       = 1e45
core_temperature      = 1e5
core_radius           = 1e14
core_n_emit           = 2e4

-- output spectrum
spectrum_nu_grid   = {1,1,1}
spectrum_time_grid = {0,110*days,2*days}

-- opacity information
opacity_grey_opacity = 0.2
opacity_epsilon      = 1.0
opacity_free_free    = 0
opacity_bound_free   = 0
opacity_use_nlte     = 0






