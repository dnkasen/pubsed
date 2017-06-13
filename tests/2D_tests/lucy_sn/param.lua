defaults_file = "../../defaults.lua"

grid_type    = "grid_2D_cyln"        -- grid geometry; match input model
model_file   = "lucy_2d.h5"    -- input model file
hydro_module = "homologous"
data_atomic_file     = "../../../data/ASD_atomdata.hdf5"

-- time stepping
days = 3600.0*24
tstep_max_steps  = 1000
tstep_time_stop  = 70.0*days
tstep_max_dt     = 0.5*days
tstep_min_dt     = 0.0
tstep_max_delta  = 0.05

output_write_times  = 50*days  -- how often to write out grid data

-- emission parameters
particles_n_emit_radioactive = 1e3

nu1 = 0.5e14
nu2 = 1.5e15
transport_nu_grid   = {nu1,nu2,0.001,1}
--spectrum_nu_grid   =  {nu1,nu2,0.1,1}

-- output spectrum
spectrum_time_grid = {-5*days,100*days,0.5*days}
spectrum_name = "optical_spectrum"
gamma_name    = "gamma_spectrum"
spectrum_n_mu      = 10

-- opacity parameters
opacity_grey_opacity = 0.1

transport_steady_iterate = 0
transport_radiative_eq   = 1


