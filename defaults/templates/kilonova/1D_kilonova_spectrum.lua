-- Template parameter file to run a snapshot spectrum 
-- of a 1D kilonova model at a certain time 

sedona_home        = os.getenv("SEDONA_HOME")

grid_type    = "grid_1D_sphere"
model_file   = "knova.mod"
hydro_module = "homologous"

defaults_file      = sedona_home.."/defaults/sedona_defaults.lua"
data_atomic_file   = sedona_home.."/data/ASD_atomdata.hdf5"


-- helper variable
days = 3600.0*24

-- transport settings
transport_steady_iterate        = 4
transport_radiative_equilibrium = 1

tstep_time_start               = 2.5*days
particles_n_emit_radioactive   = 1e5
particles_last_iter_pump       = 5
force_rprocess_heating         = 1

-- grid to calculate and store opacities
nu1  = 0.3e14
nu2  = 2e16
transport_nu_grid   = {nu1,nu2,0.005,1}
spectrum_nu_grid    = {nu1,nu2,0.005,1}

spectrum_time_grid = {2*days,2*days,0.05*days}

transport_store_Jnu          = 0
output_write_radiation       = 0

-- opacity settings
opacity_grey_opacity         = 0
opacity_epsilon              = 1.0
opacity_electron_scattering  = 1
opacity_free_free            = 1
opacity_line_expansion       = 1
opacity_bound_free           = 0
opacity_bound_bound          = 0
opacity_fuzz_expansion       = 0
opacity_user_defined         = 0


