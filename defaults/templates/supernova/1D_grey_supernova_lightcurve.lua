
-- model type and file
grid_type    = "grid_1D_sphere"
model_file   = "mymodel.h5"
hydro_module = "homologous"

-- defaults and atomic data files
sedona_home        = os.getenv('SEDONA_HOME')
defaults_file      = sedona_home.."/defaults/sedona_defaults.lua"
data_atomic_file = sedona_home.."/data/ASD_atomdata.hdf5"

-- transport properites
transport_radiative_equilibrium  = 1
-- time stepping
days = 3600.0*24
tstep_max_steps  = 1000
tstep_time_stop  = 70.0*days
tstep_max_dt     = 0.5*days
tstep_min_dt     = 0.0
tstep_max_delta  = 0.05

-- output spectrum frequency grid
spectrum_time_grid = {-0.5*days,100*days,0.5*days}

-- radioactive particle emission
particles_n_emit_radioactive = 1e4

-- opacity information
opacity_grey_opacity         = 0.1

-- output files
output_write_radiation = 0
