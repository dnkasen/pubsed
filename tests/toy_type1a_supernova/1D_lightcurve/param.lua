
-- model type and file
grid_type    = "grid_1D_sphere"
model_file   = "../models/toy_SNIa_1D.mod"
hydro_module = "homologous"

-- defaults and atomic data files
sedona_home        = os.getenv('SEDONA_HOME')
defaults_file      = sedona_home.."/defaults/sedona_defaults.lua"
data_atomic_file   = sedona_home.."/data/cmfgen_levelcap100.hdf5"

-- transport properites
transport_nu_grid  = {0.8e14,1.0e16,0.001,1}  -- frequency grid
transport_radiative_equilibrium  = 1
-- time stepping
days = 3600.0*24
tstep_max_steps  = 1000
tstep_time_stop  = 70.0*days
tstep_max_dt     = 0.5*days
tstep_min_dt     = 0.0
tstep_max_delta  = 0.05

-- output spectrum frequency grid
spectrum_nu_grid   = {0.8e14,1.0e16,0.002,1}
spectrum_time_grid = {-0.5*days,100*days,0.5*days}

-- radioactive particle emission
particles_n_emit_radioactive = 1e5

-- opacity information
opacity_grey_opacity         = 0
opacity_electron_scattering  = 1
opacity_fuzz_expansion       = 0
opacity_line_expansion       = 1
opacity_bound_bound          = 0
opacity_epsilon              = 1

-- output files
output_write_radiation = 0
