
-- model type and file
grid_type    = "grid_2D_cyln"
model_file   = "../models/toy_SNIa_2D.h5"
hydro_module = "homologous"

-- defaults and atomic data files
sedona_home        = os.getenv('SEDONA_HOME')
defaults_file      = sedona_home.."/defaults/sedona_defaults.lua"
data_atomic_file   = sedona_home.."/data/cmfgen_levelcap100.hdf5"

-- transport properites
transport_nu_grid  = {0.8e14,1.0e16,0.001,1}  -- frequency grid
transport_radiative_equilibrium  = 1
transport_steady_iterate         = 5

-- output spectrum frequency grid
spectrum_nu_grid   = {0.8e14,1.0e16,0.002,1}
spectrum_n_mu = 10

-- time of spectrum calculation
tstep_time_start = 20*3600.0*24.0

-- radioactive particle emission
particles_n_emit_radioactive = 1e6
particles_last_iter_pump     = 10

-- opacity information
opacity_grey_opacity         = 0
opacity_electron_scattering  = 1
opacity_fuzz_expansion       = 0
opacity_line_expansion       = 1
opacity_bound_bound          = 0
opacity_epsilon              = 1

-- output files
output_write_radiation = 0

limits_temp_min = 1e3
limits_temp_max = 1e6