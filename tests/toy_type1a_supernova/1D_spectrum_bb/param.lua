
-- model type and file
grid_type    = "grid_1D_sphere"
model_file   = "../models/toy_SNIa_1D_d20.mod" -- toy_SNIa_1D.mod"
hydro_module = "homologous"

-- defaults and atomic data files
sedona_home        = os.getenv('SEDONA_HOME')
defaults_file      = sedona_home.."/defaults/sedona_defaults.lua"
data_atomic_file   = sedona_home.."/data/cmfgen_levelcap100.hdf5"

-- transport properites
transport_nu_grid  = {0.8e14,1.0e16,0.0003,1}  -- frequency grid
transport_radiative_equilibrium  = 1
transport_steady_iterate         = 1

-- output spectrum frequency grid
spectrum_nu_grid   = {0.8e14,1.0e16,0.002,1}

-- time of spectrum calculation
tstep_time_start = 20*3600.0*24.0

-- radioactive particle emission
particles_n_emit_radioactive = 1e5

-- opacity information
opacity_grey_opacity         = 0
opacity_electron_scattering  = 1
opacity_fuzz_expansion       = 0
opacity_line_expansion       = 0
opacity_bound_free           = 1
opacity_bound_bound          = 1
opacity_free_free            = 1
opacity_epsilon              = 1.0
line_velocity_width          = 0.002*3e10

-- output files
output_write_radiation = 1
