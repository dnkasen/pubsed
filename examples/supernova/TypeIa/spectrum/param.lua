grid_type    = "grid_1D_sphere"   
model_file   = "w7.mod" 
hydro_module = "homologous"


sedona_home        = os.getenv('SEDONA_HOME')
defaults_file      = sedona_home.."/defaults/sedona_defaults.lua"
data_atomic_file   = sedona_home.."/data/ASD_atomdata.hdf5"
data_fuzzline_file = sedona_home.."/data/kurucz_cd23_fuzz.hdf5"

-- transport properites
transport_nu_grid  = {0.8e14,1.0e16,0.001,1}  -- frequency grid
transport_radiative_equilibrium  = 1
transport_steady_iterate         = 6
spectrum_nu_grid   = {0.8e14,1.0e16,0.005,1}


-- inner source emission
texp             = 20*3600.0*24.0
core_n_emit      = 1e5
core_radius      = 4e8*texp
core_luminosity  = 1.0e43
particles_n_emit_radioactive = 1e4
particles_last_iter_pump     = 10

tstep_time_start = texp

-- output spectrum

-- opacity information
opacity_grey_opacity         = 0
opacity_electron_scattering  = 1
opacity_fuzz_expansion       = 1
opacity_epsilon              = 1

-- output files
output_write_radiation = 0
