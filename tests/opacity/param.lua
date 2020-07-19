sedona_home   = os.getenv('SEDONA_HOME')
defaults_file    = sedona_home.."/defaults/sedona_defaults.lua"
data_atomic_file = sedona_home.."/data/atom_data_weizmann_SNIA_complete_lc1000.hdf5"
--data_fuzzline_file = sedona_home.."/data/kurucz_cd23_fuzz.hdf5"

model_file    = "onezone.mod"

transport_nu_grid  = {1.e14,3.e15,0.003,1}  -- frequency grid
transport_steady_iterate         = 1
transport_radiative_equilibrium  = 0

-- opacity information
opacity_grey_opacity  = 0
opacity_use_nlte      = 0
opacity_bound_bound   = 0
opacity_bound_free    = 0
opacity_free_free     = 0
opacity_line_expansion = 1
opacity_fuzz_expansion = 0
opacity_electron_scattering = 0

output_write_radiation = 1
