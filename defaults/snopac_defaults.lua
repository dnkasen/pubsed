-----------------------------------------------------
-- This file sets the defaults for all parameters
-- used in the snopac opacity module
-----------------------------------------------------

-- output file
write_planck_mean   = 0
-- set this to some filename to output mesa file
output_mesa_file    = "" 
-- set this to somefilename output opacity as a function of frequency (Hz)
output_frequency_file = "" 
-- set this to somefilename output opacity as a function of wavelength (Angs)
output_wavelength_file = ""
-- set this to somefilename output gas state (density, ionization, etc...)
output_gas_state   = ""
-- set this to somefilename to output mean opacities
output_mean_opacities = ""

-- fuzz line data file
sedona_home = os.getenv("SEDONA_HOME")
data_fuzzline_file   = ""
data_atomic_file     = sedona_home.."/data/cmfgen_atomdata.hdf5"

-- if this is set = 1, density array will really set logR
use_logR = 0
-- properties of the gas
density     =  {-13}
temperature =  {5e3}
time        =   3600.0*24.0

-- transport frequency grid
transport_nu_grid  = {1,1,1}

-- list of isotopes to use: format = Z.A
-- with Z = atomic number, A = atomic weight
elements_Z  = {1}
elements_A  = {1}
-- mass fractions of each isotope above
mass_fractions = {1.0}


-- opacities to use
line_velocity_width         = 0
opacity_electron_scattering = 0
opacity_bound_free          = 0
opacity_epsilon             = 1.0
opacity_free_free           = 0
opacity_line_expansion      = 0
opacity_bound_bound         = 0
opacity_fuzz_expansion      = 0