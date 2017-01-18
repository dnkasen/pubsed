defaults_file = "defaults.lua"

-- atomic data files (both atomic levels and lines)
data_atomic_file     = "../data/cmfgen_atomdata.hdf5"

-- another file with different lines, 
data_fuzzline_file   = ""

-- output file
write_planck_mean   = 1
-- set this to some filename to output mesa file
output_mesa_file    = "" 
-- set this to somefilename output opacity as a function of frequency (Hz)
output_frequency_file = "" --line_opac.dat"
-- set this to somefilename output opacity as a function of wavelength (Angs)
output_wavelength_file = "line_opac.dat"

-- list of isotopes to use: format = Z.A
-- with Z = atomic number, A = atomic weight
elements       = {26.56, 27.56, 28.56} --1.1, 2.4, 6.12,     8.16,     26.56}
-- mass fractions of each isotope above
mass_fractions = {0.1, 0.8, 0.1} -- 0.7, 0.3, 2.78E-03, 7.56E-03,  1.3e-3  }

-----------------------------------------------------------
-- array format is: (x_start, x_stop, x_delta, log?)
-- log = 0 does uniform spacing:     x_{i+1} = x_i + x_delta
-- log = 1 does logarithmic spacing: x_{i+1} = x_i*(1 + x_delta)
-- only 3 elements in array assumes log = 0
-- only 1 element in array uses only 1 value, that one
-----------------------------------------------------------

-- frequency grid to do calculations on (in Hz)
transport_nu_grid  = {1e14,1e18,0.01,1} 

-- array of log10 density in g/cm^3 (or logR array, see below)
density =  {-13}
-- if this is set = 1, density array will really set logR
use_logR = 1
-- array of temperature in K
temperature = {1e4}
-- time in seconds
time     = 10*60.0*60.0*24.0 -- 10 days

-- opacities to use
opacity_electron_scattering = 0
opacity_bound_free          = 0
opacity_free_free           = 0
opacity_line_expansion      = 1
opacity_bound_bound         = 0

-- another way to include lines
opacity_fuzz_expansion      = 0