-- Parameter file to plot mean opacity for a range of 
-- and temperatures and densities

sedona_home        = os.getenv("SEDONA_HOME")
defaults_file      = sedona_home.."/defaults/snopac_defaults.lua"
data_atomic_file   = sedona_home.."/data/cmfgen_atomdata.hdf5"

-- set this to somefilename to output mean opacities
output_mean_opacities = "mean_opacities.dat"

-- list of isotopes to use: format = Z.A
-- with Z = atomic number, A = atomic weight
elements_Z  = {26} 
elements_A  = {56}
-- mass fractions of each isotope above
mass_fractions = {1.0}

-- frequency grid to do calculations on (in Hz)
nu1 = 3e10/(1e4*1e-8) 
nu2 = 3e10/(990*1e-8) 
-- logarthimic frequency spacing by 0.01
transport_nu_grid     = {nu1,nu2,0.01,1}

-- range of (log) densities and temperatures
-- format is {start, stop, delta} 
density     =  {-15,-13,1}
temperature =  {4e3,2e4,1e3}

days = 3600.0*24.0
time = 20*days

-- opacities to use
opacity_electron_scattering = 0
opacity_bound_free          = 0
opacity_epsilon             = 1.0
opacity_free_free           = 0
opacity_line_expansion      = 1
opacity_bound_bound         = 0
opacity_fuzz_expansion      = 0