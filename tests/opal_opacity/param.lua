-- Parameter file to calculatet mean opacities for a range of
-- and temperatures and densities (really logR)

sedona_home        = os.getenv("SEDONA_HOME")
defaults_file      = sedona_home.."/defaults/snopac_defaults.lua"
data_atomic_file   = sedona_home.."/data/cmfgen_newdata.hdf5"
--data_atomic_file   = sedona_home.."/data/cmfgen_levelcap100.hdf5"

-- set this to somefilename to output mean opacities
output_mean_opacities = "mean_opacities.dat"

-- list of isotopes to use: format = Z.A
-- with Z = atomic number, A = atomic weight
elements_Z  = {1,2}
elements_A  = {1,4}
-- mass fractions of each isotope above
mass_fractions = {0.9,0.1}

-- frequency grid to do calculations on (in Hz)
nu1 = 3e10/(1e4*1e-8)
nu2 = 3e10/(0.0001*1e-8)
-- logarthimic frequency spacing by 0.005
transport_nu_grid     = {nu1,nu2,0.001,1}

-- range of temperatures to calculate
-- format is {start, stop, delta}
temperature =  {5e3,5e8,0.1,1}

-- setting use_logR = 1 means that 
-- the "density" array really sets logR
-- where R = rho*(T/1e6)**(-3)
use_logR = 1

-- range of logR to calculate
density     =  {-4.5,-4.5,1}

-- opacities to use
opacity_electron_scattering = 1
opacity_bound_free          = 1
opacity_epsilon             = 1.0
opacity_free_free           = 1
opacity_line_expansion      = 0
opacity_bound_bound         = 0
opacity_fuzz_expansion      = 0
