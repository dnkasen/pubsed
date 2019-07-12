sedona_home   = os.getenv('SEDONA_HOME')
defaults_file = "../../defaults/sedona_defaults.lua"
model_file    = "constant.mod"       
--data_atomic_file = "atom_data.hdf5"
data_atomic_file = sedona_home.."/data/H_3lev_plusC_atomdata.hdf5"

transport_radiative_equilibrium  = 1
transport_steady_iterate         = 50

-- inner source emission
h       = 6.6260755e-27      -- planck's constant (ergs-s)
nu      = 3.32e15            -- emission frequency
Q       = 5e48               -- photons/s emitted

core_luminosity       = Q*h*nu
core_spectrum_file    = "source_spectrum.dat"
core_radius           = 0
core_n_emit           = 1e5

-- frequency grid

nu1 = 3.5e11
nu2 = 3.e19
transport_nu_grid  = {nu1,nu2,0.003,1}  

-- output spectrum
spectrum_nu_grid   = {nu1,nu2,0.003,1}  

output_write_radiation = 1
output_write_atomic_levels = 1

-- opacity information
opacity_grey_opacity  = 0

opacity_use_nlte     = 1
opacity_atoms_in_nlte = {1,2,6,7,8}
opacity_free_free    = 1
opacity_bound_bound  = 1
opacity_bound_free   = 1
opacity_electron_scattering = 0
line_velocity_width  = 3e8

limits_temp_max = 1e6
limits_temp_min = 100





