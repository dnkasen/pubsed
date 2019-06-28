-- model file to read in
init_file = "./Radshock.mod"

-- time stepping for radiation module (all in seconds)
tstep_min   =    1.e-15           -- smallest allowed time step
tstep_max   =    10.
tstep_del   =    2.0           -- largest allowed step fraction of current time
n_times     =    1.e8           -- maximum nuber of steps taken
t_stop      =    1.4e4          -- time to stop calculation. 
write_out   =    5.e2           -- maximum time spacing to write out grid files
output_interval = 1.e15          -- maximum number of timesteps between outputs
log_time    =    0

-- simulation options

use_transport = 1
use_hydro     = 1
use_gravity   = 0
radiative_eq  = 0
step_size     = 0.3

-- inner spherical luminosity "core" source
L_core   = 0.              -- luminosity of core (ergs/sec). Eddington luminosity 4 pi G M mp c /sigmaT. 
r_core   = 0.              -- radius of core (cm). If particles are emitted from core, this must be > inner boundary or program will crash!
T_core   = 1e5               -- blackbody temperature of core
 -- number of particles to emit from core each time step
inject_photons = 0 

-- particle creation parameters (for radioactive/thermal emission)
init_photons   =  10      -- number of particles per zone to initialize calculation
emit_min       =  1      -- max # of particles to emit per zone per time step
emit_max       =  1000     -- min # of particles to emit per zone per time step
E_particle     = 1e42*1e3
T_init   = 1e1

-- wavelength grid for opacities; format = {start, stop, delta} in Angstroms
wave_grid    = {0.,1.,2.}

-- wavelength grid for output spectrum in Angstroms
spec_grid    = {0,1.e5,1.e4}
-- time grid for output spectrum
spec_time    = {0,1,2}
-- number of theta and phi bins for output spectrum
n_mu   = 1
n_phi  = 1

-- opacity parameters
grey_opacity = 0.4 -- set to zero for realistic opacity calculations. This is essentially 1 / mfp
epsilon      = 1.     -- ratio of absorptive to total (absorb + scatter) opacity

tau_diffuse  = 0     -- parameter for using discrete diffusion
fleck_alpha  = 1.0     -- parameter for implicit monte carlo

 -- units that we are using, "Angs" (default) or "MeV"
wave_units   = ""  

-- hydro parameters

cfl         = 0.5
gamfac      = 5./3.
c_s         = 0. -- only matters if gamfac is identically 1.
M_gravity   = 0.   -- used for lagrangian grid only
v_piston    = 2.e6   -- used for lagrangian grid only
p_ext       = 0.   -- only matters if set to nonzero value and using lagrangian grid 
C_q         = 0.5 -- for artificial visocsity
L_q         = 0.   -- linear term for artificial visocsity
mu          = 0.5
nghost      = 0

-- boundary conditons

boundary_xL = "reflect"
boundary_xR = "dirichlet"
boundary_yL = "dirichlet"
boundary_yR = "dirichlet"
boundary_zL = "dirichlet"
boundary_zR = "dirichlet"

transport_xL = "reflect"
transport_xR = "outflow"
transport_yL = "outflow"
transport_yR = "outflow"
transport_zL = "outflow"
transport_zR = "outflow"



