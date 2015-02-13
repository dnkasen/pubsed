#include <mpi.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <sstream>
#include <limits>
#include <vector>
#include <cassert>
#include <list>
#include <algorithm>

#include "transport.h"
#include "ParameterReader.h"
#include "VoigtProfile.h"
#include "physical_constants.h"

using std::cout;
namespace pc = physical_constants;


//----------------------------------------------------------------------------
// Initialize the transport module
// Includes setting up the grid, particles,
// and MPI work distribution
//----------------------------------------------------------------------------
void transport::init(ParameterReader* par, grid_general *g)
{ 
  params_  = par;
  grid = g;

  // get mpi rank
  MPI_Comm_size( MPI_COMM_WORLD, &MPI_nprocs );
  MPI_Comm_rank( MPI_COMM_WORLD, &MPI_myID  );
  MPI_real = ( sizeof(real)==4 ? MPI_FLOAT : MPI_DOUBLE );
  verbose = (MPI_myID==0);

  // setup and seed random number generator
  const gsl_rng_type * TypeR;
  gsl_rng_env_setup();
  gsl_rng_default_seed = (unsigned int)time(NULL) + MPI_myID;
  TypeR = gsl_rng_default;
  rangen = gsl_rng_alloc (TypeR);
  
  // read relevant parameters
  step_size_          = params_->getScalar<double>("particles_step_size");
  max_total_particles = params_->getScalar<int>("particles_max_total");
  radiative_eq   = params_->getScalar<int>("transport_radiative_equilibrium");
  steady_state   = (params_->getScalar<int>("transport_steady_iterate") > 0);

  // initialize the frequency grid  
  std::vector<double> nu_dims = params_->getVector<double>("transport_nu_grid");
  if (nu_dims.size() != 3) {
    cout << "# improperly defined nu_grid; need {nu_1, nu_2, dnu}; exiting\n";
    exit(1); }
  nu_grid.init(nu_dims[0],nu_dims[1],nu_dims[2]);
    
  // intialize output spectrum
  std::vector<double>stg = params_->getVector<double>("spectrum_time_grid");
  std::vector<double>sng = params_->getVector<double>("spectrum_nu_grid");
  int nmu  = params_->getScalar<int>("spectrum_n_mu");
  int nphi = params_->getScalar<int>("spectrum_n_phi");
  optical_spectrum.init(stg,sng,nmu,nphi);
  std::vector<double>gng = params_->getVector<double>("gamma_nu_grid");  
  gamma_spectrum.init(stg,sng,nmu,nphi);


  // initialize nlte_gas class
  std::string atomdata = params_->getScalar<string>("data_atomic_file");  
  gas.initialize(atomdata,grid->elems_Z,grid->elems_A,nu_grid);
  std::string fuzzfile = params_->getScalar<string>("data_fuzzline_file");  
  int nl = gas.read_fuzzfile(fuzzfile);
  if (verbose) std::cout << "# From fuzzfile \"" << fuzzfile << "\" " << 
		 nl << " lines used\n";
  
  // set gas opacity flags and parameters
  gas.epsilon_      = params_->getScalar<double>("opacity_epsilon");
  gas.grey_opacity_ = params_->getScalar<double>("opacity_grey_opacity");
  gas.use_electron_scattering_opacity 
    = params_->getScalar<int>("opacity_electron_scattering");
  gas.use_line_expansion_opacity  
    = params_->getScalar<int>("opacity_line_expansion");
  gas.use_fuzz_expansion_opacity  
    = params_->getScalar<int>("opacity_fuzz_expansion");
  gas.use_bound_free_opacity  
    = params_->getScalar<int>("opacity_bound_free");
  gas.use_free_free_opacity  
    = params_->getScalar<int>("opacity_free_free");

  // parameters for treatment of detailed lines
  use_detailed_lines_  = params_->getScalar<double>("opacity_lines");
  line_velocity_width_ = params_->getScalar<double>("line_velocity_width");
  
  // get line frequencies and ion masses
  line_nu_ = gas.get_line_frequency_list();
  line_sqrt_Mion_ = gas.get_line_ion_mass_list();
  for (int i=0;i<line_sqrt_Mion_.size();i++)
    line_sqrt_Mion_[i] = sqrt(line_sqrt_Mion_[i]);

  // allocate and initalize space for opacities
  n_lines_ = gas.get_number_of_lines();
  line_opacity_.resize(grid->n_zones);
  abs_opacity_.resize(grid->n_zones);
  scat_opacity_.resize(grid->n_zones);
  emissivity_.resize(grid->n_zones);

  for (int i=0; i<grid->n_zones;  i++)
  {
    line_opacity_[i].resize(n_lines_);
    abs_opacity_[i].resize(nu_grid.size());
    scat_opacity_[i].resize(nu_grid.size());
    emissivity_[i].resize(nu_grid.size());
  }
  compton_opac.resize(grid->n_zones);
  photoion_opac.resize(grid->n_zones);
  

  // initialize time
  t_now_ = g->t_now;

  // initialize particles
  int n_parts = params_->getScalar<int>("particles_n_initialize");
  initialize_particles(n_parts);

  // allocate memory for core emission
  this->core_emis.resize(nu_grid.size());
 
}


//------------------------------------------------------------
// take a transport time step 
//------------------------------------------------------------
void transport::step(double dt)
{
  // nominal time for iterative calc is 1
  if (this->steady_state) dt = 1;
  
  // calculate opacities
  set_opacity();

  // emit new particles
  emit_particles(dt);

  // clear the tallies of the radiation quantities in each zone
  wipe_radiation();

  // Propagate the particles
  int n_active = particles.size();
  int n_escape = 0;
  std::list<particle>::iterator pIter = particles.begin();
  while (pIter != particles.end())
  {
    ParticleFate fate = propagate(*pIter,dt);
    if (fate == escaped) n_escape++;
    if ((fate == escaped)||(fate == absorbed)) pIter = particles.erase(pIter);
    else pIter++;
  }

  // calculate percent particles escaped, and rescale if wanted
  double per_esc = (1.0*n_escape)/(1.0*n_active);
  if ((verbose)&&(steady_state)) {
    cout << "# Percent particles escaped = " << 100.0*per_esc << "\n";
    optical_spectrum.rescale(1.0/per_esc); }

  // normalize and MPI combine radiation tallies
  reduce_radiation(dt);

  // solve for T_gas structure if radiative eq. applied
  if (radiative_eq) solve_eq_temperature();
   
  // advance time step
  if (!steady_state) t_now_ += dt;
}



//--------------------------------------------------------
// Propagate a single monte carlo particle until
// it  escapes, is absorbed, or the time step ends
//--------------------------------------------------------
ParticleFate transport::propagate(particle &p, double dt)
{
  enum ParticleEvent {scatter, boundary, tstep};
  ParticleEvent event;
  
  VoigtProfile voigt;

  // To be sure, get initial position of the particle 
  ParticleFate  fate = moving;
  p.ind = grid->get_zone(p.x);
  if (p.ind == -1) {return absorbed;}
  if (p.ind == -2) {return  escaped;}
  
  // time of end of timestep
  double tstop = t_now_ + dt;

  // pointer to current zone
  zone *zone = &(grid->z[p.ind]);

  // propagate until this flag is set
  while (fate == moving)
  {
    // set pointer to current zone
    zone = &(grid->z[p.ind]);
    
    // maximum step size inside zone
    double d_bn = step_size_*grid->zone_min_length(p.ind);

    // determine the doppler shift from comoving to lab
    //double dshift = dshift_comoving_to_lab(&p);
    double dshift = dshift_lab_to_comoving(&p);
    double nu_cmf = p.nu*dshift;

    // get continuum opacity and absorption fraction (epsilon)
    double continuum_opac_cmf, eps_absorb_cmf;
    get_opacity(p,dshift,continuum_opac_cmf,eps_absorb_cmf);

    // ------------------------------------------------
    // get total line opacity
    // ------------------------------------------------
    double line_opac_cmf = 0; 
    if (use_detailed_lines_)
    {
      // line width beta = v/c
      double beta_line;
      if (line_velocity_width_ != 0) 
	beta_line = line_velocity_width_/pc::c;
      else
	beta_line = sqrt(2*pc::k*zone->T_gas/pc::m_p)/pc::c;
      
      // find the nearest line
      int iLine =  upper_bound(line_nu_.begin(), line_nu_.end(), nu_cmf)
	- line_nu_.begin();
 
      // look for lines in front, add to opacity
      for (int i = iLine;i < n_lines_;i++)
      {
	double dnu_t   = line_nu_[i]*beta_line;
	if (!line_velocity_width_) dnu_t = dnu_t/line_sqrt_Mion_[i];
	double line_x  = (line_nu_[i] - nu_cmf)/dnu_t;
	if (line_x > 5) break;
	double line_profile = voigt.getProfile(line_x,1e-4)/dnu_t;
	line_opac_cmf += line_profile*line_opacity_[p.ind][i];
      }
      // look for lines behind, add to opacity
      for (int i = iLine-1;i >=0; i--)
      {
	double dnu_t   = line_nu_[i]*beta_line;
	if (!line_velocity_width_) dnu_t = dnu_t/line_sqrt_Mion_[i];
	double line_x  = (nu_cmf - line_nu_[i])/dnu_t*line_sqrt_Mion_[i];
	if (line_x > 5) break;
	double line_profile = voigt.getProfile(line_x,1e-4)/dnu_t;
	line_opac_cmf += line_profile*line_opacity_[p.ind][i];
      }


    }

    // convert opacity from comoving to lab frame for the purposes of 
    // determining the interaction distance in the lab frame
    // This corresponds to equation 90.8 in Mihalas&Mihalas. You multiply 
    // the comoving opacity by nu_0 over nu, which is why you
    // multiply by dshift instead of dividing by dshift here
    double tot_opac_cmf      = continuum_opac_cmf + line_opac_cmf;
    double tot_opac_labframe = tot_opac_cmf*dshift;
    
    // random optical depth to next interaction
    double tau_r = -1.0*log(1 - gsl_rng_uniform(rangen));
    
    // step size to next interaction event
    double d_sc  = tau_r/tot_opac_labframe;
    if (tot_opac_labframe == 0) d_sc = std::numeric_limits<double>::infinity();
    if (d_sc < 0) cout << "ERROR: negative interaction distance!\n";
  
    // find distance to end of time step
    double d_tm = (tstop - p.t)*pc::c;
    // if iterative calculation, let all particles escape
    if (this->steady_state) d_tm = std::numeric_limits<double>::infinity();

    // find out what event happens (shortest distance)
    double this_d;
    if ((d_sc < d_bn)&&(d_sc < d_tm))
      {event = scatter;    this_d = d_sc;}
    else if (d_bn < d_tm)
      {event = boundary;   this_d = d_bn;}
    else 
      {event = tstep;      this_d = d_tm; }

    //    std::cout << d_sc << " " << d_bn << " " << nu_cmf << " " << line_opac_cmf << "\n";

    // tally in contribution to zone's radiation energy (both *lab* frame)
    double this_E = p.e*this_d; 
    zone->e_rad  += this_E; 

    // store absorbed energy in *comoving* frame 
    // (will turn into rate by dividing by dt later)
    // Extra dshift definitely needed here (two total)
    zone->e_abs  += this_E*dshift*(continuum_opac_cmf)*eps_absorb_cmf*dshift; 

    // put back in radiation force tally here
    // fx_rad =

    // move particle the distance
    p.x[0] += this_d*p.D[0];
    p.x[1] += this_d*p.D[1];
    p.x[2] += this_d*p.D[2]; 
    // advance the time
    p.t = p.t + this_d/pc::c;

    // Find position of the particle now
    p.ind = grid->get_zone(p.x);
    if (p.ind == -1) fate = absorbed;
    if (p.ind == -2) fate = escaped;

    // check for inner boundary absorption
    if (p.r() < r_core_)  {fate = absorbed;}

    // ---------------------------------
    // do an interaction event
    // ---------------------------------
    if (fate == moving) 
    {
      if (event == scatter)   
	fate = do_scatter(&p,eps_absorb_cmf);  
      else if (event == tstep) 
	fate = stopped;
    }
  }

  // Add escaped photons to output spectrum
  if (fate == escaped) 
  {
    // account for light crossing time, relative to grid center
    double xdot = p.x[0]*p.D[0] + p.x[1]*p.D[1] + p.x[2]*p.D[2];
    double t_obs = p.t - xdot/pc::c;
    if (p.type == photon)   optical_spectrum.count(t_obs,p.nu,p.e,p.D);
    if (p.type == gammaray) gamma_spectrum.count(t_obs,p.nu,p.e,p.D);
  }
  return fate;
}



//--------------------------------------------------------------
// Finalize the spectrum and output it
//--------------------------------------------------------------
void transport::output_spectrum(int it)
{
  std::stringstream ss;
  if (it < 0)  ss << "_0" << it;
  else ss << "_" << it;
  string base = ss.str();

  string specname = params_->getScalar<string>("spectrum_name");
  if (specname != "") 
    {
    optical_spectrum.set_name(specname + base + ".dat");
    optical_spectrum.MPI_average();
    if (verbose) optical_spectrum.print();
  }
  optical_spectrum.wipe();

  string gamname = params_->getScalar<string>("gamma_name");
  if (gamname != "") 
  {
    gamma_spectrum.set_name(gamname + base + ".dat");
    gamma_spectrum.MPI_average();
    if (verbose) gamma_spectrum.print();
  }
  gamma_spectrum.wipe();

}

void transport::output_spectrum()
{
  output_spectrum(0);
}


void transport::wipe_radiation()
{
  for (int i=0;i<grid->n_zones;i++) 
  {
    grid->z[i].e_rad  = 0;
    grid->z[i].e_abs  = 0;
    //grid->z[i].fx_rad = 0;
    //grid->z[i].fy_rad = 0;
    //grid->z[i].fz_rad = 0;
  }
}


//------------------------------------------------------------
// Combine the radiation tallies in all zones
// from all processors using MPI 
//------------------------------------------------------------
 void transport::reduce_radiation(double dt)
{
  // properly normalize the radiative quantities
  for (int i=0;i<grid->n_zones;i++) 
  {
    double vol = grid->zone_volume(i);
    grid->z[i].e_rad   /= vol*pc::c*dt;
    grid->z[i].e_abs   /= vol*dt; 
    //grid->z[i].fx_rad  /= vol*pc::c*dt; 
    //grid->z[i].fy_rad  /= vol*pc::c*dt;
    //grid->z[i].fz_rad  /= vol*pc::c*dt;
  }


//   vector<real> send, receive;
//   int my_begin, my_end, size;

//   //-- EACH PROCESSOR GETS THE REDUCTION INFORMATION IT NEEDS
//   for(int proc=0; proc<MPI_nprocs; proc++){

//     // set the begin and end indices so a process covers range [begin,end)
//     my_begin = ( proc==0 ? 0 : my_zone_end[proc-1] );
//     my_end = my_zone_end[proc];

//     // set the computation size and create the send/receive vectors
//     size = my_end - my_begin;
//     send.resize(size);
//     receive.resize(size);

//     // reduce e_rad
//     for(int i=my_begin; i<my_end; i++) send[i-my_begin] = grid->z[i].e_abs;
//     MPI_Reduce(&send.front(), &receive.front(), size, MPI_real, MPI_SUM, proc, MPI_COMM_WORLD);
//     for(int i=my_begin; i<my_end; i++) grid->z[i].e_abs = receive[i-my_begin] / (real)MPI_nprocs;
//     }

//     // TODO - need to put in other quantities...
//   }
}
