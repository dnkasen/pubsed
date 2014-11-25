#include <mpi.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <limits>
#include <vector>
#include <cassert>
#include "transport.h"
#include "physical_constants.h"
#include "Lua.h"

using std::cout;
namespace pc = physical_constants;


//----------------------------------------------------------------------------
// Initialize the transport module
// Includes setting up the grid, particles,
// and MPI work distribution
//----------------------------------------------------------------------------
void transport::init(Lua* l, grid_general *g)
{ 
  this->lua  = l;
  this->grid = g;

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
  this->radiative_eq = lua->scalar<int>("radiative_eq");
  this->steady_state = (lua->scalar<int>("steady_iterate") > 0);
  this->step_size    = lua->scalar<double>("step_size");
  this->max_total_particles = lua->scalar<int>("max_total_particles");

  // initialize the frequency grid  
  std::vector<double> nu_dims = lua->vector<double>("nu_grid");
  if (nu_dims.size() != 3) {
    cout << "# improperly defined nu_grid; need {nu_1, nu_2, dnu}; exiting\n";
    exit(1); }
  nu_grid.init(nu_dims[0],nu_dims[1],nu_dims[2]);
    
  // intialize output spectrum
  std::vector<double>stg = lua->vector<double>("spec_time_grid");
  std::vector<double>sng = lua->vector<double>("spec_nu_grid");
  int nmu  = lua->scalar<int>("spec_n_mu");
  int nphi = lua->scalar<int>("spec_n_phi");
  optical_spectrum.init(stg,sng,nmu,nphi);
  std::vector<double>gng = lua->vector<double>("gamma_nu_grid");  
  gamma_spectrum.init(stg,sng,nmu,nphi);

  // initalize and allocate space for opacities
  this->epsilon   = lua->scalar<double>("epsilon");
  this->grey_opac = lua->scalar<double>("grey_opacity");
  abs_opac.resize(grid->n_zones);
  scat_opac.resize(grid->n_zones);
  emis.resize(grid->n_zones);
  for (int i=0; i<abs_opac.size();  i++)
  {
    abs_opac[i].resize(nu_grid.size());
    scat_opac[i].resize(nu_grid.size());
    emis[i].resize(nu_grid.size());
  }
  compton_opac.resize(grid->n_zones);
  photoion_opac.resize(grid->n_zones);
  
  // initialize nlte_gas
  std::string atomdata = lua->scalar<string>("atomic_data");  
  gas.init(atomdata,grid->elems_Z,grid->elems_A,nu_grid);
  std::string fuzzfile = lua->scalar<string>("fuzzline_file");  
  int nl = gas.read_fuzzfile(fuzzfile);
  if (verbose) std::cout << "# read fuzzfile " << fuzzfile << "; " << 
		 nl << " lines used\n";

  this->t_now = g->t_now;

  // initialize particles
  int n_parts = lua->scalar<int>("init_particles");
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
  std::vector<particle>::iterator pIter = particles.begin();
  int n_active = particles.size();
  int n_escape = 0;
  while (pIter != particles.end())
  {
    ParticleFate fate = propagate(*pIter,dt);
    if (fate == escaped) n_escape++;
    if ((fate == escaped)||(fate == absorbed)) particles.erase(pIter);
    else pIter++;
  }
  double per_esc = (100.0*n_escape)/n_active;
  if ((verbose)&&(steady_state)) {
    cout << "# Percent escaped = " << per_esc << "\n";
    optical_spectrum.rescale(1.0/per_esc); }

  // normalize and MPI combine radiation tallies
  reduce_radiation(dt);

  // solve for T_gas structure if radiative eq. applied
  if (radiative_eq) solve_eq_temperature();
   
  // advance time step
  if (!steady_state) t_now += dt;
}



//--------------------------------------------------------
// Propagate a single monte carlo particle until
// it  escapes, is absorbed, or the time step ends
//--------------------------------------------------------
ParticleFate transport::propagate(particle &p, double dt)
{
  enum ParticleEvent {scatter, boundary, tstep};
  ParticleEvent event;
  
  // To be sure, get initial position of the particle 
  ParticleFate  fate = moving;
  p.ind = grid->get_zone(p.x);
  if (p.ind == -1) {return absorbed;}
  if (p.ind == -2) {return  escaped;}
  
  // time of end of timestep
  double tstop = t_now + dt;

  // pointer to current zone
  zone *zone = &(grid->z[p.ind]);

  // propagate until this flag is set
  while (fate == moving)
  {
    // set pointer to current zone
    zone = &(grid->z[p.ind]);
    
    // maximum step size inside zone
    double d_bn = this->step_size*grid->zone_min_length(p.ind);

    // determine the doppler shift from comoving to lab
    //double dshift = dshift_comoving_to_lab(&p);
    double dshift = dshift_lab_to_comoving(&p);

    // get local opacity and absorption fraction (epsilon)
    double opac, eps;
    get_opacity(p,dshift,opac,eps);
    
    // convert opacity from comoving to lab frame for the purposes of 
    // determining the interaction distance in the lab frame
    // This corresponds to equation 90.8 in Mihalas&Mihalas. You multiply 
    // the comoving opacity by nu_0 over nu, which is why you
    // multiply by dshift instead of dividing by dshift here
    opac = opac*dshift;

    // random optical depth to next interaction
    double tau_r = -1.0*log(1 - gsl_rng_uniform(rangen));
    
    // step size to next interaction event
    double d_sc  = tau_r/opac;
    if (opac == 0) d_sc = std::numeric_limits<double>::infinity();
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

    // tally in contribution to zone's radiation energy (both *lab* frame)
    double this_E = p.e*this_d; 
    zone->e_rad  += this_E; 

    // shift opacity back to comoving frame for energy and momentum exchange. 
    // Radiation energy is still lab frame
    opac = opac / dshift;

    // store absorbed energy in *comoving* frame 
    // (will turn into rate by dividing by dt later)
    // Extra dshift definitely needed here (two total)
    zone->e_abs  += this_E*dshift*(opac)*eps*dshift; 

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
    if (p.r() < this->r_core)  {fate = absorbed;}

    // ---------------------------------
    // do an interaction event
    // ---------------------------------
    if (fate == moving) 
    {
      if (event == scatter)  fate = do_scatter(&p,eps);
      else if (event == tstep) fate = stopped;
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
  string base = "";
  if (it > 0) base = "_" + std::to_string(it);
  string specname = lua->scalar<string>("spectrum_name");
  optical_spectrum.set_name(specname + "_optical" + base + ".spec");
  optical_spectrum.MPI_average();
  if (verbose) optical_spectrum.print();
  optical_spectrum.wipe();

  gamma_spectrum.set_name(specname + "_gamma" + base + ".spec");
  gamma_spectrum.MPI_average();
  if (verbose) gamma_spectrum.print();
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
