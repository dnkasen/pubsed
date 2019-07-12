
#include <math.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>
#include <vector>
#include <cassert>
#include <list>
#include <algorithm>
#include <ctime>

#include "transport.h"
#include "ParameterReader.h"
#include "physical_constants.h"

using std::cout;
using std::cerr;
using std::endl;
namespace pc = physical_constants;

//------------------------------------------------------------
// take a transport time step
//------------------------------------------------------------
void transport::step(double dt)
{
  // nominal time for iterative calc is 1
  if (this->steady_state) dt = 1;

  double tend,tstr;
  double get_system_time(void);

  // calculate the opacities
  set_opacity(dt); // need to pass dt for computing implicit monte carlo factors

  // mpi combine the opacities calculated
  tstr = get_system_time();
  reduce_opacities();
  reduce_Lthermal();
  tend = get_system_time();
  if (verbose) cout << "# Communicated opacities (" << (tend-tstr) << " secs) \n";

  if (use_ddmc_) compute_diffusion_probabilities(dt);

  // clear the tallies of the radiation quantities in each zone
  wipe_radiation();

  // emit new particles
  tstr = get_system_time();
  emit_particles(dt);

  // Propagate the particles
  int n_active = particles.size();
  int n_particles = particles.size();

  #pragma omp parallel for schedule(guided)
  for(int i=0; i<n_particles; i++)
  {
    // propagate particles
    particles[i].fate = propagate(particles[i],dt);

    // Add escaped photons to output spectrum
    if (particles[i].fate == escaped)
    {
      // account for light crossing time, relative to grid center
        double t_obs = particles[i].t - particles[i].x_dot_d()/pc::c;
        if (particles[i].type == photon)
          optical_spectrum.count(t_obs,particles[i].nu,particles[i].e,particles[i].D);
        if (particles[i].type == gammaray)
          gamma_spectrum.count(t_obs,particles[i].nu,particles[i].e,particles[i].D);
      }
  }

  // Remove escaped and absorbed particles from the particle vector
  int n_escaped = clean_up_particle_vector();

  // calculate percent particles escaped, and rescale if wanted
  if (steady_state)
  {
    double per_esc = (1.0*n_escaped)/(1.0*n_active);
    if (core_fix_luminosity_)
    {
      if (verbose)
        cout << "# Percent particles escaped = " << 100.0*per_esc << " (rescaling)\n";

      double fac = 1.0/per_esc;
      optical_spectrum.rescale(fac);
      for (int i=0;i<grid->n_zones;++i)
      {
        grid->z[i].e_rad *= fac;
        if (store_Jnu_)
         for (size_t j=0;j<nu_grid.size();++j)
            J_nu_[i][j] *= fac;
      }
    }
    else {
      if (verbose)
        cout << "# Percent particles escaped = " << 100.0*per_esc << " (not rescaling)\n";}
  }

  tend = get_system_time();
  if (verbose) cout << "# Propagated particles   (" << (tend-tstr) << " secs) \n";

  // normalize and MPI combine radiation tallies
  tstr = get_system_time();
  reduce_radiation(dt);
  tend = get_system_time();
  if (verbose) cout << "# Communicated radiation (" << (tend-tstr) << " secs) \n";

  if ((fix_Tgas_during_transport_ == 0 && solve_Tgas_with_updated_opacities_ == 0))
  {
    tstr = get_system_time();
    solve_eq_temperature();
    tend = get_system_time();
    if (verbose) cout << "# Calculated temperature (" << (tend-tstr) << " secs) \n";
  }


  // advance time step
  if (!steady_state) t_now_ += dt;

  if (first_step_) first_step_ = 0.;

}



//--------------------------------------------------------
// little local helper function to get the current
// time for timing
//--------------------------------------------------------
double get_system_time()
{
#ifdef MPI_PARALLEL
  return MPI_Wtime();
#else
  return  ((double)clock())  / (double)CLOCKS_PER_SEC;
#endif

}


//--------------------------------------------------------
// Loop over the vector of particles
// and remove those that are either escaped or absorbed
// Returns the number of particles that escaped
//--------------------------------------------------------
int transport::clean_up_particle_vector()
{

  // do nothing to an empty particle vector
  if (particles.size() == 0) return 0;

  int n_escaped = 0;
  int i=0;
  while (true)
  {
    // remove particles from back until we have one that is allive
    while ((particles.back().fate == escaped)||(particles.back().fate == absorbed))
    {
      if (particles.back().fate == escaped) n_escaped++;
        particles.pop_back();
        if (particles.size() == 0) break;
    }

   // see if we've finished with all particles
   if (i >= particles.size()) break;

   // check if we should remove this particle
   if (particles[i].fate == escaped) n_escaped++;
   if ((particles[i].fate == escaped)||(particles[i].fate == absorbed)){
     particles[i] = particles.back();
     particles.pop_back();
   }
   i = i+1;
 }
  return n_escaped;
}

//--------------------------------------------------------
// Propagate a particle until either the
// time step ends at a time tstop
// or the particle escapes or is absorbed.
// Returns this fate of the particle
//--------------------------------------------------------
ParticleFate transport::propagate(particle &p, double dt)
{
  // To be sure, get initial position of the particle
  p.ind = grid->get_zone(p.x);

  if (p.ind == -1) {return absorbed;}
  if (p.ind == -2) {return  escaped;}

  // time of end of timestep
  double tstop = t_now_ + dt;

  ParticleFate  fate = moving;
  while (fate == moving)
  {
    // check if we are in DDMC zone
    int in_ddmc_zone = 0;
    if (use_ddmc_)
      if ((ddmc_use_in_zone_[p.ind])&&(p.type == photon))
        in_ddmc_zone = 1;

    if (in_ddmc_zone){
      if(use_ddmc_ == 1)
	fate = discrete_diffuse_IMD(p, tstop);
      else if(use_ddmc_ == 2)
	fate = discrete_diffuse_DDMC(p, tstop);
      else if(use_ddmc_ == 3)
	fate = discrete_diffuse_RandomWalk(p, tstop);
      else{
	cout << "Invalid diffusion method" << endl;
	exit(1);
      }
    }
    else
        fate = propagate_monte_carlo(p, tstop);
  }

return fate;

}

//--------------------------------------------------------
// Propagate a single monte carlo particle until
// it  escapes, is absorbed, or the time step ends
//--------------------------------------------------------
ParticleFate transport::propagate_monte_carlo(particle &p, double tstop)
{
  enum ParticleEvent {scatter, boundary, tstep};
  ParticleEvent event;

  ParticleFate  fate = moving;
  while (fate == moving)
  {
    // set pointer to current zone
    assert(p.ind >= 0);
    zone *zone = &(grid->z[p.ind]);

    // check if we have moved into a DDMC zone
    if (use_ddmc_)
      if ((ddmc_use_in_zone_[p.ind])&&(p.type == photon))
        return moving;

    // printout for debug
    //std::cout << " p = " << sqrt(p.x[0]*p.x[0] + p.x[1]*p.x[1]);
    //std::cout << " x = " << p.x[0] << " y = " << p.x[1] << "; z = " << p.x[2]"\n";
    //std::cout << "D = " << p.D[0] << ", " << p.D[1] << ", " << p.D[2] << "\n";

    // get distance and index to the next zone boundary
    double d_bn = 0;
    int new_ind = grid->get_next_zone(p.x,p.D,p.ind,r_core_,&d_bn);

    // determine the doppler shift from comoving to lab
    double dshift = dshift_lab_to_comoving(&p);

    // get continuum opacity and absorption fraction (epsilon)
    double continuum_opac_cmf,eps_absorb_cmf;
    int i_nu = get_opacity(p,dshift,continuum_opac_cmf,eps_absorb_cmf);

    // check for distance to next frequency bin
    // nushift = nu*(dvds*l)/c --> l = nushift/nu*c/dvds
    double d_nu = nu_grid.delta(i_nu)/p.nu*pc::c/p.dvds;
    if (p.dvds == 0) d_nu = std::numeric_limits<double>::infinity();
    if (d_nu < 0) d_nu = -1*d_nu;
    if (d_nu < d_bn)
    {
      d_bn = d_nu;
      new_ind = p.ind;
    }

    if (d_bn == 0) std::cout << "zerob\n";

    // convert opacity from comoving to lab frame for the purposes of
    // determining the interaction distance in the lab frame
    // This corresponds to equation 90.8 in Mihalas&Mihalas. You multiply
    // the comoving opacity by nu_0 over nu, which is why you
    // multiply by dshift instead of dividing by dshift here
    double tot_opac_cmf      = continuum_opac_cmf;
    double tot_opac_labframe = tot_opac_cmf*dshift;

    // random optical depth to next interaction
    double tau_r = -1.0*log(1 - rangen.uniform());

    // step size to next interaction event
    double d_sc  = tau_r/tot_opac_labframe;
    if (tot_opac_labframe == 0) d_sc = std::numeric_limits<double>::infinity();
    if (d_sc < 0)
      cerr << "ERROR: negative interaction distance! " << d_sc << " " << p.nu << " " << dshift << " " <<
        tot_opac_labframe  << endl;
    //std::cout << "dists: "<< d_sc << "\t" << d_bn << "\t" << tot_opac_labframe << "\n";

    // find distance to end of time step
    double d_tm = (tstop - p.t)*pc::c;
    // if iterative calculation, give infinite time for particle escape
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

    // store absorbed energy in *comoving* frame
    // (will turn into rate by dividing by dt later)
    // Extra dshift definitely needed here (two total)
    // don't add gamma-rays here (they would be separate)
    if (p.type == photon)
    {
      #pragma omp atomic
      zone->e_abs  += this_E*dshift*(continuum_opac_cmf)*eps_absorb_cmf*dshift * zone->eps_imc;
      if (store_Jnu_)
	     #pragma omp atomic
	      J_nu_[p.ind][i_nu] += this_E;
      else
	     #pragma omp atomic
	      J_nu_[p.ind][0] += this_E;
    }

     // tally radiation force
     // Extra dshift definitely needed here (two total)
    #pragma omp atomic
    zone->fx_rad += this_E*dshift*continuum_opac_cmf*p.D[0] * dshift;
    #pragma omp atomic
    zone->fy_rad += this_E*dshift*continuum_opac_cmf*p.D[1] * dshift;
    #pragma omp atomic
    zone->fz_rad += this_E*dshift*continuum_opac_cmf*p.D[2] * dshift;
    // radial radiation force
    double rr = sqrt(p.x[0]*p.x[0] + p.x[1]*p.x[1] + p.x[2]*p.x[2]);
    double xdotD = p.x[0]*p.D[0] + p.x[1]*p.D[1] + p.x[2]*p.D[2];
    #pragma omp atomic
    zone->fr_rad += this_E*dshift*continuum_opac_cmf*xdotD/rr * dshift;

    // move particle the distance
    p.x[0] += this_d*p.D[0];
    p.x[1] += this_d*p.D[1];
    p.x[2] += this_d*p.D[2];
    // advance the time
    p.t = p.t + this_d/pc::c;

    // ---------------------------------
    // do a boundary event
    // ---------------------------------
    if (event == boundary)
    {
      // inner boundary hit
      if (new_ind == -1)
      {
        if (boundary_in_reflect_)
        {
          // flip direction
          p.D[0] *= -1;
          p.D[1] *= -1;
          p.D[2] *= -1;
        }
        else
        {
          p.ind = new_ind;
          fate = absorbed;
        }
      }
      // outer boundary hit
      else if (new_ind == -2)
      {
        if (boundary_out_reflect_)
        {
          // flip direction
          p.D[0] *= -1;
          p.D[1] *= -1;
          p.D[2] *= -1;
        }
        else
        {
          p.ind = new_ind;
          fate = escaped;
        }
      }
      // just another cell
      else p.ind = new_ind;
    }

    // ---------------------------------
    // do an interaction event
    // ---------------------------------
    else if (event == scatter)
    {
       //eps_absorb_cmf = 0;

      if (fate == moving)
      {
 //        fate = do_scatter(&p,eps_absorb_cmf);
   //     if (rangen.uniform() > 0.38)
    //    fate = absorbed;
        //else fate = do_scatter(&p,0);
         // debug
        fate = do_scatter(&p,eps_absorb_cmf);
        //if (p.nu*pc::h < pc::rydberg) fate = escaped;
       }
    }

    // ---------------------------------
    // do an end of timestep event
    // ---------------------------------
    else if (event == tstep)
    {
       fate = stopped;
    }
   }

  return fate;
}

transport::~transport() {
  if (src_MPI_block)
    delete[] src_MPI_block;
  if (src_MPI_zones)
    delete[] src_MPI_zones;
  if (dst_MPI_block)
    delete[] dst_MPI_block;
  if (dst_MPI_zones)
    delete[] dst_MPI_zones;
}
