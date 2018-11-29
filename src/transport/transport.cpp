
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

  tstr = get_system_time();
  set_opacity();
  tend = get_system_time();
  if (verbose) cout << "# Calculated opacities   (" << (tend-tstr) << " secs) \n";

  // mpi combine the opacities calculated
  tstr = get_system_time();
  reduce_opacities();
  reduce_Lthermal();
  tend = get_system_time();
  if (verbose) cout << "# Communicated opacities (" << (tend-tstr) << " secs) \n";

  if (use_ddmc_) compute_diffusion_probabilities(dt);

  // clear the tallies of the radiation quantities in each zone
  wipe_radiation();

  set_eps_imc();

  // emit new particles
  tstr = get_system_time();
  emit_particles(dt);

  // Propagate the particles
  int n_active = particles.size();
  int n_escape = 0;
  std::list<particle>::iterator pIter = particles.begin();
  while (pIter != particles.end())
  {
    int ddmc_zone = 0;
    if (use_ddmc_)
    {
      double dx;
      grid->get_zone_size(pIter->ind,&dx);
      double ztau = planck_mean_opacity_[pIter->ind]*dx;
      if (ztau > ddmc_tau_) ddmc_zone = 1;
    }
    ParticleFate fate;
    if (ddmc_zone) fate = discrete_diffuse(*pIter,dt);
    else fate = propagate(*pIter,dt);

    if (fate == escaped) n_escape++;
    if ((fate == escaped)||(fate == absorbed)) pIter = particles.erase(pIter);
    else pIter++;
  }

  // calculate percent particles escaped, and rescale if wanted
  if (steady_state)
  {
    double per_esc = (1.0*n_escape)/(1.0*n_active);
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
  if (verbose) cout << "# Propagated particles          (" << (tend-tstr) << " secs) \n";

  // normalize and MPI combine radiation tallies
  tstr = get_system_time();
  reduce_radiation(dt);
  tend = get_system_time();
  if (verbose) cout << "# Reduced radiation (" << (tend-tstr) << " secs) \n";

  // solve for T_gas structure if radiative eq. applied
  if (radiative_eq)
  {
    tstr = get_system_time();
    solve_eq_temperature();
    tend = get_system_time();
    if (verbose) cout << "# Calculated Radiative Equilib (" << (tend-tstr) << " secs) \n";
  }
  // advance time step
  if (!steady_state) t_now_ += dt;

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
  double tstop = t_now_ + dt;

  // pointer to current zone
  zone *zone = &(grid->z[p.ind]);

  // propagate until this flag is set
  while (fate == moving)
  {
    // set pointer to current zone
    assert(p.ind >= 0);
    zone = &(grid->z[p.ind]);

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
    double tau_r = -1.0*log(1 - gsl_rng_uniform(rangen));

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
      zone->e_abs  += this_E*dshift*(continuum_opac_cmf)*eps_absorb_cmf*dshift;
      if (store_Jnu_) J_nu_[p.ind][i_nu] += this_E;
      else J_nu_[p.ind][0] += this_E;
      //std::cout << p.ind << " " << i_nu << " " << p.e << " " << this_E << " " << J_nu_[p.ind][i_nu] << "\n";
    }

     // radiation force
     zone->fx_rad += this_E*dshift*continuum_opac_cmf*p.D[0] * dshift; // Extra dshift definitely needed here (two total)
     zone->fy_rad += this_E*dshift*continuum_opac_cmf*p.D[1] * dshift;
     zone->fz_rad += this_E*dshift*continuum_opac_cmf*p.D[2] * dshift;

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
   //     if (gsl_rng_uniform(rangen) > 0.38)
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
//      std::cout << "dists: "<< d_sc << "\t" << d_bn << "\t" << d_tm << "\t" << d_nu << "\n";
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
