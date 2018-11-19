#include <math.h>
#include <gsl/gsl_rng.h>
#include "transport.h"
#include "radioactive.h"
#include "physical_constants.h"

namespace pc = physical_constants;


// ------------------------------------------------------
// Propagate a particle using discrete diffusion
// This is only implemented for 1D spherical so far
// ------------------------------------------------------
ParticleFate transport::discrete_diffuse(particle &p, double dt)
{
  int stop = 0;

  double dx;
  grid->get_zone_size(p.ind,&dx);
  
  // set this zone
  while (!stop)
  {
    // pointer to current zone
    zone *zone = &(grid->z[p.ind]);

    // add in tally of absorbed and total radiation energy
    zone->e_abs += p.e*ddmc_P_abs_[p.ind];
    //zone->e_rad += p.e*ddmc_P_stay_[p.ind];
    J_nu_[p.ind][0] += p.e*ddmc_P_stay_[p.ind]*dt*pc::c;
//    std::cout <<  p.ind << " " << p.e*ddmc_P_stay_[p.ind] << "\n";

    // total probability of diffusing in some direction
    double P_diff = ddmc_P_up_[p.ind]  + ddmc_P_dn_[p.ind];
    double P_stay = ddmc_P_abs_[p.ind] + ddmc_P_stay_[p.ind];

    // see if diffuse
    double r1 = gsl_rng_uniform(rangen);
    if (r1 < P_diff)
    {
      // find dimension to diffuse in...
      // move up or down in this dimension
      double r2 = gsl_rng_uniform(rangen);
      if (r2 < ddmc_P_up_[p.ind]/P_diff)
      {
        p.ind++;
        double rr = p.r();
        p.x[0] += p.x[0]/rr*dx;
        p.x[1] += p.x[1]/rr*dx;
        p.x[2] += p.x[2]/rr*dx;
      }
      else
      {
        p.ind--;
        double rr = p.r();
        p.x[0] -= p.x[0]/rr*dx;
        p.x[1] -= p.x[1]/rr*dx;
        p.x[2] -= p.x[2]/rr*dx;
      }
    }
    // don't diffuse
    else
    {
      // check for absorption
      double r2 = gsl_rng_uniform(rangen);
      double f_abs = ddmc_P_abs_[p.ind]/P_stay;
      // see if absorbed
      if (r2 < f_abs) { return absorbed;}

      // advect it
      double zone_vel[3], dvds;
      grid->get_velocity(p.ind,p.x,p.D,zone_vel, &dvds);
      p.x[0] += zone_vel[0]*dt;
      p.x[1] += zone_vel[1]*dt;
      p.x[2] += zone_vel[2]*dt;
      p.ind = grid->get_zone(p.x);

      // adiabatic loses
      //p.e *= (1 - zone->diff_v);
      stop = 1;
    }

    p.t += dt;

    // check for escape
    if (p.ind < 0)         {return absorbed;}
    if (p.ind > grid->n_zones - 1)	  {return escaped; }

  }
  return stopped;

}



// ------------------------------------------------------
// Calcuate the probabilities of diffusion
// for now this only works in 1D spherical coords
// ------------------------------------------------------
void transport::compute_diffusion_probabilities(double dt)
{
  int nz = grid->n_zones;

  // debug hard coded for now
  double grey_opac = 1.0;

  for (int i=0;i<nz;i++)
  {
    double dx;
    grid->get_zone_size(i,&dx);

    // indices of adjacent zones
    int ip = i+1;
    if (ip == nz) ip = i;
    int im = i-1;
    if (im < 0)   im = 0;

    // diffusion probability in zone and adjacent zones
    double Dj0 = pc::c/(3.0*grey_opac*grid->z[i].rho);
    double Djp = pc::c/(3.0*grey_opac*grid->z[ip].rho);
    double Djm = pc::c/(3.0*grey_opac*grid->z[im].rho);

    double Dh_up = 2*dx*(Dj0*Djp)/(Dj0*dx + Djp*dx);
    double Dh_dn = 2*dx*(Dj0*Djm)/(Dj0*dx + Djm*dx);

    ddmc_P_up_[i] = (dt/dx)*(Dh_up/dx);
    ddmc_P_dn_[i] = (dt/dx)*(Dh_dn/dx);
    // boundary condition -- don't diffuse inward at innermost cell
    if (i==0)
      ddmc_P_dn_[i] = 0;

    // advection probability
    ddmc_P_adv_[i] = 0.0;
    ddmc_P_abs_[i] = 0.0; //pc::c*grey_opac*grid->z[i].rho; //*epsilon*grid->z[i0].eps_imc;
    ddmc_P_stay_[i] = 0.0;

    // get normalization
    double norm = 1 + ddmc_P_adv_[i] + ddmc_P_abs_[i];
    norm += ddmc_P_up_[i] +  ddmc_P_dn_[i];
    ddmc_P_adv_[i]  /= norm;
    ddmc_P_abs_[i]  /= norm;
    ddmc_P_up_[i]   /= norm;
    ddmc_P_dn_[i]   /= norm;
    ddmc_P_stay_[i] = 1.0/norm;

    //std::cout <<   ddmc_P_up_[i] << "\t" <<   ddmc_P_dn_[i] << "\t" << ddmc_P_adv_[i] << "\t" << ddmc_P_abs_[i] << "\t" << ddmc_P_stay_[i] << "\n";
  }
}
