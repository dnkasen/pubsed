#include <math.h>
#include <gsl/gsl_rng.h>
#include "transport.h"
#include "radioactive.h"
#include "physical_constants.h"

namespace pc = physical_constants;


// ------------------------------------------------------
// Propagate a particle using the
// Implicit Monte Carlo Diffusion (IMD) approach.
// Reference: Gentile, J. of Comput. Physics 172, 543–571 (2001)
// This is only implemented for 1D spherical so far
// ------------------------------------------------------
ParticleFate transport::discrete_diffuse_IMD(particle &p, double dt)
{
  int stop = 0;

  double dx;
  grid->get_zone_size(p.ind,&dx);

  // for now incrementing whole time-stepping
  // this is not really correct
  p.t += dt;

  // set this zone
  while (!stop)
  {
    // find current zone and check for escape
    p.ind = grid->get_zone(p.x);
    if (p.ind == -1) {return absorbed;}
    if (p.ind == -2) {return escaped;}

    // pointer to current zone
    zone *zone = &(grid->z[p.ind]);

    // add in tally of absorbed and total radiation energy
    #pragma omp atomic
    zone->e_abs += p.e*ddmc_P_abs_[p.ind];
    //zone->e_rad += p.e*ddmc_P_stay_[p.ind];
    #pragma omp atomic
    J_nu_[p.ind][0] += p.e*ddmc_P_stay_[p.ind]*dt*pc::c;

    // total probability of diffusing in some direction
    double P_diff = ddmc_P_up_[p.ind]  + ddmc_P_dn_[p.ind];
    double P_stay = ddmc_P_abs_[p.ind] + ddmc_P_stay_[p.ind];

    // randomly choose whether to diffuse
    double r1 = rangen.uniform();
    if (r1 < P_diff)
    {
      // find dimension to diffuse in...
      // move up or down in this dimension
      double r2 = rangen.uniform();
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
    // else don't diffuse and stop
    else
    {
      // check for absorption
      double r2 = rangen.uniform();
      double f_abs = ddmc_P_abs_[p.ind]/P_stay;
      // see if absorbed
      if (r2 < f_abs) { return absorbed;}

      // advect it
      double zone_vel[3], dvds;
      grid->get_velocity(p.ind,p.x,p.D,zone_vel, &dvds);
      p.x[0] += zone_vel[0]*dt;
      p.x[1] += zone_vel[1]*dt;
      p.x[2] += zone_vel[2]*dt;

      // adiabatic loss (assumes small change in volume I think)
      p.e *= (1 - dvds*dt);
      stop = 1;
    }
  }

  // find current zone and check for escape
  p.ind = grid->get_zone(p.x);
  if (p.ind == -1) {return absorbed;}
  if (p.ind == -2) {return escaped;}
  return stopped;

}

// ------------------------------------------------------
// Propagate a particle using the
// Discrete Diffusion Monte Carlo (DDMC) approach.
// Reference: Densmore+, J. of Comput. Physics 222, 485-503 (2007)
// This is only implemented for 1D spherical also.
// ------------------------------------------------------
ParticleFate transport::discrete_diffuse_DDMC(particle &p, double tstop)
{
  // pointer to current zone
  zone *zone = &(grid->z[p.ind]);
  int nz = grid->n_zones;

  // initialize particle's timestamp
  double dt_remaining = tstop - p.t;

  // indices of adjacent zones
  int ii = p.ind;
  int ip = ii + 1;
  int im = ii - 1;

  if (ip == nz) ip = ii;
  if (im < 0)   im = 0;

  double dx;
  grid->get_zone_size(ii,&dx);

  double dxp1, dxm1;
  grid->get_zone_size(ip,&dxp1);
  grid->get_zone_size(im,&dxm1);

  // Gathering gas opacity
  double sigma_i   = planck_mean_opacity_[ii];
  double sigma_ip1 = planck_mean_opacity_[ip];
  double sigma_im1 = planck_mean_opacity_[im];

  double rcoords[3];
  grid->coordinates(ip,rcoords);
  double r_p = rcoords[0];
  grid->coordinates(im,rcoords);
  double r_m = rcoords[0];
  if (ii == 0) r_m = 0;
  double r_0 = 0.5*(r_p + r_m);

  // Compute left/right leakage opacity
  double sigma_leak_left  = (2.0/3.0/dx) * (1.0 / (sigma_i*dx + sigma_im1*dxm1))*r_m*r_m/r_0/r_0;
  double sigma_leak_right = (2.0/3.0/dx) * (1.0 / (sigma_i*dx + sigma_ip1*dxp1))*r_p*r_p/r_0/r_0;
  double sigma_leak_tot   = sigma_leak_left + sigma_leak_right;



  // While loop to sample
  while (dt_remaining > 0.0)
  {
    double xi = rangen.uniform();

    double d_stay, d_leak;

    // Distance until end of time step
    d_stay = pc::c * dt_remaining;

    // Distance to leakage event
    d_leak = -log(xi) / (sigma_leak_left + sigma_leak_right);

    // Perform the event with a smaller distance
    if (d_stay < d_leak)  // Stay in current zone
    {
      // Tally mean intensity
      #pragma omp atomic
      J_nu_[p.ind][0] += p.e*dt_remaining*pc::c;

      // advect it
      double zone_vel[3], dvds;
      grid->get_velocity(p.ind,p.x,p.D,zone_vel, &dvds);
      p.x[0] += zone_vel[0]*dt_remaining;
      p.x[1] += zone_vel[1]*dt_remaining;
      p.x[2] += zone_vel[2]*dt_remaining;
      p.ind = grid->get_zone(p.x);

      // adiabatic loss (assumes small change in volume I think)
      p.e *= (1 - dvds*dt_remaining);

      // stay in this zone the remainder of the time step
      p.t += dt_remaining;
      dt_remaining = -1.0;
    }
    else // Leak to adjacent zone
    {
      double P_leak_left = sigma_leak_left / sigma_leak_tot;
      double xi2 = rangen.uniform();

      double zone_vel[3], dvds;
      grid->get_velocity(p.ind,p.x,p.D,zone_vel, &dvds);
      double dt_change = d_leak / pc::c;
      // Tally mean intensity
      #pragma omp atomic
      J_nu_[p.ind][0] += p.e*dt_change*pc::c;

      if (xi2 <= P_leak_left)  // leak left
      {
        p.ind--;
        double rr = p.r();
        p.x[0] += -p.x[0]/rr*dx + zone_vel[0]*dt_change;
        p.x[1] += -p.x[1]/rr*dx + zone_vel[1]*dt_change;
        p.x[2] += -p.x[2]/rr*dx + zone_vel[2]*dt_change;
      }
      else // leak right
      {
        p.ind++;
        double rr = p.r();
        p.x[0] += p.x[0]/rr*dx + zone_vel[0]*dt_change;
        p.x[1] += p.x[1]/rr*dx + zone_vel[1]*dt_change;
        p.x[2] += p.x[2]/rr*dx + zone_vel[2]*dt_change;
      }
      p.ind = grid->get_zone(p.x);
      if (p.ind == -1) {return absorbed;}
      if (p.ind == -2) {return escaped;}
      
      // adiabatic loss (assumes small change in volume I think)
      p.e *= (1 - dvds*dt_change);
      p.t += dt_change;



      dt_remaining -= dt_change;
    }

    // determine current zone and check for escape
    p.ind = grid->get_zone(p.x);
    if (p.ind == -1) {return absorbed;}
    if (p.ind == -2) {return escaped;}

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

  for (int i=0;i<nz;i++)
  {
    double dx;
    grid->get_zone_size(i,&dx);

    // determine if we will use ddmc in this zone
    double ztau = planck_mean_opacity_[i]*dx;
    if (ztau > ddmc_tau_) ddmc_use_in_zone_[i] = 1;
    else ddmc_use_in_zone_[i] = 0;

    // indices of adjacent zones
    int ip = i+1;
    if (ip == nz) ip = i;
    int im = i-1;
    if (im < 0)   im = 0;

    // diffusion probability in zone and adjacent zones
    double Dj0 = pc::c/(3.0*planck_mean_opacity_[i]);
    double Djp = pc::c/(3.0*planck_mean_opacity_[ip]);
    double Djm = pc::c/(3.0*planck_mean_opacity_[im]);

    double Dh_up = 2*dx*(Dj0*Djp)/(Dj0*dx + Djp*dx);
    double Dh_dn = 2*dx*(Dj0*Djm)/(Dj0*dx + Djm*dx);

    double rcoords[3];
    grid->coordinates(i,rcoords);
    double r_p = rcoords[0];
    grid->coordinates(im,rcoords);
    double r_m = rcoords[0];
    if (i == 0) r_m = 0;
    double r_0 = 0.5*(r_p + r_m);

    // diffusion probabilitys up and down
    // assumes spherical symmetry so the r^2/r_0^2 factor
    ddmc_P_up_[i] = (dt/dx)*(Dh_up/dx)*r_p*r_p/r_0/r_0;
    ddmc_P_dn_[i] = (dt/dx)*(Dh_dn/dx)*r_m*r_m/r_0/r_0;
    // boundary condition -- don't diffuse inward at innermost cell
    if (i==0) ddmc_P_dn_[i] = 0;

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
