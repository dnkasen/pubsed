// DIFFUSION_METHOD 1=IMD 2=DDMC 3=RW
#include <math.h>
#include <gsl/gsl_rng.h>
#include <cassert>
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
  enum ParticleEvent {scatter, boundary, tstep};
  ParticleEvent event;

  ParticleFate fate = moving;

  // Constant for interfacing MC and DDMC
  double lambda_ddmc = 0.7104;
  double ddmc_sml_push = 1.0e-8;

  // pointer to current zone
  zone *zone = &(grid->z[p.ind]);
  int nz = grid->n_zones;

  // initialize particle's timestamp
  double dt_remaining = tstop - p.t;
  double dt_hydro = tstop - t_now_;

  while (fate == moving)
  {
    // indices of current and adjacent zones
    int ii = p.ind;
    int ip = ii + 1;
    int im = ii - 1;

    // Boundary zones
    if (ip == nz) ip = ii;
    if (im < 0)   im = 0;

    double dx;
    grid->get_zone_size(ii,&dx);

    double dxp1, dxm1;
    grid->get_zone_size(ip,&dxp1);
    grid->get_zone_size(im,&dxm1);

    // Use comoving nu to get the total/transport opacity, (abs + scat).
    // The Planck mean could be used for tallying energy absorption,
    // but use actual nu-dependent opacity for leakage opacity calculation.
    double k_p = planck_mean_opacity_[ii];
    double sigma_i, sigma_im1, sigma_ip1;

    int i_nu;
    double dshift, eps_i_cmf, eps_ip1_cmf, eps_im1_cmf;
    dshift = dshift_lab_to_comoving(&p);
    i_nu = get_opacity(p,dshift,sigma_i,eps_i_cmf);
    // Sorry about this super hacky way of getting neighbor's opacities
    p.ind = ip;
    i_nu = get_opacity(p,dshift,sigma_ip1,eps_ip1_cmf);
    p.ind = im;
    i_nu = get_opacity(p,dshift,sigma_im1,eps_im1_cmf);
    p.ind = ii; // set it back to original index

    // return if no long a DDMC particle
    double ztau = sigma_i * dx;
    if (ztau <= ddmc_tau_) {return moving;}

    // Getting radii at zone boundaries and center
    double rcoords[3];
    grid->coordinates(ii,rcoords);
    double r_p = rcoords[0]; // outer edge of zone ii
    grid->coordinates(im,rcoords);
    double r_m = rcoords[0]; // inner edge of zone ii
    if (ii == 0) {grid->get_r_out_min(&r_m);}
    double r_0 = 0.5*(r_p + r_m); // zone center

    // Compute left/right leakage opacity, including the
    // geometric factors for spherical coordinates
    double Xi, geo_factor;
    Xi = 1.0 + dx*dx/(12.0*r_0*r_0);
    geo_factor = r_m*r_m/r_0/r_0/Xi;
    double tau_im1 = sigma_im1*dxm1;
    int ddmc_on_im1 = 1;
    if (tau_im1 <= ddmc_tau_) ddmc_on_im1 = 0;
    double sigma_leak_left;
    if (ddmc_on_im1) sigma_leak_left  = (2.0/3.0/dx) * (1.0 / (sigma_i*dx + sigma_im1*dxm1)) * geo_factor;
    else sigma_leak_left  = (2.0/3.0/dx) * (1.0 / (sigma_i*dx + 2.0*lambda_ddmc)) * geo_factor;

    geo_factor = r_p*r_p/r_0/r_0/Xi;
    double tau_ip1 = sigma_ip1*dxp1;
    int ddmc_on_ip1 = 1;
    if (tau_ip1 <= ddmc_tau_) ddmc_on_ip1 = 0;
    double sigma_leak_right;
    if (ddmc_on_ip1) sigma_leak_right = (2.0/3.0/dx) * (1.0 / (sigma_i*dx + sigma_ip1*dxp1)) * geo_factor;
    else sigma_leak_right = (2.0/3.0/dx) * (1.0 / (sigma_i*dx + 2.0*lambda_ddmc)) * geo_factor;

    // Leakage away from outer boundary of problem domain
    if (ip == ii)
      {sigma_leak_right = (2.0/3.0/dx) * (1.0 / (sigma_i*dx + 2.0*lambda_ddmc)) * geo_factor;}
    // no left leakage in the innermost zone
    if (ii == 0) sigma_leak_left = 0.0;

    double sigma_leak_tot = sigma_leak_left + sigma_leak_right;

    double xi = rangen.uniform();
    double d_stay, d_leak;
    bool leaked2imc = false;

    // Distance until end of time step
    d_stay = pc::c * dt_remaining;

    // Distance to leakage event
    d_leak = -log(xi) / (sigma_leak_left + sigma_leak_right);

    // In non-grey RT, we also need to account for physical and
    // effective scattering that could change the particle's frequency.
    double d_sc, k_es_inelastic;

    if (nu_grid_.size() == 1) // no extra scattering for grey case
      {d_sc = std::numeric_limits<double>::infinity();}
    else  // non-grey case
    {
      double tau_r = -1.0*log(rangen.uniform());
      if (eps_i_cmf > 0.0)
      {
        k_es_inelastic = sigma_i*eps_i_cmf;
        // setting elastic_frac from emissivity_
        dshift = dshift_lab_to_comoving(&p);
        i_nu = get_opacity(p,dshift,sigma_i,eps_i_cmf);
        // emissivity_ has been normalized in transport_opacity.cpp
        double elastic_frac = emissivity_[p.ind].get_value(i_nu);
        double inelastic_frac = 1.0 - elastic_frac;
        k_es_inelastic *= inelastic_frac;
        d_sc = tau_r / k_es_inelastic;
        if (d_sc < 0.0) std::cout << "Negative d_sc!" << std::endl;
      }
      else
      {d_sc = std::numeric_limits<double>::infinity();}
    }
    // find out which event actually happens (one with shortest distance)
    double this_d, this_dt;
    if ((d_sc < d_leak) && (d_sc < d_stay))
      {event = scatter;    this_d = d_sc;}
    else if (d_leak < d_stay)
      {event = boundary;   this_d = d_leak;}
    else
      {event = tstep;      this_d = d_stay;}
    this_dt = this_d / pc::c;
    p.t += this_dt;
    dt_remaining -= this_dt;

    // tally the contribution of zone's radiation energy
    // only one factor of dshift above because opacity is in cmf,
    // just need to covert p.e from lab to cmf.
    #pragma omp atomic
    J_nu_[p.ind][0] += p.e*this_d;
    #pragma omp atomic
    grid->z[p.ind].e_abs  += (p.e*dshift)*this_d*sigma_i*eps_i_cmf;

    // Perform the event with a smaller distance
    if (event == scatter)  // effective scattering
    {
      // Ensure the effective scattering is inelastic
      int now_i_nu = i_nu;
      while (now_i_nu == i_nu)
      {
        fate = do_scatter(&p,1.0);
        dshift = dshift_lab_to_comoving(&p);
        now_i_nu = get_opacity(p,dshift,sigma_i,eps_i_cmf);
      }
    }
    else if (event == tstep)  // Stay in current zone
    {
      // stay in this zone the remainder of the time step
      dt_remaining = -1.0;
      fate = stopped;
    }
    else if (event == boundary) // Leak to adjacent zone
    {
      double P_leak_left = sigma_leak_left / sigma_leak_tot;
      double xi2 = rangen.uniform();

      // Step 1: leakage to the neighboring zone
      double dr;
      double mu, phi, smu;

      if (xi2 <= P_leak_left)  // leak left
      {
        p.ind--;
        double rr = p.r();
        dr = rr - r_m;
        if (ddmc_on_im1) dr += rangen.uniform()*dxm1;
        else dr *= (1.0 + ddmc_sml_push);
        //dr += rangen.uniform()*dxm1;

        p.x[0] -= p.x[0]/rr*dr;
        p.x[1] -= p.x[1]/rr*dr;
        p.x[2] -= p.x[2]/rr*dr;

        if (!ddmc_on_im1)
        {
          // need to place the particle at the interface
          rr = p.r();
          dr = r_m - rr;
          dr *= (1.0-ddmc_sml_push);
          leaked2imc = true;

          // Sample velocity from face of blackbody
          // in case of DDMC-to-IMC leakage
          transform_lab_to_comoving(&p);
          sample_dir_from_blackbody_surface(&p);

          mu = (p.x[0]*p.D[0] + p.x[1]*p.D[1] + p.x[2]*p.D[2]) / p.r();
          if (mu > 0.0)  // make sure it points toward -ve radial direction
            {p.D[0] = -p.D[0]; p.D[1] = -p.D[1]; p.D[2] = -p.D[2];}

          //return moving;
          transform_comoving_to_lab(&p);
        }
      }
      else // leak right
      {
        p.ind++;
        double rr = p.r();
        dr = r_p - rr;
        if (ddmc_on_ip1) dr += rangen.uniform()*dxp1;
        else dr *= (1.0 + ddmc_sml_push);
        //dr += rangen.uniform()*dxp1;

        p.x[0] += p.x[0]/rr*dr;
        p.x[1] += p.x[1]/rr*dr;
        p.x[2] += p.x[2]/rr*dr;

        if (!ddmc_on_ip1)
        {
          rr = p.r();
          dr = rr - r_p;
          dr *= (1.0-ddmc_sml_push);
          leaked2imc = true;

          // Sample velocity from face of blackbody
          // in case of DDMC-to-IMC leakage
          transform_lab_to_comoving(&p);
          sample_dir_from_blackbody_surface(&p);

          mu = (p.x[0]*p.D[0] + p.x[1]*p.D[1] + p.x[2]*p.D[2]) / p.r();
          if (mu < 0.0)  // make sure it points towards +ve radial direction
           {p.D[0] = -p.D[0]; p.D[1] = -p.D[1]; p.D[2] = -p.D[2];}
          //return moving;
          transform_comoving_to_lab(&p);
        }
      }
    }  // end (event == boundary)

    // Get zone velocity and velocity gradient
    double zone_vel[3], dvds;
    grid->get_velocity(p.ind,p.x,p.D,zone_vel, &dvds);

    // Compute the Dln(rho)/Dt term
    double vel = sqrt(zone_vel[0]*zone_vel[0] + zone_vel[1]*zone_vel[1] + zone_vel[2]*zone_vel[2]);
    double v_over_r = vel / p.r();
    double dlnrhodt = dvds + 2.0*v_over_r; // For homology it boils down to 3*dv/dr
    double doppler_shift_ddmc = (1.0 - dlnrhodt*this_dt/3.0);

    // Step 2: adiabatic loss in comoving frame
    // Apply adiabatic loss (assumes small change in Dln(rho)/Dt)
    transform_lab_to_comoving(&p);
    p.e *= doppler_shift_ddmc;
    p.nu *= doppler_shift_ddmc;
    transform_comoving_to_lab(&p);

    // Step 3: advection for particles left in DDMC zones
    if (!leaked2imc)
    {
      p.x[0] += zone_vel[0]*this_dt;
      p.x[1] += zone_vel[1]*this_dt;
      p.x[2] += zone_vel[2]*this_dt;
    }

    // determine current zone and check for escape
    p.ind = grid->get_zone(p.x);
    if (p.ind == -1) {fate = absorbed;}
    if (p.ind == -2) {fate = escaped;}

  }  // end while(fate == moving)

  return fate;
}


// ------------------------------------------------------
// Interfacing with DDMC zones, allowing IMC-to-DDMC conversion
// If the neighbor is in DDMC, there is a probability the particle
// gets converted into DDMC. If the particle is not converted,
// it is returned to the MC region.
// ------------------------------------------------------
int transport::move_across_DDMC_interface(particle &p, int new_ind, double sigma_i, double dr)
{
  // gather information for neighboring zone
  int ip = p.ind + 1;
  int im = p.ind - 1;

  int nz = grid->n_zones;
  if (ip == nz) ip = p.ind;
  if (im < 0)   im = 0;

  // Getting radii at zone boundaries and center
  double rcoords[3];
  grid->coordinates(p.ind,rcoords);
  double r_p = rcoords[0]; // outer edge of zone ii
  grid->coordinates(im,rcoords);
  double r_m = rcoords[0]; // inner edge of zone ii
  if (p.ind == 0) {grid->get_r_out_min(&r_m);} // Getting r_out.min

  double r_interface, r_0;
  if (new_ind == ip)
  {
    r_interface = r_p;
    r_0 = r_interface + 0.5*dr;
  }
  else if (new_ind == im)
  {
    r_interface = r_m;
    r_0 = r_interface - 0.5*dr;
  }
  else std::cerr << "transport.cpp: Unknown boundary crossing type!  "\
                   << new_ind << " " << im << " " << ip  << std::endl;

  // Compute IMC-to-DDMC conversion probability
  transform_lab_to_comoving(&p);
  double rr = p.r();
  // co-moving frame velocity vector
  double mu = (p.x[0]*p.D[0] + p.x[1]*p.D[1] + p.x[2]*p.D[2]) / rr;
  mu = fabs(mu); // get normal

  double k_p = planck_mean_opacity_[new_ind];
  double Xi, geo_factor;
  Xi = 1.0 + dr*dr/(12.0*r_0*r_0);
  geo_factor = r_interface*r_interface/r_0/r_0/Xi;
  double p_convert;
  // Asymptotic diffusion limit
  p_convert = 4.0 * (1.0 + 1.5*mu) / (3.0*sigma_i*dr + 6.0*0.7104);
  // Alternative formalism in Densmore, Evans, and Buksas (2008)
  //p_convert = 4.0 * (0.91 + 1.635*mu) / (3.0*sigma_i*dr + 6.0*0.7104);

  double xi = rangen.uniform();
  std::vector<double> rand;
  rand.push_back(rangen.uniform());
  rand.push_back(rangen.uniform());
  rand.push_back(rangen.uniform());
  double new_r[3];
  double ddmc_sml_push = 1.0e-8;

  if (xi <= p_convert) // converted to DDMC
  {
    p.ind = new_ind;

    grid->sample_in_zone(p.ind,rand,new_r);
    p.x[0] = new_r[0];
    p.x[1] = new_r[1];
    p.x[2] = new_r[2];

    return 1; // no longer transport with MC
  }
  else // returned to original MC zone
  {
    grid->sample_in_zone(p.ind,rand,new_r);
    p.x[0] = new_r[0];
    p.x[1] = new_r[1];
    p.x[2] = new_r[2];

    double dr;
    if (new_ind == ip)
    {
      dr = r_interface - p.r();
      dr *= (1.0-ddmc_sml_push);

      p.x[0] += p.x[0]/p.r()*dr;
      p.x[1] += p.x[1]/p.r()*dr;
      p.x[2] += p.x[2]/p.r()*dr;
    }
    else if (new_ind == im)
    {
      dr = p.r() - r_interface;
      dr *= (1.0-ddmc_sml_push);
      p.x[0] -= p.x[0]/p.r()*dr;
      p.x[1] -= p.x[1]/p.r()*dr;
      p.x[2] -= p.x[2]/p.r()*dr;
    }

    // emit from blackbody face in comoving frame
    sample_dir_from_blackbody_surface(&p);

    // make sure the particle moves away from the zone face
    double mu = (p.x[0]*p.D[0] + p.x[1]*p.D[1] + p.x[2]*p.D[2]) / p.r();
    if (new_ind == ip)  // returned from i+1 zone
    {
      if (mu > 0.0)
        {p.D[0] = -p.D[0]; p.D[1] = -p.D[1]; p.D[2] = -p.D[2];}
    }
    else if (new_ind == im)  // returned from i-1 zone
    {
      if (mu < 0.0)
        {p.D[0] = -p.D[0]; p.D[1] = -p.D[1]; p.D[2] = -p.D[2];}
    }
  }

  transform_comoving_to_lab(&p);
  return 0;
}


// ------------------------------------------------------
// Propagate a particle using the
// Random Walk Monte Carlo approach.
// Reference: Fleck & Canfield, J. of Comput. Physics 54, 508-523 (1984)
// ------------------------------------------------------
void transport::setup_RandomWalk(){
  //============//
  // RANDOMWALK //
  //============//
  // calculate randomwalk diffusion time probability array
  int sumN = params_->getScalar<int>("randomwalk_sumN");
  int npoints = params_->getScalar<int>("randomwalk_npoints");
  double randomwalk_max_x = params_->getScalar<double>("randomwalk_max_x");

  randomwalk_x.init(0, randomwalk_max_x, npoints);
  randomwalk_Pescape.resize(npoints);

  #pragma omp parallel for
  for(int i=1; i<=npoints; i++){
    double x = randomwalk_x.right(i);

    double sum = 0;
    for(int n=1; n<=sumN; n++){
      double tmp = exp(-x * (n*pc::pi)*(n*pc::pi));
      if(n%2 == 0) tmp *= -1;
      sum += tmp;
    }

    randomwalk_Pescape[i] = 1.0-2.*sum;
  }

  // normalize the results
  for(int i=0; i<npoints; i++)
    randomwalk_Pescape[i] /= randomwalk_Pescape[npoints-1];
}
void random_direction(double dir[3], thread_RNG& rangen){
  // double magnitude = 0;
  // for(int i=0; i<3; i++){
  //   dir[i] = rangen.uniform();
  //   magnitude += dir[i]*dir[i];
  // }
  // magnitude = sqrt(magnitude);
  // for(int i=0; i<3; i++) dir[i] /= magnitude;
  double costheta = 2.*rangen.uniform() - 1.;
  double sintheta = sqrt(1.-costheta*costheta);
  double phi = 2.*M_PI * rangen.uniform();
  dir[0] = sintheta * cos(phi);
  dir[1] = sintheta * sin(phi);
  dir[2] = costheta;
}
double sample_CDF(const vector<double>& CDF, const locate_array& x, const double u){
  int index = upper_bound(CDF.begin(), CDF.end(), u) - CDF.begin();
  assert(index<CDF.size());
  assert(index>=0);

  double x1 = x.left(index);
  double x2 = x.right(index);
  double P1 = index==0 ? 0 : CDF[index-1];
  double P2 = CDF[index];
  assert(u>=P1);
  assert(u<=P2);
  double result = x1 + (x2-x1)/(P2-P1) * (u-P1);
  assert(result<=x2);
  assert(result>=x1);
  return result;
}
double interpolate_CDF(const vector<double>& CDF, const locate_array& x, const double xval){
  int index = x.locate(xval);
  assert(index<CDF.size());
  assert(index>=0);

  double x1 = x.left(index);
  double x2 = x.right(index);
  assert(xval>=x1);
  assert(xval<=x2);
  double P1 = index==0 ? 0 : CDF[index-1];
  double P2 = CDF[index];
  double result = P1 + (P2-P1)/(x2-x1) * (xval-x1);
  assert(result>=P1);
  assert(result<=P2);
  return result;
}
ParticleFate transport::discrete_diffuse_RandomWalk(particle &p, double t_stop)
{
  int stop = 0;
  if (steady_state) t_stop = 1e99;

  double dx;

  // set this zone
  while (!stop)
  {
    double dt_remaining = t_stop - p.t;
    if (steady_state) dt_remaining = 1e99;
    assert(dt_remaining > 0);

    // find current zone and check for escape
    p.ind = grid->get_zone(p.x);
    grid->get_zone_size(p.ind,&dx);

    if (p.ind == -1) {return absorbed;}
    if (p.ind == -2) {return escaped;}

    // total probability of diffusing to the edge of the sphere
    double D = pc::c/(3.0*planck_mean_opacity_[p.ind]);// * 3./4.;
    double X = dt_remaining*D/(dx*dx);

    double dt_step, R_diffuse;
    double u = rangen.uniform();
    double sampled_X = sample_CDF(randomwalk_Pescape, randomwalk_x, u);
    if (sampled_X < X)
    { // particle reaches surface before census
      // use sampled X to determine in-flight time
      R_diffuse = dx;
      dt_step = sampled_X*dx*dx/D;
      stop = 0;
      //assert(p.t + dt_step <= t_stop + dt_step*1e-6);
    }
    else
    { // particle still in sphere at census
      // use sampled X to determine travel displacement
      dt_step = t_stop - p.t;
      R_diffuse = sqrt((dt_step*D) / sampled_X);
      assert(R_diffuse <= dx);
      stop = 1;
    }
    p.t += dt_step;
    if(p.t >= t_stop) stop = 1;
    assert(R_diffuse <= dx);

    // add in tally of absorbed and total radiation energy
    //#pragma omp atomic
    //zone->e_abs += p.e*ddmc_P_abs_[p.ind];
    //zone->e_rad += p.e*ddmc_P_stay_[p.ind];
    #pragma omp atomic
    J_nu_[p.ind][0] += p.e*dt_step*pc::c;
    #pragma omp atomic
    grid->z[p.ind].e_abs  += p.e*dt_step*pc::c*planck_mean_opacity_[p.ind];

    // move the particle a distance R_diffuse
    double diffuse_dir[3];
    random_direction(diffuse_dir, rangen);
    for(int i=0; i<3; i++) p.x[i] += diffuse_dir[i] * R_diffuse;

    // get outgoing direction
    do{
      random_direction(p.D, rangen);
    } while (p.D[0]*diffuse_dir[0] + p.D[1]*diffuse_dir[1] + p.D[2]*diffuse_dir[2] < 0);

    // advect it
    double zone_vel[3], dvds;
    grid->get_velocity(p.ind,p.x,p.D,zone_vel, &dvds);
    p.x[0] += zone_vel[0]*dt_step;
    p.x[1] += zone_vel[1]*dt_step;
    p.x[2] += zone_vel[2]*dt_step;

    // adiabatic loss (assumes small change in volume I think)
    p.e *= (1 - dvds*dt_step);
  }

  // find current zone and check for escape
  p.ind = grid->get_zone(p.x);

  if (p.ind == -1) {return absorbed;}
  if (p.ind == -2) {return escaped;}
  return stopped;
}


// ------------------------------------------------------
// Calculate the probabilities of diffusion
// for now this only works in 1D spherical coords
// ------------------------------------------------------
void transport::compute_diffusion_probabilities(double dt)
{
  int nz = grid->n_zones;

  double dtau_ddmc = 0.0, dtau_mc = 0.0;
  for (int i=0;i<nz;i++)
  {
    double dx;
    grid->get_zone_size(i,&dx);

    // determine if we will use ddmc in this zone
    double ztau = planck_mean_opacity_[i]*dx;
    if (ztau > ddmc_tau_) ddmc_use_in_zone_[i] = 1;
    else ddmc_use_in_zone_[i] = 0;

    // Accumulating tau in mc and ddmc regions
    if (ztau > ddmc_tau_) dtau_ddmc+= ztau;
    else dtau_mc += ztau;

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
  if (MPI_myID == 0)
  {
    std::cout << "# Integrated tau in DDMC regions: " << dtau_ddmc << std::endl;
    std::cout << "# Integrated tau in MC regions: " << dtau_mc << std::endl;
  }
}

// Sample particle's direction from the face of a blackbody,
// which has a normal component proportional to sqrt(rand()).
void transport::sample_dir_from_blackbody_surface(particle* p)
{
  double v1, v2, v3;
  double mu, phi, smu;
  mu = sqrt(rangen.uniform()); // using sqrt(rand)
  phi = 2.0*pc::pi*rangen.uniform();
  smu = sqrt(1.0 - mu*mu);

  v1 = smu*cos(phi);
  v2 = smu*sin(phi);
  v3 = mu;

  double r1, r2, r3, sintheta, costheta;
  // unit r_hat vector of the particle
  r1 = p->x[0] / p->r();
  r2 = p->x[1] / p->r();
  r3 = p->x[2] / p->r();
  costheta = r3;
  sintheta = sqrt(1 - costheta*costheta);

  double v1p, v2p, v3p;
  // rotate the direction vector to the particle's location
  // Rodrigues' rotation formula
  v1p = v1*costheta + r1*v3*sintheta - r2*(1.0-costheta)*(-r2*v1 + r1*v2);
  v2p = v2*costheta + r2*v3*sintheta - r1*(1.0-costheta)*(-r2*v1 + r1*v2);
  v3p = v3*costheta - (r1*v1+r2*v2)*sintheta;
  double norm = sqrt(v1p*v1p + v2p*v2p + v3p*v3p);
  v1p /= norm; // renormalize
  v2p /= norm;
  v3p /= norm;

  p->D[0] = v1p;
  p->D[1] = v2p;
  p->D[2] = v3p;
}
