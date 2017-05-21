#include <math.h>
#include <gsl/gsl_rng.h>
#include <cassert>
#include "transport.h"
#include "physical_constants.h"

namespace pc = physical_constants;

//------------------------------------------------------------
// interaction physics
//------------------------------------------------------------
ParticleFate transport::do_scatter(particle *p, double eps)
{
  zone *zone = &(grid->z[p->ind]);
  ParticleFate fate = moving;

  // do photon interaction physics
  if (p->type == photon)
  {
    // debug
    if (gsl_rng_uniform(rangen) > eps) isotropic_scatter(p,0);
    else isotropic_scatter(p,1);
    return moving;

    // see if scattered 
    if (gsl_rng_uniform(rangen) > eps) isotropic_scatter(p,0);
    else
    {
      // check for effective scattering
      double z2;
      if (radiative_eq) z2 = 2;
      else z2 = gsl_rng_uniform(rangen);
      // do an effective scatter
      if (z2 > zone->eps_imc) isotropic_scatter(p,1);
      // otherwise really absorb (kill) it
      else  fate = absorbed; 
    }
  }
    
  // gamma ray interaction physics
  if (p->type == gammaray)
  {
    // see if scattered 
    if (gsl_rng_uniform(rangen) > eps) compton_scatter(p);
    // or if absorbed, turn it into a photon
    else 
    {
      grid->z[p->ind].L_radio_dep += p->e;
      p->type = photon;
      // isotropic emission in comoving frame
      double mu  = 1 - 2.0*gsl_rng_uniform(rangen);
      double phi = 2.0*pc::pi*gsl_rng_uniform(rangen);
      double smu = sqrt(1 - mu*mu);
      p->D[0] = smu*cos(phi);
      p->D[1] = smu*sin(phi);
      p->D[2] = mu;
      sample_photon_frequency(p);

      // lorentz transform back to lab frame
      transform_comoving_to_lab(p);
      // debug - didn't put in transform in
    } 
  }
  return fate;

}
      

//------------------------------------------------------------
// physics of compton scattering
//------------------------------------------------------------
void transport::compton_scatter(particle *p)
{
  assert(p->ind >= 0);

  // get doppler shift from lab to comoving frame
  double dshift_in = dshift_lab_to_comoving(p);
  
  // transform energy, frequency into comoving frame
  p->e  *= dshift_in;
  p->nu *= dshift_in;
  
  // sample new direction by rejection method
  double E_ratio;
  double D_new[3];
  while (true)
  {
    // isotropic new direction
    double mu  = 1 - 2.0*gsl_rng_uniform(rangen);
    double phi = 2.0*pc::pi*gsl_rng_uniform(rangen);
    double smu = sqrt(1 - mu*mu);
    D_new[0] = smu*cos(phi);
    D_new[1] = smu*sin(phi);
    D_new[2] = mu;
    
    // angle between old and new directions
    double cost = p->D[0]*D_new[0] + p->D[1]*D_new[1] + p->D[2]*D_new[2];
    // new energy ratio (E_new/E_old) at this angle (assuming lambda in MeV)
    E_ratio = 1/(1 + p->nu/pc::m_e_MeV*(1 - cost));
    // klein-nishina differential cross-section
    double diff_cs = 0.5*(E_ratio*E_ratio*(1/E_ratio + E_ratio - 1 + cost*cost));
    // see if this scatter angle OK
    if (gsl_rng_uniform(rangen) < diff_cs) break;
  }
  
  // new frequency
  p->nu = p->nu*E_ratio;
  
  // add in gamma-ray energy deposition
  //if (p->type == gammaray) grid->z[p->ind].L_radio_dep += p->e*(1 - E_ratio);

  // sample whether we stay alive, if not become a photon
  if (gsl_rng_uniform(rangen) > E_ratio) 
  {
    grid->z[p->ind].L_radio_dep += p->e;
    p->type = photon;
    // isotropic emission in comoving frame
    double mu  = 1 - 2.0*gsl_rng_uniform(rangen);
    double phi = 2.0*pc::pi*gsl_rng_uniform(rangen);
    double smu = sqrt(1 - mu*mu);
    D_new[0] = smu*cos(phi);
    D_new[1] = smu*sin(phi);
    D_new[2] = mu;
    sample_photon_frequency(p);
  }

  // set new direction
  p->D[0] = D_new[0];
  p->D[1] = D_new[1];
  p->D[2] = D_new[2];

  // lorentz transform back to lab frame
  transform_comoving_to_lab(p);
}


//------------------------------------------------------------
// physics of isotropic scattering
//------------------------------------------------------------
// void transport::isotropic_scatter(particle &p, int redistribute)
// {
//   assert(p->ind >= 0);
//   // get doppler shift from lab to comoving frame
//   double dshift_in = dshift_lab_to_comoving(&p);
  
//   // transform energy, frequency into comoving frame
//   p->e  *= dshift_in; 
//   p->nu *= dshift_in;

//   // Randomly generate new direction isotropically in comoving frame
//   double mu  = 1 - 2.0*gsl_rng_uniform(rangen);
//   double phi = 2.0*pc::pi*gsl_rng_uniform(rangen);
//   double smu = sqrt(1 - mu*mu);
//   p->D[0] = smu*cos(phi);
//   p->D[1] = smu*sin(phi);
//   p->D[2] = mu;

//   // change wavelength, if needed
//   if (redistribute)
//   {
//     // sample frequency from local emissivity
//     int ilam  = emis[p->ind].sample(gsl_rng_uniform(rangen));
//     p->nu = nu_grid.sample(ilam,gsl_rng_uniform(rangen));
//   }
  
//   // lorentz transform back to lab frame
//   transform_comoving_to_lab(&p);

// }





//------------------------------------------------------------
// FAKE physics of non-isotropic Compton scattering
//------------------------------------------------------------
void transport::isotropic_scatter(particle *p, int redist)
{
  double V[3], dvds;
  grid->get_velocity(p->ind,p->x,p->D,V,&dvds);

  // local velocity vector
  double vdotD  = V[0]*p->D[0] + V[1]*p->D[1] + V[2]*p->D[2];
  double beta2 = (V[0]*V[0] + V[1]*V[1] + V[2]*V[2])/pc::c/pc::c;
  double gamma = 1.0/sqrt(1 - beta2);
  double dshift_in  = gamma*(1 - vdotD/pc::c);

  // transform quantities into comoving frame
  p->e  = p->e*dshift_in;
  p->nu = p->nu*dshift_in;
  
  // sample new direction by rejection method
  double D_new[3];

  // choose new isotropic direction in comoving frame
  double mu  = 1 - 2.0*gsl_rng_uniform(rangen);
  double phi = 2.0*pc::pi*gsl_rng_uniform(rangen);
  double smu = sqrt(1 - mu*mu);
  D_new[0] = smu*cos(phi);
  D_new[1] = smu*sin(phi);
  D_new[2] = mu;
  
  // choose new wavelength if redistributed
  if (redist) sample_photon_frequency(p);

  // outgoing velocity vector
  for (int i=0;i<3;i++) V[i] = -1*V[i];
  
  // doppler shifts outgoing
  double vdp = (D_new[0]*V[0] + D_new[1]*V[1] + D_new[2]*V[2]);
  double vd_out = gamma*(1 - vdp/pc::c);
  
  // transformation of direction vector back into lab frame
  p->D[0] = 1.0/vd_out*(D_new[0] - gamma*V[0]/pc::c*(1 - gamma*vdp/pc::c/(gamma+1)));
  p->D[1] = 1.0/vd_out*(D_new[1] - gamma*V[1]/pc::c*(1 - gamma*vdp/pc::c/(gamma+1)));
  p->D[2] = 1.0/vd_out*(D_new[2] - gamma*V[2]/pc::c*(1 - gamma*vdp/pc::c/(gamma+1)));
  double norm = sqrt(p->D[0]*p->D[0] + p->D[1]*p->D[1] + p->D[2]*p->D[2]);
  p->D[0] = p->D[0]/norm;
  p->D[1] = p->D[1]/norm;
  p->D[2] = p->D[2]/norm;
  
  // transformation of energy/wavelength into lab frame
  p->e  = p->e*vd_out;
  p->nu = p->nu*vd_out;
}
  



