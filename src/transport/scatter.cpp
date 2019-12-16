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

  // Update position of last interaction
  p->x_interact[0] = p->x[0];
  p->x_interact[1] = p->x[1];
  p->x_interact[2] = p->x[2];

  // do photon interaction physics
  if (p->type == photon)
  {
    // see if scattered
    if (rangen.uniform() > eps)
    {
      if (compton_scatter_photons_)
        compton_scatter_photon(p);
      else
        isotropic_scatter(p,0);
    }
    else
    {
      // check for effective scattering
      double z2 = rangen.uniform();
      // enforced radiative equilibrium always effective scatters
      if ((z2 > zone->eps_imc)||(radiative_eq))
        isotropic_scatter(p,1);
      else fate = absorbed;
    }
  }

  // gamma ray interaction physics
  if (p->type == gammaray)
  {
    // see if scattered
    if (rangen.uniform() > eps) compton_scatter(p);
    // or if absorbed, turn it into a photon
    else
    {
      grid->z[p->ind].L_radio_dep += p->e;
      p->type = photon;
      // isotropic emission in comoving frame
      double mu  = 1 - 2.0*rangen.uniform();
      double phi = 2.0*pc::pi*rangen.uniform();
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
// physics of compton scattering for gamma-rays
//------------------------------------------------------------
void transport::compton_scatter(particle *p)
{
  assert(p->ind >= 0);

  transform_lab_to_comoving(p);

  // sample new direction by rejection method
  double E_ratio;
  double D_new[3];
  while (true)
  {
    // isotropic new direction
    double mu  = 1 - 2.0*rangen.uniform();
    double phi = 2.0*pc::pi*rangen.uniform();
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
    if (rangen.uniform() < diff_cs) break;
  }

  // new frequency
  p->nu = p->nu*E_ratio;

  // add in gamma-ray energy deposition
  //if (p->type == gammaray) grid->z[p->ind].L_radio_dep += p->e*(1 - E_ratio);

  // sample whether we stay alive, if not become a photon
  if (rangen.uniform() > E_ratio)
  {
    grid->z[p->ind].L_radio_dep += p->e;
    p->type = photon;
    // isotropic emission in comoving frame
    double mu  = 1 - 2.0*rangen.uniform();
    double phi = 2.0*pc::pi*rangen.uniform();
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


void transport::sample_MB_vector(double T, double* v_e, double* p_d)
{

  while (true)
  {

    // if you prefer, you could also rejection sample to get v_tot.

    double v_tot = sqrt(2. * pc::k * T /pc::m_e) * mb_dv * (mb_cdf_.sample(rangen.uniform()) + rangen.uniform() );

    double mu  = 1. - 2.0*rangen.uniform();
    double phi = 2.0*pc::pi*rangen.uniform();
    double smu = sqrt(1 - mu*mu);
    double ed0 = smu*cos(phi);
    double ed1 = smu*sin(phi);
    double ed2 = mu;


    double omega = ed0 * p_d[0] + ed1 * p_d[1] + ed2 * p_d[2];

    // could tighten this bound if you think you know how small v_tot/C will be
    if (rangen.uniform() < 0.5 * (1. - omega * v_tot/pc::c)) // this is crucial. For the more relativistic case, the formula gets more complicated. See the discussion at the top of pdf page 135 (journal page 323) of the Pozdnyakov 1983 paper, which references a formula for sigma-hat four pages earlier
    {

      v_e[0] = v_tot * ed0;
      v_e[1] = v_tot * ed1;
      v_e[2] = v_tot * ed2;
      return;
    }
  }
}



//------------------------------------------------------------
// physics of compton scattering for photons of arbitrary energies
// samples from thermal velocity distribution of non-relativistic electrons
//------------------------------------------------------------
void transport::compton_scatter_photon(particle *p)
{

  assert(p->ind >= 0);

  double E_initial = p->e;

  transform_lab_to_comoving(p);

  zone *zone = &(grid->z[p->ind]);

  // Find random thermal electron velocity
  double v_sc[3];
  sample_MB_vector(zone->T_gas,v_sc,p->D);

  //Transform into rest frame of electon
  double v_tot = sqrt(v_sc[0] * v_sc[0] + v_sc[1] * v_sc[1] + v_sc[2] * v_sc[2]);
  double beta   = v_tot/pc::c;
  double gamma  = 1.0/sqrt(1 - beta*beta);
  double vdd = (v_sc[0] * p->D[0] + v_sc[1] * p->D[1] + v_sc[2] * p->D[2] );

  double dshift_into_scatterer = gamma * (1. - vdd/pc::c); // this is simplified since we defined beta and gamma based on the parallel velocity component

  // transform the 0th component (energy and frequency)
  p->e  *= dshift_into_scatterer;
  p->nu *= dshift_into_scatterer;

  // transform the 1-3 components (direction)
  // See Mihalas & Mihalas eq 89.8
  p->D[0] = 1.0/dshift_into_scatterer * (p->D[0] - gamma*v_sc[0]/pc::c * (1. - gamma*vdd/pc::c/(gamma+1)) );
  p->D[1] = 1.0/dshift_into_scatterer * (p->D[1] - gamma*v_sc[1]/pc::c * (1. - gamma*vdd/pc::c/(gamma+1)) );
  p->D[2] = 1.0/dshift_into_scatterer * (p->D[2] - gamma*v_sc[2]/pc::c * (1. - gamma*vdd/pc::c/(gamma+1)) );


  // sample new direction by rejection method
  double E_ratio;
  double D_new[3];

  while (true)
  {
    // isotropic new direction
    double mu  = 1. - 2.0*rangen.uniform();
    double phi = 2.0*pc::pi*rangen.uniform();
    double smu = sqrt(1. - mu*mu);
    D_new[0] = smu*cos(phi);
    D_new[1] = smu*sin(phi);
    D_new[2] = mu;

    // angle between old and new directions
    double cost = p->D[0]*D_new[0] + p->D[1]*D_new[1] + p->D[2]*D_new[2];
    // new energy ratio (E_new/E_old) at this angle (assuming lambda in MeV)
    E_ratio = 1./(1. + pc::h * p->nu / (pc::m_e_MeV * pc:: Mev_to_ergs)*(1. - cost));
    // klein-nishina differential cross-section
    double diff_cs = 0.5*(E_ratio*E_ratio*(1./E_ratio + E_ratio - 1. + cost*cost));
    // see if this scatter angle OK
    if (rangen.uniform() < diff_cs) break;
  }

  // new frequency
  p->nu = p->nu * E_ratio;
  p->e = p->e * E_ratio;


  //Transform out of rest frame of electon into comoving frame
  // shift back to comoving frame
  for (int i=0;i<3;i++) v_sc[i] = -1.*v_sc[i];
  vdd = (v_sc[0] * D_new[0] + v_sc[1] * D_new[1] + v_sc[2] * D_new[2] );

  double dshift_out_scatterer = gamma * (1. - vdd/pc::c);

  // transform the 0th component (energy and frequency)
  p->e  *= dshift_out_scatterer;
  p->nu *= dshift_out_scatterer;

  // transform the 1-3 components (direction)
  // See Mihalas & Mihalas eq 89.8
  p->D[0] = 1.0/dshift_out_scatterer * (D_new[0] - gamma*v_sc[0]/pc::c * (1. - gamma*vdd/pc::c/(gamma+1)) );
  p->D[1] = 1.0/dshift_out_scatterer * (D_new[1] - gamma*v_sc[1]/pc::c * (1. - gamma*vdd/pc::c/(gamma+1)) );
  p->D[2] = 1.0/dshift_out_scatterer * (D_new[2] - gamma*v_sc[2]/pc::c * (1. - gamma*vdd/pc::c/(gamma+1)) );

  // lorentz transform back to lab frame
  transform_comoving_to_lab(p);

  double E_final = p->e;

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
//   double mu  = 1 - 2.0*rangen.uniform();
//   double phi = 2.0*pc::pi*rangen.uniform();
//   double smu = sqrt(1 - mu*mu);
//   p->D[0] = smu*cos(phi);
//   p->D[1] = smu*sin(phi);
//   p->D[2] = mu;

//   // change wavelength, if needed
//   if (redistribute)
//   {
//     // sample frequency from local emissivity
//     int ilam  = emis[p->ind].sample(rangen.uniform());
//     p->nu = nu_grid.sample(ilam,rangen.uniform());
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
  double mu  = 1 - 2.0*rangen.uniform();
  double phi = 2.0*pc::pi*rangen.uniform();
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
