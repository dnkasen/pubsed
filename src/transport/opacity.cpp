#include <math.h>
#include <cassert>
#include "transport.h"
#include "physical_constants.h"

namespace pc = physical_constants;

//-----------------------------------------------------------------
// set opacity and emissivities
//-----------------------------------------------------------------
void transport::set_opacity()
{
  // tmp vector to hold emissivity
  vector<double> emis(nu_grid.size());
  
  // tmp vector to hold line stuff
  vector<double> lopac(n_lines_);

  // loop over all zones
  for (int i=0;i<grid->n_zones;i++)
  {
    // pointer to current zone for easy access
    zone* z = &(grid->z[i]);
      
    //------------------------------------------------------
    // optical photon opacities
    //------------------------------------------------------

    // set up the state of the gas in this zone
    gas.dens = z->rho;
    gas.temp = z->T_gas;
    gas.time = t_now_;
    gas.set_mass_fractions(z->X_gas);
    
    // solve for the state (assume LTE now)
    if (!gas.grey_opacity_) gas.solve_state(1);

    // calculate the opacities/emissivities
    gas.computeOpacity(abs_opacity_[i],scat_opacity_[i],emis);
    
    // save and normalize emissivity cdf
    for (int j=0;j<nu_grid.size();j++)
      emissivity_[i].set_value(j,emis[j]*nu_grid.delta(j));
    emissivity_[i].normalize();
  
    // set line opacities
    if (use_detailed_lines_)
      gas.get_line_opacities(line_opacity_[i]);
    
    //------------------------------------------------------
    // gamma-ray opacity (compton + photo-electric)
    //------------------------------------------------------
    compton_opac[i]  = 0;
    photoion_opac[i] = 0;
    for (int k=0;k<grid->n_elems;k++)
    {
      double dens  = z->X_gas[k]*z->rho;
      double ndens = dens/(pc::m_p*grid->elems_A[k]);
      // compton scattering opacity
      compton_opac[i] += ndens*pc::thomson_cs*grid->elems_Z[k];
      // photoelectric opacity
      double photo = pow(pc::alpha_fs,4.0)*4.0*sqrt(2.0);
      photo *= pow(1.0*grid->elems_Z[k],5.0);
      photo *= pow(pc::m_e_MeV,3.5);
      photoion_opac[i] += ndens*2.0*pc::thomson_cs*photo;
    }

  }
}


//-----------------------------------------------------------------
// get comoving opacity at the frequency
// returns the frequency index of the photon in the
// comoving frame
//-----------------------------------------------------------------
int transport::get_opacity(particle &p, double dshift, double &opac, double &eps)
{
  assert(p.ind >= 0);


  // comoving frame frequency
  double nu = p.nu*dshift;
  int i_nu = 0;
  
  // get opacity if it is an optical photon. 
  if (p.type == photon)
  {
    // interpolate opacity at the local comving frame frequency
    i_nu = nu_grid.locate(nu);
    double a_opac = nu_grid.value_at(nu,abs_opacity_[p.ind],i_nu);
    double s_opac = nu_grid.value_at(nu,scat_opacity_[p.ind],i_nu);
    opac = a_opac + s_opac;
    eps  = a_opac/opac;
  }

  // get opacity if it is a gamma-ray
  if (p.type == gammaray)
  {
    double c_opac = compton_opac[p.ind]*klein_nishina(p.nu);
    double p_opac = photoion_opac[p.ind]*pow(p.nu,-3.5);
    opac = c_opac + p_opac;
    eps  = p_opac/(c_opac + p_opac);
  }

  return i_nu;
}


//-----------------------------------------------------------------
// Klein_Nishina correction to the Compton cross-section
// assumes energy x is in MeV
//-----------------------------------------------------------------
double transport::klein_nishina(double x)
{
  // divide by m_e c^2 = 0.511 MeV
  x = x/pc::m_e_MeV;
  double logfac = log(1 + 2*x);
  double term1 = (1+x)/x/x/x*(2*x*(1+x)/(1+2*x) - logfac);
  double term2 = 1.0/2.0/x*logfac;
  double term3 = -1.0*(1 + 3*x)/(1+2*x)/(1+2*x);
  double KN    = .75*(term1 + term2 + term3);
  return KN;
}


//-----------------------------------------------------------------
// calculate planck function in frequency units
//-----------------------------------------------------------------
double transport::blackbody_nu(double T, double nu)
{
  double zeta = pc::h*nu/pc::k/T;
  return 2.0*nu*nu*nu*pc::h/pc::c/pc::c/(exp(zeta)-1);
}





