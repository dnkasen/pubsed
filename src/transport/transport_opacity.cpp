#include <math.h>
#include <cassert>
#include <ctime>
#include "transport.h"
#include "physical_constants.h"
#include "radioactive.h"

using std::cout;
using std::cerr;
using std::endl;
namespace pc = physical_constants;

//-----------------------------------------------------------------
// set opacity and emissivities
//-----------------------------------------------------------------
void transport::set_opacity(double dt)
{

  double tend,tstr;
  double get_system_time(void);
  
  // tmp vector to hold emissivity
  vector<OpacityType> emis(nu_grid.size());
  vector<OpacityType> scat(nu_grid.size());
  emis.assign(emis.size(),0.0);

  // always do LTE on first step
  int nlte = gas_state_.use_nlte_;
  if (first_step_) gas_state_.use_nlte_ = 0;

  // zero out opacities, etc...
  for (int i=0;i<grid->n_zones;i++)
  {
    compton_opac[i]  = 0;
    photoion_opac[i] = 0;
    rosseland_mean_opacity_[i] = 0;
    planck_mean_opacity_[i]    = 0;
    emissivity_[i].wipe();
    for (int j=0;j<nu_grid.size();j++)
    {
      abs_opacity_[i][j] = 0;
      if (!omit_scattering_) scat_opacity_[i][j] = 0;
    }
  }

  radioactive radio;
  vector<double> X_now(grid->n_elems);

  // loop over my zones to calculate

  tstr = get_system_time();
  for (int i=my_zone_start_;i<my_zone_stop_;i++)
  {
    // pointer to current zone for easy access
    zone* z = &(grid->z[i]);

    //------------------------------------------------------
    // optical photon opacities
    //------------------------------------------------------

    // set up the state of the gas in this zone
    gas_state_.dens_ = z->rho;
    gas_state_.temp_ = z->T_gas;
    gas_state_.time_ = t_now_;
    if (gas_state_.temp_ < temp_min_value_) gas_state_.temp_ = temp_min_value_;
    if (gas_state_.temp_ > temp_max_value_) gas_state_.temp_ = temp_max_value_;

    // radioactive decay the composition
    for (size_t j=0;j<X_now.size();j++) X_now[j] = z->X_gas[j];
    if (!omit_composition_decay_) {
      radio.decay_composition(grid->elems_Z,grid->elems_A,X_now,t_now_);
    }
    gas_state_.set_mass_fractions(X_now);

    gas_state_.total_grey_opacity_ = gas_state_.smooth_grey_opacity_ + z->grey_opacity;

    if (first_step_)
      {
	  zone* z = &(grid->z[i]);
	  gas_state_.dens_ = z->rho;
	  gas_state_.temp_ = z->T_gas;

	  int solve_error = 0;
  
	  if ( (gas_state_.smooth_grey_opacity_ == 0) && (gas_state_.use_zone_dependent_grey_opacity_ == 0) )
	    {
	      solve_error = gas_state_.solve_state(J_nu_[i]);
	    }
      }
    

    else
      {

	if (radiative_eq)
	  {
	    solve_eq_temperature(i); // NLTE solution will take place in here
	  }
	else
	  {
	    int solve_error = 0;
	    if ( (gas_state_.smooth_grey_opacity_ == 0) && (gas_state_.use_zone_dependent_grey_opacity_ == 0) ){
	      solve_error = gas_state_.solve_state(J_nu_[i]);
	    }

	  }
      }

    
    // solve for the state
    // if (!gas_state_.grey_opacity_) solve_error = gas_state_.solve_state(J_nu_[i]);

  
    //gas_state_.print();

    if(write_levels) gas_state_.write_levels(i);

    // calculate the opacities/emissivities
    gas_state_.computeOpacity(abs_opacity_[i],scat,emis);

    double max_extinction = maximum_opacity_* z->rho;

    // save and normalize emissivity cdf
    grid->z[i].L_thermal = 0;
    if (nu_grid.size() == 1)
    {
      double bb_int = pc::sb*pow(grid->z[i].T_gas,4)/pc::pi;
      grid->z[i].L_thermal += 4*pc::pi*abs_opacity_[i][0]*bb_int;
      emissivity_[i].set_value(0,1);
      if (!omit_scattering_) scat_opacity_[i][0] = scat[0];
    }
    else for (int j=0;j<nu_grid.size();j++)
    {
      double ednu = emis[j]*nu_grid.delta(j);
      emissivity_[i].set_value(j,ednu);
      grid->z[i].L_thermal += 4*pc::pi * ednu;
      if (!omit_scattering_) scat_opacity_[i][j] = scat[j];

      // check for maximum opacity
      if (!omit_scattering_)
      {
        if (scat_opacity_[i][j] > max_extinction)
          scat_opacity_[i][j]   = max_extinction;
      }
      if (abs_opacity_[i][j]  > max_extinction)
        abs_opacity_[i][j]    = max_extinction;
    }
    emissivity_[i].normalize();

    // calculate mean opacities
    planck_mean_opacity_[i] =
      gas_state_.get_planck_mean(abs_opacity_[i],scat_opacity_[i]);
    rosseland_mean_opacity_[i] =
      gas_state_.get_rosseland_mean(abs_opacity_[i],scat_opacity_[i]);

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

  tend = get_system_time();
  if (verbose) cout << "# Calculated Radiative Equilib (" << (tend-tstr) << " secs) \n";
  // mpi reduce the results
  reduce_Tgas();


  //------------------------------------------------------------
  // Calcuate eps_imc...
  //------------------------------------------------------------
  for (int i=0;i<grid->n_zones;i++)
  {
    if (radiative_eq)
    {
      grid->z[i].eps_imc = 1.;
    }
    else
    {
       // Not distinguishing between lab frame density and comoving frame density
      double fleck_beta  = 4.0*pc::a*pow(grid->z[i].T_gas,4)/(grid->z[i].e_gas*grid->z[i].rho);
       // here planck mean opac has units cm^-1 .
       // When grey opacity is used, planck_mean_opacity should just be the correct grey opacity
       double tfac = pc::c*planck_mean_opacity_[i]*dt;
       double f_imc = fleck_alpha_*fleck_beta*tfac;

      // make sure to avoid divide by zero if fleck alpha is zero and not computing e_gas through hydro
      if (fleck_alpha_ == 0)
	f_imc = 0.;
       
       grid->z[i].eps_imc = 1.0/(1.0 + f_imc);
    }
  }

  // turn nlte back on after first step, if wanted
  if (first_step_) {
    gas_state_.use_nlte_ = nlte;
    first_step_ = 0;    }

  /*  
  // flag any error
  if (verbose)
  {
    if (solve_error == 1) std::cerr << "# Warning: root not bracketed in n_e solve\n";
    if (solve_error == 2) std::cerr << "# Warning: max iterations hit in n_e solve\n";
  }
  */
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
    i_nu = nu_grid.locate_within_bounds(nu);
    double a_opac = nu_grid.value_at(nu,abs_opacity_[p.ind],i_nu);
    double s_opac = 0;
    if (!omit_scattering_) s_opac = nu_grid.value_at(nu,scat_opacity_[p.ind],i_nu);
    opac = a_opac + s_opac;
    if (opac == 0) eps = 0;
    else eps  = a_opac/opac;
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

