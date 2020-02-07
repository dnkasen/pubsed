#include <stdio.h>
#include <math.h>
#include <iostream>
#include "transport.h"
#include "physical_constants.h"
#include "grid_general.h"

namespace pc = physical_constants;

//-------------------------------------------------------------
//  Solve for the temperature assuming radiative equilibrium
//-------------------------------------------------------------

int transport::solve_state_and_temperature(GasState* gas_state_ptr, int i)
{

  vector<OpacityType> emis(nu_grid_.size());
  vector<OpacityType> scat(nu_grid_.size());
  emis.assign(emis.size(),0.0);

  int solve_error = 0;

  // Simple option to set gas temp based on radiation energy density
  if (set_Tgas_to_Trad_ == 1)
  {
    grid->z[i].T_gas = pow(grid->z[i].e_rad/pc::a,0.25);
    return 0;
  }
  // Otherwise set gas temperature from balancing heating-cooling
  else
  {
    zone* z = &(grid->z[i]);
    gas_state_ptr->dens_ = z->rho;
    gas_state_ptr->temp_ = z->T_gas;

    // For LTE, do an initial solve of the gas state
    if (gas_state_ptr->use_nlte_ == 0)
    {
      solve_error = gas_state_ptr->solve_state();
      gas_state_ptr->computeOpacity(abs_opacity_[i],scat,emis);
    }

    // Calculate equilibrium temperature.
    // Additional gas_state solve may also happen here
    grid->z[i].T_gas = temp_brent_method(gas_state_ptr, i,1,solve_error);

    if (gas_state_ptr->use_nlte_ == 0)
    {
      // For LTE, do a final solve
      solve_error = gas_state_ptr->solve_state();
    }

    if (gas_state_ptr->use_nlte_)
    {
      bf_heating[i] = gas_state_ptr->bound_free_heating_rate(grid->z[i].T_gas,J_nu_[i]);
      ff_heating[i] = gas_state_ptr->free_free_heating_rate(grid->z[i].T_gas,J_nu_[i]);
      bf_cooling[i] = gas_state_ptr->bound_free_cooling_rate(grid->z[i].T_gas);
      ff_cooling[i] = gas_state_ptr->free_free_cooling_rate(grid->z[i].T_gas);
     coll_cooling[i] = gas_state_ptr->collisional_net_cooling_rate(grid->z[i].T_gas);
    }
    return solve_error;
  }

}

void transport::solve_eq_temperature()
{
  int solve_error = 0;
  GasState* gas_state_ptr = &(gas_state_vec_[0]);
#pragma omp parallel for default(none) firstprivate(gas_state_ptr, solve_error)
  for (int i=my_zone_start_;i<my_zone_stop_;i++)
  {
    if (set_Tgas_to_Trad_ == 1)
    grid->z[i].T_gas = pow(grid->z[i].e_rad/pc::a,0.25);
    else
  {
      // solve_error won't be updated here because that's for the gas_state solve which isn't happening here
     grid->z[i].T_gas = temp_brent_method(gas_state_ptr, i,0,solve_error);

	  if (gas_state_ptr->use_nlte_)
    {
	    bf_heating[i] = gas_state_ptr->bound_free_heating_rate(grid->z[i].T_gas,J_nu_[i]);
	    ff_heating[i] = gas_state_ptr->free_free_heating_rate(grid->z[i].T_gas,J_nu_[i]);
	    bf_cooling[i] = gas_state_ptr->bound_free_cooling_rate(grid->z[i].T_gas);
	    ff_cooling[i] = gas_state_ptr->free_free_cooling_rate(grid->z[i].T_gas);
	    coll_cooling[i] = gas_state_ptr->collisional_net_cooling_rate(grid->z[i].T_gas);
	  }

	}
  }
  reduce_Tgas();
}


//***************************************************************/
// This is the function that expresses radiative equillibrium
// in a cell (i.e. E_absorbed = E_emitted).  It is used in
// The Brent solver below to determine the temperature such
// that radiative equilibrium holds
//
// This LTE version assumes the cooling emission is just given
// by Kirchofs law: j_nu = k_abs*B_nu where k_abs is the absorptive
// opacity and B_nu the Planck function.  This is integrated
// over frequency to get the total cooling radiative
//
// The NLTE version below carries out a more detailed
// treatment of the cooling rate based on microphysics
//
// solve_flag = 0 use the old opacities kappa(T_old) to
//              compute emission rates
//            = 1 recalculate the opacities kappa(T_new) to
//              compute emission rates
//
//************************************************************/
double transport::rad_eq_function_LTE(GasState* gas_state_ptr, int c,double T, int solve_flag, int & solve_error)
{
  zone* z = &(grid->z[c]);
  gas_state_ptr->dens_ = z->rho;
  gas_state_ptr->temp_ = T;

  // helper variables need for call (will not be used)
  vector<OpacityType> emis(nu_grid_.size());
  vector<OpacityType> scat(nu_grid_.size());
  emis.assign(emis.size(),0.0);

  // recalculate opacities based on current T if desired
  if (solve_flag)
  {
    // solve_error = gas_state_ptr->solve_state();
    gas_state_ptr->computeOpacity(abs_opacity_[c],scat,emis);
  }

  // total energy emitted (to be calculated)
  double E_emitted = 0.;

  // total energy absorbed in zone
  double E_absorbed = 0.;
  // if not solve_flag, use the absorption rate estimated during transport
  if (solve_flag == 0)
    E_absorbed = grid->z[c].e_abs;

  // Calculate total emission assuming no frequency (grey) opacity
  if (nu_grid_.size() == 1)
  {
    E_emitted = 4.0*pc::pi*abs_opacity_[c][0]*pc::sb/pc::pi*pow(T,4);

    if (solve_flag && solve_error == 0)
      E_absorbed = pc::c *abs_opacity_[c][0] * grid->z[c].e_rad;
  }

  // integrate emisison over frequency (angle
  // integration gives the 4*PI) to get total
  // ergs/sec/cm^3 radition emitted. Opacities are
  // held constant for this (assumed not to change
  // much from the last time step).
  else for (int i=0;i<nu_grid_.size();i++)
  {
    double dnu  = nu_grid_.delta(i);
    double nu   = nu_grid_.center(i);
    double B_nu = blackbody_nu(T,nu);
    double kappa_abs = abs_opacity_[c][i];
    E_emitted += 4.0*pc::pi*kappa_abs*B_nu*dnu;
    if (solve_flag == 1)
      E_absorbed += 4.0*pc::pi*kappa_abs*J_nu_[c][i]*dnu;
  }

  // radiative equillibrium condition: "emission equals absorbtion"
  // return to Brent function to iterate this to zero
  return (E_emitted - E_absorbed);
}


//***************************************************************/
// This is the function that expresses radiative equillibrium
// in a cell (i.e. E_absorbed = E_emitted).  It is used in
// The Brent solver below to determine the temperature such
// that radiative equilibrium holds
//
// The NLTE version below carries out a more detailed
// treatment of the cooling rate than the NLTE one above
//
// solve_flag = 0 use the old opacities kappa(T_old) to
//              compute emission rates
//            = 1 recalculate the opacities kappa(T_new) to
//              compute emission rates
//
//************************************************************/
double transport::rad_eq_function_NLTE(GasState* gas_state_ptr, int c,double T, int solve_flag, int &solve_error)
{

  zone* z = &(grid->z[c]);
  gas_state_ptr->dens_ = z->rho;
  gas_state_ptr->temp_ = T;

  // make sure grey_opacity is not being used
  if ( (gas_state_ptr->smooth_grey_opacity_ == 1) || (gas_state_ptr->use_zone_dependent_grey_opacity_ == 1) )
  {
	std::cerr << "# ERROR: NLTE solve should not be used with grey opacity\n";
	exit(1);
  }

  // if flag set, recompute the entire NLTE problem for this iteration
  if (solve_flag)
	solve_error = gas_state_ptr->solve_state(J_nu_[c]);

  // total energy absorbed
  double E_absorbed = gas_state_ptr->free_free_heating_rate(T,J_nu_[c]) +
        gas_state_ptr->bound_free_heating_rate(T,J_nu_[c]) ;

  // total energy emitted
  double E_emitted =  E_emitted= gas_state_ptr->free_free_cooling_rate(T) +
        gas_state_ptr->bound_free_cooling_rate(T);

  if (gas_state_ptr->use_collisions_nlte_)
      E_emitted += gas_state_ptr->collisional_net_cooling_rate(T);

  // radiative equillibrium condition: "emission equals absorbtion"
  // return to Brent function to iterate this to zero
  return (E_emitted - E_absorbed);
}


//-----------------------------------------------------------
// Brents method (from Numerical Recipes) to solve
// non-linear equation for T in rad equillibrium
//-----------------------------------------------------------
// definitions used for temperature solver

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
double transport::temp_brent_method(GasState* gas_state_ptr, int cell, int solve_flag, int &solve_error)
{
  double brent_solve_tolerance = 1.0e-2;
  double temp_range_min = temp_min_value_;
  double temp_range_max = temp_max_value_;

  int ITMAX = 100;
  double EPS = 3.0e-8;
  int iter;

  // Initial guesses
  double a=temp_range_min;
  double b=temp_range_max;
  double c=b;
  double d,e,min1,min2;
  double fa,fb = 0.;
  if (gas_state_ptr->use_nlte_ == 0)
    {
      fa=rad_eq_function_LTE(gas_state_ptr, cell,a,solve_flag,solve_error);
      fb=rad_eq_function_LTE(gas_state_ptr, cell,b,solve_flag,solve_error);
    }
  else
    {
      fa=rad_eq_function_NLTE(gas_state_ptr, cell,a,solve_flag,solve_error);
      fb=rad_eq_function_NLTE(gas_state_ptr, cell,b,solve_flag,solve_error);
    }

  double fc,p,q,r,s,tol1,xm;

  //if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
  //  printf("Root must be bracketed in zbrent");
  fc=fb;
  for (iter=1;iter<=ITMAX;iter++) {
    if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
      c=a;
      fc=fa;
      e=d=b-a;
    }
    if (fabs(fc) < fabs(fb)) {
      a=b;
      b=c;
      c=a;
      fa=fb;
      fb=fc;
      fc=fa;
    }
    tol1=2.0*EPS*fabs(b)+0.5*brent_solve_tolerance;
    xm=0.5*(c-b);
    if (fabs(xm) <= tol1 || fb == 0.0) return b;
    if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
      s=fb/fa;
      if (a == c) {
         p=2.0*xm*s;
         q=1.0-s;
      } else {
	q=fa/fc;
	r=fb/fc;
	p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
	q=(q-1.0)*(r-1.0)*(s-1.0);
      }
      if (p > 0.0) q = -q;
      p=fabs(p);
      min1=3.0*xm*q-fabs(tol1*q);
      min2=fabs(e*q);
      if (2.0*p < (min1 < min2 ? min1 : min2)) {
	e=d;
	d=p/q;
      } else {
	d=xm;
	e=d;
      }
    } else {
      d=xm;
      e=d;
    }
    a=b;
    fa=fb;
    if (fabs(d) > tol1)
      b += d;
    else
      b += SIGN(tol1,xm);

    if (gas_state_ptr->use_nlte_ == 0)
      {
	fb=rad_eq_function_LTE(gas_state_ptr, cell,b,solve_flag,solve_error);
      }
    else
      {
	fb=rad_eq_function_NLTE(gas_state_ptr, cell,b,solve_flag,solve_error);
      }
  }
  printf("Maximum number of iterations exceeded in zbrent");
  return 0.0;
}
