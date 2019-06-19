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

int transport::solve_state_and_temperature(int i)
{
  
  int solve_error = 0;
  if (set_gas_temp_to_rad_temp_ == 1)
    {
      grid->z[i].T_gas = pow(grid->z[i].e_rad/pc::a,0.25);
      return 0;
    }
  else
    {
      grid->z[i].T_gas = temp_brent_method(i,1,solve_error); // gas_state solve will happen in here. The solve error is for the gas_state solve

        if (gas_state_.use_nlte_)
	  {
	    bf_heating[i] = gas_state_.bound_free_heating_rate(grid->z[i].T_gas,J_nu_[i]);
	    ff_heating[i] = gas_state_.free_free_heating_rate(grid->z[i].T_gas,J_nu_[i]);
	    bf_cooling[i] = gas_state_.bound_free_cooling_rate(grid->z[i].T_gas);
	    ff_cooling[i] = gas_state_.free_free_cooling_rate(grid->z[i].T_gas);
	    coll_cooling[i] = gas_state_.collisional_net_cooling_rate(grid->z[i].T_gas);
	  }
      return solve_error;
    }


}

void transport::solve_eq_temperature()
{
  int solve_error = 0;
  for (int i=my_zone_start_;i<my_zone_stop_;i++)
    {
      if (set_gas_temp_to_rad_temp_ == 1)
	grid->z[i].T_gas = pow(grid->z[i].e_rad/pc::a,0.25);
      else
	{
	  grid->z[i].T_gas = temp_brent_method(i,0,solve_error); // solve_error won't be updated here because that's for the gas_state solve which isn't happening here

	  if (gas_state_.use_nlte_)
	    {
	      bf_heating[i] = gas_state_.bound_free_heating_rate(grid->z[i].T_gas,J_nu_[i]);
	      ff_heating[i] = gas_state_.free_free_heating_rate(grid->z[i].T_gas,J_nu_[i]);
	      bf_cooling[i] = gas_state_.bound_free_cooling_rate(grid->z[i].T_gas);
	      ff_cooling[i] = gas_state_.free_free_cooling_rate(grid->z[i].T_gas);
	      coll_cooling[i] = gas_state_.collisional_net_cooling_rate(grid->z[i].T_gas);
	  }
	
	}

      
    }
}


//***************************************************************/
// This is the function that expresses radiative equillibrium
// in a cell (i.e. E_absorbed = E_emitted).  It is used in
// The Brent solver below to determine the temperature such
// that RadEq holds
//**************************************************************/
double transport::rad_eq_function_LTE(int c,double T, int solve_flag, int & solve_error)
{

  zone* z = &(grid->z[c]);
  gas_state_.dens_ = z->rho;
  gas_state_.temp_ = T;
  
  if (solve_flag)
    {
	solve_error = gas_state_.solve_state(J_nu_[c]);
    }
  // total energy absorbed in zone
  double E_absorbed = grid->z[c].e_abs; // debug + grid->z[c].L_radio_dep;
  // total energy emitted (to be calculated)
  double E_emitted = 0.;

  // Calculate total emission assuming no frequency (grey) opacity
  if (nu_grid.size() == 1)
  {
    E_emitted = 4.0*pc::pi*abs_opacity_[c][0]*pc::sb/pc::pi*pow(T,4);
  }
  // integrate emisison over frequency (angle
  // integration gives the 4*PI) to get total
  // ergs/sec/cm^3 radition emitted. Opacities are
  // held constant for this (assumed not to change
  // much from the last time step).
  else for (int i=0;i<nu_grid.size();i++)
  {
    double dnu  = nu_grid.delta(i);
    double nu   = nu_grid.center(i);
    double B_nu = blackbody_nu(T,nu);
    double kappa_abs = abs_opacity_[c][i];
    E_emitted += 4.0*pc::pi*kappa_abs*B_nu*dnu;
    //Eab += 4.0*pc::pi*kappa_abs*J_nu_[c][i]*dnu;
  }
  //if (verbose) std::cout << c << " " << E_emitted << " " << Eab << " " << grid->z[c].e_abs << " " << grid->z[c].L_radio_dep << "\n";

  //std::cout << E_emitted << " " << E_absorbed << "\n";
  // radiative equillibrium condition: "emission equals absorbtion"
  // return to Brent function to iterate this to zero
  return (E_emitted - E_absorbed);
  

}


double transport::rad_eq_function_NLTE(int c,double T, int solve_flag, int &solve_error)
{

  zone* z = &(grid->z[c]);
  gas_state_.dens_ = z->rho;
  gas_state_.temp_ = T;
  
  if ( (gas_state_.smooth_grey_opacity_ == 1) || (gas_state_.use_zone_dependent_grey_opacity_ == 1) )
      {
	std::cerr << "# ERROR: NLTE solve should not be used with grey opacity\n";
	exit(1);
      }
      
  if (solve_flag)
    {
	solve_error = gas_state_.solve_state(J_nu_[c]);
    }
  
  double E_absorbed = gas_state_.free_free_heating_rate(T,J_nu_[c]) + gas_state_.bound_free_heating_rate(T,J_nu_[c]) ;

      // total energy emitted (to be calculated)
  double E_emitted =  E_emitted= gas_state_.free_free_cooling_rate(T) + gas_state_.bound_free_cooling_rate(T);

  if (gas_state_.use_collisions_nlte_)
    {
      printf("collisions\n");
      E_emitted += gas_state_.collisional_net_cooling_rate(T);
    }
  
    //std::cout << E_emitted << " " << E_absorbed << "\n";
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
double transport::temp_brent_method(int cell, int solve_flag, int &solve_error)
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
  if (gas_state_.use_nlte_ == 0)
    {
      fa=rad_eq_function_LTE(cell,a,solve_flag,solve_error);
      fb=rad_eq_function_LTE(cell,b,solve_flag,solve_error);
    }
  else 
    {
      fa=rad_eq_function_NLTE(cell,a,solve_flag,solve_error);
      fb=rad_eq_function_NLTE(cell,b,solve_flag,solve_error);
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

    if (gas_state_.use_nlte_ == 0)
      {
	fb=rad_eq_function_LTE(cell,b,solve_flag,solve_error);
      }
    else
      {
	fb=rad_eq_function_NLTE(cell,b,solve_flag,solve_error);
      }
  }
  printf("Maximum number of iterations exceeded in zbrent");
  return 0.0;
}
