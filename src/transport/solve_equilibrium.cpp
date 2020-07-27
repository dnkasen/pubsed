#include <stdio.h>
#include <math.h>
#include <iostream>
#include "transport.h"
#include "physical_constants.h"
#include "grid_general.h"
#include "brent.h"
#include "brent.cpp"

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

    transportMemFn f;

    // For LTE, do an initial solve of the gas state
    if (gas_state_ptr->use_nlte_ == 0)
    {
      solve_error = gas_state_ptr->solve_state();
      gas_state_ptr->computeOpacity(abs_opacity_[i],scat,emis);
      f = &transport::rad_eq_wrapper_LTE;
    }
    else f = &transport::rad_eq_wrapper_NLTE;

    // Calculate equilibrium temperature.
    // Additional gas_state solve may also happen here
    //    grid->z[i].T_gas = temp_brent_method(gas_state_ptr, i,1,solve_error);

    brent_args.gas_state_ptr = gas_state_ptr;
    brent_args.c = i;
    brent_args.solve_flag = 1;
    brent_args.solve_error = &solve_error;

    brent_solver<transport> solver;

    int n; // will store number of brent solver iterations
    // still using hard-coded eps; that could be set here
    // lower bracket and uppr bracket have been set in .lua files
    // Calculate equilibrium temperature.
    // Additional gas_state solve may also happen here
    grid->z[i].T_gas = solver.solve(*this, f, temp_min_value_,temp_max_value_,0.001, &n);

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
#pragma omp parallel for schedule(dynamic, 16) default(none) firstprivate(gas_state_ptr, solve_error)
  for (int i=my_zone_start_;i<my_zone_stop_;i++)
  {
    if (set_Tgas_to_Trad_ == 1)
    grid->z[i].T_gas = pow(grid->z[i].e_rad/pc::a,0.25);
    else
  {

    transportMemFn f;
    if (gas_state_ptr->use_nlte_ == 0)
    {
      f = &transport::rad_eq_wrapper_LTE;
    }
    else f = &transport::rad_eq_wrapper_NLTE;

    
    brent_args.gas_state_ptr = gas_state_ptr;
    brent_args.c = i;
    brent_args.solve_flag = 0;
    brent_args.solve_error = &solve_error;  // solve_error won't actually be updated here because that's for the gas_state solve which isn't happening here

    brent_solver<transport> solver;
    int n; // will store number of brent solver iterations
    // still using hard-coded eps; that could be set here
    // lower bracket and uppr bracket have been set in .lua files
    double T_solution = solver.solve(*this, f, temp_min_value_,temp_max_value_,0.001, &n);
    //    solve_error = *(brent_args.solve_error); // see abovee

    grid->z[i].T_gas = T_solution;

    // solve_error won't be updated here because that's for the gas_state solve which isn't happening here
    //     grid->z[i].T_gas = temp_brent_method(gas_state_ptr, i,0,solve_error);

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

double transport::rad_eq_wrapper_LTE(double T)
{
  return rad_eq_function_LTE(brent_args.gas_state_ptr, brent_args.c, T, brent_args.solve_flag, brent_args.solve_error);
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
double transport::rad_eq_function_LTE(GasState* gas_state_ptr, int c,double T, int solve_flag, int *solve_error)
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
    // if you want to resolve the LTE state on each T iteration:
    // *solve_error = gas_state_ptr->solve_state();
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


double transport::rad_eq_wrapper_NLTE(double T)
{
  return rad_eq_function_NLTE(brent_args.gas_state_ptr, brent_args.c, T, brent_args.solve_flag, brent_args.solve_error);
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
double transport::rad_eq_function_NLTE(GasState* gas_state_ptr, int c,double T, int solve_flag, int *solve_error)
{

  zone* z = &(grid->z[c]);
  gas_state_ptr->dens_ = z->rho;
  gas_state_ptr->temp_ = T;

  // make sure grey_opacity is not being used
  if (gas_state_ptr->total_grey_opacity_ != 0)
  {
    std::cerr << "# ERROR: NLTE solve should not be used with grey opacity\n";
    exit(1);
  }

  // if flag set, recompute the entire NLTE problem for this iteration
  if (solve_flag)
	*solve_error = gas_state_ptr->solve_state(J_nu_[c]);

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

