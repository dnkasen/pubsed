#include <stdio.h>
#include <math.h>
#include <iostream>
#include "transport.h"
#include "physical_constants.h"
#include "grid_general.h"
#include "brent.h"

namespace pc = physical_constants;

//-------------------------------------------------------------
//  Solve for the temperature assuming radiative equilibrium
//  This function may not work
//-------------------------------------------------------------
int transport::solve_state_and_temperature(GasState* gas_state_ptr, int i)
{

  vector<OpacityType> emis(nu_grid_.size());
  vector<OpacityType> scat(nu_grid_.size());
  emis.assign(emis.size(),0.0);

  int solve_root_errors_temp = 0; // store total such errors
  int solve_iter_errors_temp = 0; // store total such errors
  int solve_root_errors_ne = 0; // store total such errors
  int solve_iter_errors_ne = 0; // store total such errors

  int gas_solve_error; // stores type of error for single gas solve

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
    if (gas_state_ptr->is_nlte_turned_on() == 0)
    {
      gas_solve_error = gas_state_ptr->solve_state();

      if (gas_solve_error == 1)
	{
	  //	  printf("root not bracketed in n_e solve\n");
	  solve_root_errors_ne++;
	}
      if (gas_solve_error == 2)
	{
	  //	  printf("max number of iterations reached in n_e solve\n");
	  solve_iter_errors_ne++;
	}
      gas_state_ptr->computeOpacity(abs_opacity_[i],scat,emis);
      f = &transport::rad_eq_function_LTE;
    }
    else
      {
	f = &transport::rad_eq_function_NLTE;
      }


    radeq_brent_args brent_args;

    brent_args.gas_state_ptr = gas_state_ptr;
    brent_args.c = i;
    brent_args.solve_flag = 1;
    brent_args.solve_error = &gas_solve_error;

    brent_solver<transport,radeq_brent_args> solver;
    double tol    = 1.e-3;
    int max_iters = 100;

    int n; // will store number of brent solver iterations

    // Calculate equilibrium temperature.
    // lower bracket and upper bracket have been set in .lua files (temperature min and max)
    grid->z[i].T_gas = solver.solve(*this, f, &brent_args,temp_min_value_,temp_max_value_,tol, max_iters, &n);

    if (n == -1)
      {
	//	printf("root not bracketed in temperature solve\n");
	solve_root_errors_temp++;
      }
    if (n == -2)
      {
	//	printf("max number of iterations reached in temperature solve\n");
	solve_iter_errors_temp++;
      }
    if (gas_solve_error == 1)
      {
	solve_root_errors_ne++;
      }
    if (gas_solve_error == 2)
      {
	//	printf("max number of iterations reached in temperature solve\n");
	solve_iter_errors_temp++;
      }

    if (gas_state_ptr->is_nlte_turned_on() == 0)
    {
      // For LTE, do a final solve
      gas_solve_error = gas_state_ptr->solve_state();
      if (gas_solve_error == -1)
	{
	  //printf("root not bracketed in n_e solve\n");
	  solve_root_errors_ne++;
	}
      if (gas_solve_error == -2)
	{
	  //printf("max number of iterations reached in n_e solve\n");
	  solve_iter_errors_ne++;
	}
    }

    if (gas_state_ptr->is_nlte_turned_on())
    {
      bf_heating[i] = gas_state_ptr->bound_free_heating_rate(grid->z[i].T_gas,J_nu_cmf[i]);
      ff_heating[i] = gas_state_ptr->free_free_heating_rate(grid->z[i].T_gas,J_nu_cmf[i]);
      bf_cooling[i] = gas_state_ptr->bound_free_cooling_rate(grid->z[i].T_gas);
      ff_cooling[i] = gas_state_ptr->free_free_cooling_rate(grid->z[i].T_gas);
      coll_cooling[i] = gas_state_ptr->collisional_net_cooling_rate(grid->z[i].T_gas);
    }

    if (solve_root_errors_temp > 0)
      std::cerr << "# WARNING: temperature not bracketed in " << solve_root_errors_temp << " zones" << std::endl;
    if (solve_iter_errors_temp > 0)
      std::cerr << "# WARNING: max iterations hit in temperature solve in " << solve_iter_errors_temp << " zones" << std::endl;
    if (solve_root_errors_ne > 0)
      std::cerr << "# WARNING: n_e not bracketed in " << solve_root_errors_ne << " zones" << std::endl;
    if (solve_iter_errors_ne > 0)
      std::cerr << "# WARNING: max iterations hit in n_e solve in " << solve_iter_errors_ne << " zones" << std::endl;

    gas_solve_error = solve_root_errors_ne +  solve_iter_errors_ne + solve_root_errors_temp + solve_iter_errors_temp ;
    return gas_solve_error;
  }

}

//-------------------------------------------------------------
//  Solve for the temperature assuming radiative equilibrium
//-------------------------------------------------------------
void transport::solve_eq_temperature()
{
  int solve_error = 0;
  int solve_root_errors_temp = 0;
  int solve_iter_errors_temp = 0;
  int solve_root_errors_ne = 0;
  int solve_iter_errors_ne = 0;
#pragma omp parallel default(none) firstprivate(solve_error) shared(solve_root_errors_temp,solve_iter_errors_temp,solve_root_errors_ne,solve_iter_errors_ne)
  {
#ifdef _OPENMP
  int tid = omp_get_thread_num();
#else
  int tid = 0;
#endif
  GasState* gas_state_ptr = &(gas_state_vec_[tid]);
#pragma omp for schedule(dynamic,16)
  for (int i=my_zone_start_;i<my_zone_stop_;i++)
  {
    if (set_Tgas_to_Trad_ == 1)
      grid->z[i].T_gas = pow(grid->z[i].e_rad/pc::a,0.25);
    else
    {
      // fill up GasState ptr and solve it's state
      solve_error = fill_and_solve_gasstate(gas_state_ptr, i);

      if (solve_error == -1)
	{
	  //	  printf("root not bracketed in n_e solve\n");
	  solve_root_errors_ne++;
	}
      if (solve_error == -2)
	{
	  //	  printf("max number of iterations reached in n_e solve\n");
	  solve_iter_errors_ne++;
	}

      transportMemFn f;
      if (gas_state_ptr->is_nlte_turned_on() == 0)
	{
	  f = &transport::rad_eq_function_LTE;
	}
      else f = &transport::rad_eq_function_NLTE;

      radeq_brent_args brent_args;

      brent_args.gas_state_ptr = gas_state_ptr;
      brent_args.c = i;
      brent_args.solve_flag = 0;
      brent_args.solve_error = &solve_error;  // solve_error won't actually be updated here because that's for the gas_state solve which isn't happening here

      brent_solver<transport,radeq_brent_args> solver;
      double tol = 1.e-3;
      int max_iters = 100;
      int n; // will store number of brent solver iterations
      // lower bracket and upper bracket have been set in .lua files (min and max temp)
      grid->z[i].T_gas = solver.solve(*this, f, &brent_args, temp_min_value_,temp_max_value_,tol, max_iters, &n);

      if (n == -1)
	{
	  //	  printf("root not bracketed in n_e solve\n");
	  solve_root_errors_temp++;
	}
      if (n == -2)
	{
	  //	  printf("max number of iterations reached in n_e solve\n");
	  solve_iter_errors_temp++;
	}

      if (gas_state_ptr->is_nlte_turned_on())
      {
	bf_heating[i] = gas_state_ptr->bound_free_heating_rate(grid->z[i].T_gas,J_nu_cmf[i]);
	ff_heating[i] = gas_state_ptr->free_free_heating_rate(grid->z[i].T_gas,J_nu_cmf[i]);
	bf_cooling[i] = gas_state_ptr->bound_free_cooling_rate(grid->z[i].T_gas);
	ff_cooling[i] = gas_state_ptr->free_free_cooling_rate(grid->z[i].T_gas);
	coll_cooling[i] = gas_state_ptr->collisional_net_cooling_rate(grid->z[i].T_gas);
      }
    }
  }
  }
    #pragma omp single
    if (verbose)
      {
  if (solve_root_errors_temp > 0)
    std::cerr << "# WARNING: temperature not bracketed in " << solve_root_errors_temp << " zones" << std::endl;
  if (solve_iter_errors_temp > 0)
    std::cerr << "# WARNING: max iterations hit in temperature solve in " << solve_iter_errors_temp << " zones" << std::endl;
  if (solve_root_errors_ne > 0)
    std::cerr << "# WARNING: n_e not bracketed in " << solve_root_errors_ne << " zones" << std::endl;
  if (solve_iter_errors_ne > 0)
    std::cerr << "# WARNING: max iterations hit in n_e solve in " << solve_iter_errors_ne << " zones" << std::endl;
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
double transport::rad_eq_function_LTE(double T, radeq_brent_args* args)
{

  zone* z = &(grid->z[args->c]);
  args->gas_state_ptr->dens_ = z->rho;
  args->gas_state_ptr->temp_ = T;

  // helper variables need for call (will not be used)
  vector<OpacityType> emis(nu_grid_.size());
  vector<OpacityType> scat(nu_grid_.size());
  emis.assign(emis.size(),0.0);

  // recalculate opacities based on current T if desired
  if (args->solve_flag)
  {
    // if you want to resolve the LTE state on each T iteration:
    // args->*solve_error = args->gas_state_ptr->solve_state();
    args->gas_state_ptr->computeOpacity(abs_opacity_[args->c],scat,emis);
  }

  // total energy emitted (to be calculated)
  double E_emitted = 0.;

  // total energy absorbed in zone
  double E_absorbed = 0.;
  // if not solve_flag, use the absorption rate estimated during transport
  if (args->solve_flag == 0)
    E_absorbed = grid->z[args->c].e_abs;

  // Calculate total emission assuming no frequency (grey) opacity
  if (nu_grid_.size() == 1)
  {
    E_emitted = 4.0*pc::pi*abs_opacity_[args->c][0]*pc::sb/pc::pi*pow(T,4);

    if (args->solve_flag == 0)
      E_absorbed = pc::c *abs_opacity_[args->c][0] * grid->z[args->c].e_rad;
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
    double kappa_abs = abs_opacity_[args->c][i];
    E_emitted += 4.0*pc::pi*kappa_abs*B_nu*dnu;
    if (args->solve_flag == 1)
      E_absorbed += 4.0*pc::pi*kappa_abs*J_nu_cmf[args->c][i]*dnu;
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
double transport::rad_eq_function_NLTE(double T, radeq_brent_args* args)
{

  zone* z = &(grid->z[args->c]);
  args->gas_state_ptr->dens_ = z->rho;
  args->gas_state_ptr->temp_ = T;

  // make sure grey_opacity is not being used
  if (z->total_grey_opacity != 0)
  {
    std::cerr << "# ERROR: NLTE solve should not be used with grey opacity\n";
    exit(1);
  }

  // if flag set, recompute the entire NLTE problem for this iteration
  if (args->solve_flag)
    {
     *(args->solve_error) = args->gas_state_ptr->solve_state(J_nu_cmf[args->c]);
    }

  // total energy absorbed
  double E_absorbed = args->gas_state_ptr->free_free_heating_rate(T,J_nu_cmf[args->c]) +
        args->gas_state_ptr->bound_free_heating_rate(T,J_nu_cmf[args->c]) ;
  // add in radioactive heating
  E_absorbed += args->gas_state_ptr->e_gamma_heat_;

  // total energy emitted
  double E_emitted= args->gas_state_ptr->free_free_cooling_rate(T) +
        args->gas_state_ptr->bound_free_cooling_rate(T);

  if (args->gas_state_ptr->use_collisions_nlte_)
      E_emitted += args->gas_state_ptr->collisional_net_cooling_rate(T);

  // radiative equillibrium condition: "emission equals absorbtion"
  // return to Brent function to iterate this to zero
  return (E_emitted - E_absorbed);
}
