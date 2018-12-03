#include <limits>
#include "AtomicSpecies.h"
#include "physical_constants.h"
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_linalg.h>
#include <iostream>

using namespace std;

// ---------------------------------------------------
// For the NLTE problem, we are solving a matrix equation
// M x = b
// where
//   x is the vector of the level population fractions
//   M is the rate matrix
//   and b is the zero vector assuming statistical equilibrium.
//
// the number density in each level is n_i = x_i*n_tot
// where n_tot is the total number density of the species
//
// Note: one of the rate equations is not independent,
// so in order for the matrix to be non-singular, we need to
// make the last equation express number conservation
//   sum_i  x_i = 1
// ---------------------------------------------------

namespace pc = physical_constants;

AtomicSpecies::AtomicSpecies()
{
  e_gamma_            = 0;
  no_ground_recomb_   = 0;
  use_betas_          = 0;
  minimum_extinction_ = 0;
  use_nlte_           = 0;

  // debug -- I hard coded this for now...
  min_level_pop_ = 1e-30;
}


//-------------------------------------------------------
// solve for the state of the atom
// determining ionization and level populations
// will use non-lte if use_nlte_ = True
// otherwise will use lte
// passed the electron density ne
//-------------------------------------------------------
int AtomicSpecies::solve_state(double ne)
{
  double status = 0;
  if (use_nlte_)
    status = solve_nlte(ne);
  else
    status = solve_lte(ne);
  return status;
}

//-------------------------------------------------------
// solve for LTE level populations and
// ionization state
// passed the electron density ne
//-------------------------------------------------------
int AtomicSpecies::solve_lte(double ne)
{

  // calculate partition functions
  for (int i=0;i<n_ions_;++i) ions_[i].part = 0;
  for (int i=0;i<n_levels_;++i)
  {
    levels_[i].n = levels_[i].g*exp(-levels_[i].E/pc::k_ev/gas_temp_);
    ions_[levels_[i].ion].part += levels_[i].n;
  }

  // thermal debroglie wavelength, lam_t**3
  double lt = pc::h*pc::h/(2.0*pc::pi*pc::m_e*pc::k*gas_temp_);
  double fac = 2/ne/pow(lt,1.5);

  // calculate saha ratios
  ions_[0].frac = 1.0;
  double norm  = 1.0;
  for (int i=1;i<n_ions_;++i)
  {
    // calculate the ratio of i to i-1
    double saha = exp(-1.0*ions_[i-1].chi/pc::k_ev/gas_temp_);
    saha = saha*(ions_[i].part/ions_[i-1].part)*fac;

    // set relative ionization fraction
    ions_[i].frac = saha*ions_[i-1].frac;

    // check for ridiculously small numbers
    if (ne < 1e-50) ions_[i].frac = 0;
    norm += ions_[i].frac;
  }
  // renormalize ionization fractions_
  for (int i=0;i<n_ions_;++i) ions_[i].frac = ions_[i].frac/norm;

  // calculate level densities (bolztmann factors)
  for (int i=0;i<n_levels_;++i)
  {
    double E = levels_[i].E;
    int    g = levels_[i].g;
    double Z = ions_[levels_[i].ion].part;
    double f = ions_[levels_[i].ion].frac;
    levels_[i].n = f*g*exp(-E/pc::k_ev/gas_temp_)/Z;
    levels_[i].n_lte = levels_[i].n;
    levels_[i].b = 1;
  }

  // this should return != 0 if failure really
  return 0;

}

//-------------------------------------------------------
// solve for non-LTE level populations and
// ionization state
// passed the electron density ne
//-------------------------------------------------------
int AtomicSpecies::solve_nlte(double ne)
{
  // initialize with LTE populations
  solve_lte(ne);

  // make sure our lte levels aren't too small
  for (int i=0;i<n_levels_;++i)
    if (levels_[i].n_lte < min_level_pop_) levels_[i].n_lte = min_level_pop_;

  // Calculate all of the transition rates
  set_rates(ne);

  // zero out matrix and vectors
  gsl_matrix_set_zero(M_nlte_);
  gsl_vector_set_zero(b_nlte_);
  gsl_vector_set_zero(x_nlte_);
  gsl_permutation_init(p_nlte_);

  // set up diagonal elements of rate matrix
  for (int i=0;i<n_levels_;++i)
  {
    double Rout = 0.0;
    // don't worry i = j rate should be zero
    for (int j=0;j<n_levels_;++j)
      Rout += rates_[i][j];
    Rout = -1*Rout;
    gsl_matrix_set(M_nlte_,i,i,Rout);
  }

  // set off diagonal elements of rate matrix
  for (int i=0;i<n_levels_;++i)
    for (int j=0;j<n_levels_;++j)
      if (i != j) gsl_matrix_set(M_nlte_,i,j,rates_[j][i]);

  // last row expresses number conservation
  for (int i=0;i<n_levels_;++i)
    gsl_matrix_set(M_nlte_,n_levels_-1,i,levels_[i].n_lte);
  gsl_vector_set(b_nlte_,n_levels_-1,1.0);

    //printf("----\n");
    //for (int i=0;i<n_levels_;++i)
    //  for (int j=0;j<n_levels_;++j)
    //   printf("%5d %5d %14.3e\n",i,j,gsl_matrix_get(M_nlte_,i,j));
    // printf("----\n");

  // solve rate matrix
  int status;
  gsl_linalg_LU_decomp(M_nlte_, p_nlte_, &status);
  gsl_linalg_LU_solve(M_nlte_, p_nlte_, b_nlte_, x_nlte_);

  // the x vector should now have the solved level
  // depature coefficients
  for (int i=0;i<n_levels_;++i)
  {
    double b = gsl_vector_get(x_nlte_,i);
    double n_nlte = b*levels_[i].n_lte;

    // make sure our solved for nlte levels_ aren't too small
    if (n_nlte < min_level_pop_) n_nlte = min_level_pop_;

    levels_[i].n = n_nlte;
    levels_[i].b = b;
  }

  // set the ionization fraction
  for (int i=0;i<n_ions_;++i)
    ions_[i].frac = 0;
  for (int i=0;i<n_levels_;++i)
    ions_[levels_[i].ion].frac += levels_[i].n;

  // this should return != 0 if failure really
  return 0;
}


//-------------------------------------------------------
// integrate up the radiation field over  lines
// to get the line J and over bound-free to get the
// photoionization rates
//-------------------------------------------------------
void AtomicSpecies::calculate_radiative_rates(std::vector<real> J_nu)
{

  // zero out recombination/photoionization rates
  for (int j=0;j<n_levels_;++j)
  {
    levels_[j].P_ic = 0;
    levels_[j].R_ci = 0;
  }

  // calculate photoionization/recombination rates
  // photoionization is Eq. 3.7 in Rutten, Stellar Atmopsheres
  // recombination is from EQ 3.16 of Rutten, stellar atmospheres
  // recombination rate includes stimulated recombination
  double fac1 = 2/pc::c/pc::c;
  int ng = nu_grid_.size();
  for (int i=1;i<ng;++i)
  {
    double nu     = nu_grid_.center(i);
    double E_ergs = pc::h*nu;
    double E_ev   = E_ergs*pc::ergs_to_ev;
    double J      = J_nu[i];
    double dnu    = nu_grid_.delta(i);
    for (int j=0;j<n_levels_;++j)
    {
      int ic = levels_[j].ic;
      if (ic == -1) continue;
      double chi = levels_[j].E_ion;
      if (E_ev < chi) continue;

      // photoionization term
      double sigma = levels_[j].s_photo.value_at_with_zero_edges(E_ev);
      double Jterm = sigma*J/E_ergs;
      levels_[j].P_ic += Jterm*dnu;

      // recombination term
      levels_[j].R_ci += (sigma*fac1*nu*nu + Jterm)*exp(-1.0*(E_ev - chi)/pc::k_ev/gas_temp_)*dnu;

    }
  }

  // multiply by overall factors
  double lam_t = sqrt(pc::h*pc::h/(2*pc::pi*pc::m_e*pc::k*gas_temp_));
  double saha_fac = lam_t*lam_t*lam_t/2.0;
  for (int j=0;j<n_levels_;++j)
  {
    int ic = levels_[j].ic;
    if (ic == -1) continue;
    double gl_o_gc = (1.0*levels_[j].g)/(1.0*levels_[ic].g);
    levels_[j].P_ic *= 4*pc::pi;
    levels_[j].R_ci *= 4*pc::pi*gl_o_gc*saha_fac;

//    std::cout << "Pic = " << levels_[j].P_ic << " ; R_ci " << levels_[j].R_ci << "\n";

  }

  //for (int j=0;j<n_levels_;++j)
   // if (levels_[j].ic != -1)
    //  levels_[j].R_ci = 2.58e-13;
    //std::cout << j << " " << levels_[j].P_ic << " " << levels_[j].R_ci << "\n";

  // calculate line J's
  double x_max = 100;
  double dx    = 0.05;

  for (int i=0;i<n_lines_;++i)
  {
    double nu0 = lines_[i].nu;
      double sum  = 0;
    double J0   = 0;

    double nu_d    = nu0*line_beta_dop_;
    double gamma   = lines_[i].A_ul;
    double a_voigt = gamma/4/pc::pi/nu_d;

    for (double x=-1*x_max;x<=x_max;x+=dx)
    {
      double phi = voigt_profile_.getProfile(x,a_voigt);
      double n = nu0 + x*nu_d;
      double J1 = nu_grid_.value_at(n,J_nu)*phi;
      sum += 0.5*(J1 + J0)*dx;
      J0 = J1;
    }
    lines_[i].J = sum;
  }
}



//-------------------------------------------------------
// Set the rates for all possible transitions
//------------------------------------------------------
void AtomicSpecies::set_rates(double ne)
{
  // zero out rate matrix
  for (int i=0;i<n_levels_;++i)
    for (int j=0;j<n_levels_;++j)
      rates_[i][j] = 0;

  // ------------------------------------------------
  // radiative bound-bound transitions
  // ------------------------------------------------

  for (int l=0;l<n_lines_;l++)
  {
    int lu     = lines_[l].lu;
    int ll     = lines_[l].ll;

    // spontaneous dexcitation + stimulated emission
    double R_ul = lines_[l].B_ul*lines_[l].J + lines_[l].A_ul;
    double R_lu = lines_[l].B_lu*lines_[l].J;

    // check for transition between degenerate levels_
    if (lines_[l].nu == 0)
      { R_ul = 0; R_lu = 0;}

    // add into rates
    rates_[ll][lu] += R_lu;
    rates_[lu][ll] += R_ul;

   // printf("RR %d %d %e\n",ll,lu,R_lu);
   // printf("RR %d %d %e\n",lu,ll,R_ul);
  }

  double norm = 0;
  for (int l=0;l<n_lines_;l++) norm   += lines_[l].f_lu;

  for (int l=0;l<n_lines_;l++)
  {
    int lu  = lines_[l].lu;
    int ll  = lines_[l].ll;

    // ------------------------------------------------
    // non-thermal (radioactive) bound-bound transitions
    // ------------------------------------------------

    double dE = (levels_[lu].E - levels_[ll].E)*pc::ev_to_ergs;
    double R_lu = e_gamma_/n_dens_/dE; //*(lines_[l].f_lu/norm);
    if (dE == 0) R_lu = 0;
    if (ll != 0) R_lu = 0;

    // add into rates
    // debug -- turning off non-thermal rates for now
    //rates_[ll][lu] += R_lu;

    // printf("GR %d %d %e %e %e\n",ll,lu,R_lu,e_gamma,dE);

    // ------------------------------------------------
    // collisional bound-bound transitions
    // ------------------------------------------------

    double zeta = dE/pc::k/gas_temp_; // note dE is in ergs
    double ezeta = exp(zeta);
    double C_up = 3.9*pow(zeta,-1.)*pow(gas_temp_,-1.5) / ezeta * ne * lines_[l].f_lu;
    if (zeta > 700) C_up = 0.; // be careful about overflow

    double C_down = 3.9*pow(zeta,-1.)*pow(gas_temp_,-1.5) * ne * lines_[l].f_lu * levels_[ll].g/levels_[lu].g;

    rates_[ll][lu] += C_up;
    rates_[lu][ll] += C_down;

  }


  // ------------------------------------------------
  // bound-free transitions
  // ------------------------------------------------
  for (int i=0;i<n_levels_;++i)
  {
    int ic = levels_[i].ic;
    if (ic == -1) continue;

    // ionization potential
    int istage  = levels_[i].ion;
    double chi  = ions_[istage].chi - levels_[i].E;
    double zeta = chi/pc::k_ev/gas_temp_;

    // collisional ionization rate
    // needs to be multiplied by number of electrons in outer shell
    double C_ion = 2.7/zeta/zeta*pow(gas_temp_,-1.5)*exp(-zeta)*ne;
    rates_[i][ic] += C_ion;

    // collisional recombination rate
    int gi = levels_[i].g;
    int gc = levels_[ic].g;
    double C_rec = 5.59080e-16/zeta/zeta*pow(gas_temp_,-3)*gi/gc*ne*ne;
    rates_[ic][i] += C_rec;

    // radiative recombination rate (debug)
    //double R_rec = ne*levels_[i].a_rec.value_at(T);
    //// suppress recombinations to ground
    //if (no_ground_recomb) if (levels_[i].E == 0) R_rec = 0;

    // photoionization and radiative recombination
    rates_[ic][i] += levels_[i].R_ci*ne;
    rates_[i][ic] += levels_[i].P_ic;

    //printf("pc::pi: %d %d %e\n",i,ic,levels_[i].P_ic);
    //printf("CI: %d %d %e %e\n",i,ic,C_rec,C_ion);

  }

   // multiply by rates by lte pop in level coming from
  // (becuase we will solve for depature coeffs)
  for (int i=0;i<n_levels_;++i)
      for (int j=0;j<n_levels_;++j)
        rates_[i][j] *= levels_[i].n_lte;

  // print out rates if you so like
  //printf("------- rates ----------\n");
  for (int i=0;i<n_levels_;++i)
    for (int j=0;j<n_levels_;++j)
      if (isnan(rates_[i][j]))
       printf("%5d %5d %14.5e\n",i,j,rates_[i][j]);
   //printf("\n");
   //for (int i=0;i<n_levels_;++i)
   // for (int j=0;j<n_levels_;++j)   {
   //     if (isnan(rates_[i][j])) std::cout << "NAN RATE\n";
   //     if (isinf(rates_[i][j])) std::cout << "INF RATE: " << i << " - " << j << "\n"; }
}


double AtomicSpecies::get_ion_frac()
{
  double x = 0;
  for (int i=0;i<n_levels_;++i)
    x += levels_[i].n*levels_[i].ion;
  return x;
}



double AtomicSpecies::Calculate_Milne(int lev, double temp)
{
  // Maxwell-Bolztmann constants
  double v_MB = sqrt(2*pc::k*temp/pc::m_e);
  double MB_A = 4/sqrt(pc::pi)*pow(v_MB,-3);
  double MB_B = pc::m_e/pc::k/2.0/temp;
  double milne_fac = pow(pc::h/pc::c/pc::m_e,2);

  // starting values
  double sum   = 0;
  double nu_t  = levels_[lev].E_ion*pc::ev_to_ergs/pc::h;
  double nu    = nu_t;
  double vel   = 0;
  double fMB   = 0;
  double sigma = 0;
  double coef  = 0;
  double old_vel  = vel;
  double old_coef = coef;

  // integrate over velocity/frequency
  for (int i=1;i<levels_[lev].s_photo.size();++i)
  {
    // recombination cross-section
    double E = levels_[lev].s_photo.x[i];
    double S = levels_[lev].s_photo.y[i];
    nu       = E*pc::ev_to_ergs/pc::h;
    vel      = sqrt(2*pc::h*(nu - nu_t)/pc::m_e);
    if (nu < nu_t) vel = 0;
    fMB   = MB_A*vel*vel*exp(-MB_B*vel*vel);
    sigma = milne_fac*S*nu*nu/vel/vel;
    coef  = vel*sigma*fMB;

    // integrate
    sum += 0.5*(coef + old_coef)*(vel - old_vel);
    // store old values
    old_vel  = vel;
    old_coef = coef;
  }

  // ionize to state
  int ic = levels_[lev].ic;
  if (ic == -1) return 0;

  // return value
  //return levels_[lev].g/levels_[ic].g*sum;
  return (1.0*levels_[lev].g)/(1.0*levels_[ic].g)*sum;
}


void AtomicSpecies::print()
{

  cout << "--------------------- ions; n = " << n_ions_ << " ---------------------\n";
  cout << "# ion \t part \t frac \t chi (eV)\n";
  cout << "#---------------------------------------------------------------\n";


  for (int i=0;i<n_ions_;++i)
    cout << "#   "<<  ions_[i].stage << "\t" << ions_[i].part << "\t" << ions_[i].frac
	 << "\t" << ions_[i].chi << endl;


  cout << "\n";
  cout << "--------------------------------------------------------------------\n";
  cout << "--------------------levels; n = " << n_levels_ << "------------------------\n";
  cout << "# lev   ion     E_ex        g      pop          b_i       ion_to\n";
  cout << "#---------------------------------------------------------------\n";

  for (int i=0;i<n_levels_;++i)
  {
    printf("%5d %4d %12.3e %5d %12.3e %12.3e %5d\n",
	   levels_[i].globalID,levels_[i].ion,
	   levels_[i].E,levels_[i].g,levels_[i].n,
	   levels_[i].b,levels_[i].ic);
  }

  printf("\n--- line data\n");

  for (int i=0;i<n_lines_;++i)
  {
    printf("%8d %4d %4d %12.3e %12.3e %12.3e %12.3e %12.3e\n",
    	   i,lines_[i].ll,lines_[i].lu,lines_[i].nu,lines_[i].f_lu,
    	   lines_[i].A_ul,lines_[i].B_ul,lines_[i].B_lu);
  }

  printf("\n--- line optical depths\n");

  for (int i=0;i<n_lines_;++i)
  {
    int    ll = lines_[i].ll;
    double nl = levels_[ll].n;

    printf("%8d %4d %4d %12.3e %12.3e %12.3e\n",
    	   i,lines_[i].ll,lines_[i].lu,lines_[i].nu,
	   lines_[i].tau,nl);
  }


}

//-----------------------------------------------------------------
// calculate planck function in frequency units
//-----------------------------------------------------------------
double AtomicSpecies::blackbody_nu(double T, double nu)
{
  double zeta = pc::h*nu/pc::k/T;
  return 2.0*nu*nu*nu*pc::h/pc::c/pc::c/(exp(zeta)-1);
}
