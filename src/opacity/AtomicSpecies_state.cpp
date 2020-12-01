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
#include <cmath>

#if defined(USE_EIGEN) && USE_EIGEN == 1
#include <Eigen/Dense>
using namespace Eigen;
#endif

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
  // loop over level to calculate partition functions
  for (int i=0;i<n_ions_;++i) ion_part_[i] = 0;
  for (int i=0;i<n_levels_;++i)
  {
    double E = adata_->get_lev_E(i);
    int    g = adata_->get_lev_g(i);
    int  ion = adata_->get_lev_ion(i);
    lev_n_[i] = g*exp(-E/pc::k_ev/gas_temp_);
    ion_part_[ion] += lev_n_[i];
  }

  // thermal debroglie wavelength, lam_t**3
  double lt = pc::h*pc::h/(2.0*pc::pi*pc::m_e*pc::k*gas_temp_);
  double fac = 2/ne/pow(lt,1.5);

  // calculate saha ratios
  ion_frac_[0] = 1.;
  double norm  = 1.;
  for (int i=1;i<n_ions_;++i)
  {
    // calculate the ratio of i to i-1
    double chi  = adata_->get_ion_chi(i-1);
    double saha = exp(-1.0*chi/pc::k_ev/gas_temp_);
    saha = saha*(ion_part_[i]/ion_part_[i-1])*fac;

    // set relative ionization fraction
    ion_frac_[i] = saha*ion_frac_[i-1];

    // check for ridiculously small numbers
    if (ne < 1e-50) ion_frac_[i] = 0;

    // ReNormalize to keep numbers reasonable
    norm = 0;
    for (int j=0; j<=i; j++) {
      norm += ion_frac_[j];
    }
    for (int j=0; j<=i; j++) {
      ion_frac_[j] = ion_frac_[j] / norm;
    }
  }
  // renormalize ionization fractions
  norm = 0.;
  for (int i=0;i<n_ions_; i++)
      norm += ion_frac_[i];
  for (int i=0;i<n_ions_;++i) ion_frac_[i] = ion_frac_[i];

  // calculate level densities (bolztmann factors)
  for (int i=0;i<n_levels_;++i)
  {
    double E = adata_->get_lev_E(i);
    int    g = adata_->get_lev_g(i);
    int  ion = adata_->get_lev_ion(i);
    double Z = ion_part_[ion];
    double f = ion_frac_[ion];
    lev_n_[i]   = f*g*exp(-E/pc::k_ev/gas_temp_)/Z;

    // make sure our LTE level pops aren't too small
    if (lev_n_[i] < min_level_pop_)
      lev_n_[i] = min_level_pop_;
    // record the LTE population value
    lev_lte_[i] = lev_n_[i];
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
    gsl_matrix_set(M_nlte_,n_levels_-1,i,lev_lte_[i]);
  gsl_vector_set(b_nlte_,n_levels_-1,1.0);

    //printf("----\n");
    //for (int i=0;i<n_levels_;++i)
    //  for (int j=0;j<n_levels_;++j)
    //   printf("%5d %5d %14.3e\n",i,j,gsl_matrix_get(M_nlte_,i,j));
    // printf("----\n");


  #if defined(USE_EIGEN) && USE_EIGEN == 1
  // solve the rate matrix with eigen
  MatrixXd eigen_nlte(n_levels_,n_levels_);
  VectorXd eigen_b(n_levels_);
  VectorXd eigen_x(n_levels_);
  for (int i=0;i<n_levels_;++i) {
    for (int j=0;j<n_levels_;++j) {
      eigen_nlte(i,j) = gsl_matrix_get(M_nlte_,i,j);
    }
    eigen_b(i) = gsl_vector_get(b_nlte_,i);
  }

  eigen_x = eigen_nlte.fullPivLu().solve(eigen_b);

  // the x vector should now have the solved level
  // depature coefficients
  for (int i=0;i<n_levels_;++i)
  {
    lev_n_[i] = eigen_x[i]*lev_lte_[i];
    // check that level populations aren't too negative or too large
    if (lev_n_[i] < -1.0e-5 || lev_n_[i] > 1.00001) {
      printf("problem with NLTE level pops\n");
      for (int j=0;j<n_levels_;++j) {
        printf("lev: %5d pop: %14.3e\n", j, lev_n_[j]);
      }
      exit(1);
    }

    if (lev_n_[i] < min_level_pop_)
      lev_n_[i] = min_level_pop_;
  }

  #else

  int status;
  gsl_linalg_LU_decomp(M_nlte_, p_nlte_, &status);
  gsl_linalg_LU_solve(M_nlte_, p_nlte_, b_nlte_, x_nlte_);

    // the x vector should now have the solved level
  // depature coefficients
  for (int i=0;i<n_levels_;++i)
  {
    double b = gsl_vector_get(x_nlte_,i);
    lev_n_[i] = b*lev_lte_[i];
    // check that level populations aren't too negative or too large
    if (lev_n_[i] < -1.0e-5 || lev_n_[i] > 1.00001) {
      printf("problem with NLTE level pops\n");
      for (int j=0;j<n_levels_;++j) {
        printf("lev: %5d pop: %14.3e\n", j, lev_n_[j]);
      }
      exit(1);
    }

    if (lev_n_[i] < min_level_pop_)
      lev_n_[i] = min_level_pop_;
  }

  #endif

  // set the ionization fraction
  for (int i=0;i<n_ions_;++i)
    ion_frac_[i] = 0;
  for (int i=0;i<n_levels_;++i)
    ion_frac_[adata_->get_lev_ion(i)] += lev_n_[i];

  // this should return != 0 if failure really
  return 0;

}

//-------------------------------------------------------
// integrate up the radiation field over  lines
// to get the line J and over bound-free to get the
// photoionization rates
//-------------------------------------------------------
void AtomicSpecies::calculate_radiative_rates(std::vector<SedonaReal> J_nu)
{
  // zero out recombination/photoionization rates
  for (int j=0;j<n_levels_;++j)
  {
    lev_Pic_[j] = 0;
    lev_Rci_[j] = 0;
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
      int ic = adata_->get_lev_ic(j);
      if (ic == -1) continue;
      double chi = adata_->get_lev_Eion(j);
      if (E_ev < chi) continue;

      // photoionization term
      double sigma = adata_->get_lev_photo_cs(j,i);
      double Jterm = sigma*J/E_ergs;
      lev_Pic_[j] += Jterm*dnu;

      // recombination term
      lev_Rci_[j] += (sigma*fac1*nu*nu + Jterm)*exp(-1.0*(E_ev - chi)/pc::k_ev/gas_temp_)*dnu;

    }
  }

  // multiply by overall factors
  double lam_t = sqrt(pc::h*pc::h/(2*pc::pi*pc::m_e*pc::k*gas_temp_));
  double saha_fac = lam_t*lam_t*lam_t/2.0;
  for (int j=0;j<n_levels_;++j)
  {
    int ic = adata_->get_lev_ic(j);
    if (ic == -1) continue;
    double gl_o_gc = (1.0*adata_->get_lev_g(j))/(1.0*adata_->get_lev_g(ic));
    lev_Pic_[j] *= 4*pc::pi;
    lev_Rci_[j] *= 4*pc::pi*gl_o_gc*saha_fac;
  }

  // printing out for debug
  //for (int j=0;j<n_levels_;++j)
   // if (levels_[j].ic != -1)
    //  levels_[j].R_ci = 2.58e-13;
    //std::cout << j << " " << levels_[j].P_ic << " " << levels_[j].R_ci << "\n";

  // calculate line J's
  double x_max = 5.;
  double dx    = 0.05;

  for (int i=0;i<n_lines_;++i)
  {
    double nu0 = adata_->get_line_nu(i);
      double sum  = 0;
    double J0   = 0;

    double nu_d    = nu0*line_beta_dop_;
    double gamma   = adata_->get_line_A(i);
    double a_voigt = gamma/4/pc::pi/nu_d;

    for (double x=-1*x_max;x<=x_max;x+=dx)
    {
      double phi = voigt_profile_.getProfile(x,a_voigt);
      double n = nu0 + x*nu_d;
      double J1 = nu_grid_.value_at(n,J_nu)*phi;
      sum += 0.5*(J1 + J0)*dx;
      J0 = J1;
    }
    line_J_[i] = sum;
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
    int ll       = adata_->get_line_l(l);
    int lu       = adata_->get_line_u(l);

    // spontaneous dexcitation + stimulated emission
    double R_ul = adata_->get_line_Bul(l)*line_J_[l] + adata_->get_line_A(l);
    double R_lu = adata_->get_line_Blu(l)*line_J_[l];

    // check for transition between degenerate levels_
    if (adata_->get_line_nu(l) == 0)
      { R_ul = 0; R_lu = 0;}

    // add into rates
    rates_[ll][lu] += R_lu;
    rates_[lu][ll] += R_ul;

   // printf("RR %d %d %e\n",ll,lu,R_lu);
   // printf("RR %d %d %e\n",lu,ll,R_ul);
  }

  double norm = 0;
  for (int l=0;l<n_lines_;l++) norm   += adata_->get_line_f(l);

  for (int l=0;l<n_lines_;l++)
  {
    int ll       = adata_->get_line_l(l);
    int lu       = adata_->get_line_u(l);
    int gl       = adata_->get_lev_g(ll);
    int gu       = adata_->get_lev_g(lu);
    double El    = adata_->get_lev_E(ll);
    double Eu    = adata_->get_lev_E(lu);
    double f_lu  = adata_->get_line_f(l);


    // ------------------------------------------------
    // non-thermal (radioactive) bound-bound transitions
    // ------------------------------------------------
    double dE = (Eu - El)*pc::ev_to_ergs;
    double R_lu = 0; // e_gamma_/n_dens_/dE; //*(lines_[l].f_lu/norm);
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

    // Rutten section 3.2.5 (page 52) points out that these van Regemorter
    // rates are only valid for permitted dipole transitions, with f_lu in the
    // range 10^-3 to 1. For forbidden lines with smaller f, he says the
    // collisional transition rates "don't drop much below the values typical
    // of permitted lines." That's not a very precise statement, but we can
    // mock it up by not letting the f_lu factor drop below 10^-3

    double effective_f_lu = 0.;
    if (f_lu < 0.01) effective_f_lu = 0.01;
    else effective_f_lu = f_lu;

    if (use_collisions_nlte_)
    {
	    double C_up = 3.9*pow(zeta,-1.)*pow(gas_temp_,-1.5) / ezeta * ne * effective_f_lu;
      // be careful about possible overflow
      if (zeta > 700) C_up = 0;

	    double C_down = 3.9*pow(zeta,-1.)*pow(gas_temp_,-1.5) * ne * effective_f_lu * gl/gu;

	    rates_[ll][lu] += C_up;
      rates_[lu][ll] += C_down;
    }
  }


  // ------------------------------------------------
  // bound-free transitions
  // ------------------------------------------------
  for (int i=0;i<n_levels_;++i)
  {
    int ic = adata_->get_lev_ic(i);
    if (ic == -1) continue;

    // ionization potential
    int istage  = adata_->get_lev_ion(i);
    double chi  = adata_->get_ion_chi(istage)- adata_->get_lev_E(i);
    double zeta = chi/pc::k_ev/gas_temp_;

    // ------------------------------------------------
    // collisional ionization and recomombination rate
    // needs to be multiplied by number of electrons in outer shell
    // ------------------------------------------------

    if (use_collisions_nlte_)
    {
	    double C_ion = 2.7/zeta/zeta*pow(gas_temp_,-1.5)*exp(-zeta)*ne;
	    rates_[i][ic] += C_ion;

	    // collisional recombination rate
	    int gi = adata_->get_lev_g(i);
	    int gc = adata_->get_lev_g(ic);
	    double C_rec = 5.59080e-16/zeta/zeta*pow(gas_temp_,-3)*gi/gc*ne*ne;
	    rates_[ic][i] += C_rec;
      //printf("pc::cl: %d %d %e %e\n",i,ic,C_ion, C_rec);
    }

    // ------------------------------------------------
    // photoionization and radiative recombination
    // ------------------------------------------------

    // suppress recombinations to ground
    if (no_ground_recomb_)
    {
	     if (adata_->get_lev_E(i) == 0) lev_Rci_[i] = 0.;
    }
    // debug
    //lev_Rci_[i] = 0;
    rates_[ic][i] += lev_Rci_[i]*ne;
    rates_[i][ic] += lev_Pic_[i];

    //printf("pc::pi: %d %d %e %e\n",i,ic,lev_Pic_[i], lev_Rci_[i]*ne);

    // ------------------------------------------------
    // non-thermal (radioactive) bound-bound transitions
    // ------------------------------------------------
    double G_lc = e_gamma_ion_*adata_->get_nonthermal_ion_cross_section(istage,1e6);
    rates_[i][ic] += G_lc;
  }

   // multiply by rates by lte pop in level coming from
  // (becuase we will solve for depature coeffs)
  for (int i=0;i<n_levels_;++i)
      for (int j=0;j<n_levels_;++j)
        rates_[i][j] *= lev_lte_[i];

  // print out rates if you so like
  //printf("------- rates ----------\n");
  for (int i=0;i<n_levels_;++i)
    for (int j=0;j<n_levels_;++j)
      if (std::isnan(rates_[i][j]))
       printf("%5d %5d %14.5e\n",i,j,rates_[i][j]);
   //printf("\n");
   //for (int i=0;i<n_levels_;++i)
   // for (int j=0;j<n_levels_;++j)   {
   //     if (isnan(rates_[i][j])) std::cout << "NAN RATE\n";
   //     if (isinf(rates_[i][j])) std::cout << "INF RATE: " << i << " - " << j << "\n"; }
}



//-----------------------------------------------------------------
// Return the non-thermal ionization deposition
// given an electron with energy E (in eV)
// returns dE/dl (ergs/cm) = energy loss per cm traveled
//-----------------------------------------------------------------
double AtomicSpecies::get_nonthermal_ionization_dep(double E)
{
  double lambda = 0;
  for (int i=0;i<n_ions_;++i)
  {
    double Q = adata_->get_nonthermal_ion_cross_section(i,E);
    double chi = adata_->get_ion_chi(i)*pc::ev_to_ergs;
    lambda += Q*chi*n_dens_*ionization_fraction(i);
  }
  return lambda;

}


//-----------------------------------------------------------------
// calculate planck function in frequency units
//-----------------------------------------------------------------
double AtomicSpecies::blackbody_nu(double T, double nu)
{
  double zeta = pc::h*nu/pc::k/T;
  return 2.0*nu*nu*nu*pc::h/pc::c/pc::c/(exp(zeta)-1);
}
