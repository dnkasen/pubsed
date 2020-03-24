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
#include "hdf5.h"
#include "hdf5_hl.h"

namespace pc = physical_constants;
using namespace std;


//----------------------------------------------------------------
// Constructor --
//----------------------------------------------------------------
AtomicSpecies::AtomicSpecies()
{
  e_gamma_            = 0;
  no_ground_recomb_   = 0;
  use_betas_          = 0;
  minimum_extinction_ = 0;
  use_nlte_           = 0;

  n_levels_            = 0;
  n_lines_             = 0;
  n_ions_              = 0;

  adata_  = NULL;

  // debug -- I hard coded this for now...
  min_level_pop_ = 1e-30;
}



//----------------------------------------------------------------
//  Destructor --
//----------------------------------------------------------------
AtomicSpecies::~AtomicSpecies()
{
}


//----------------------------------------------------------------
// Initialize
//----------------------------------------------------------------
int AtomicSpecies::initialize(int Z, AtomicData *ad)
{
  // set atomic number
  atomic_number = Z;

  // pointer to atomic data storage
  adata_ = ad->get_pointer_to_individual_atom(Z);

  // copy over some useful stuff to have locally
  nu_grid_.copy(ad->nu_grid_);
  n_ions_   = adata_->n_ions_;
  n_levels_ = adata_->n_levels_;
  n_lines_  = adata_->n_lines_;

  // allocate memory for state data
  ion_part_.resize(n_ions_);
  ion_frac_.resize(n_ions_);
  lev_n_.resize(n_levels_);
  lev_lte_.resize(n_levels_);
  lev_Pic_.resize(n_levels_);
  lev_Rci_.resize(n_levels_);
  line_J_.resize(n_lines_);

  return 0;
}

//----------------------------------------------------------------
// setup atom for solving things in non-LTE
//----------------------------------------------------------------
int AtomicSpecies::set_use_nlte()
{
  use_nlte_ = true;

  //allocate memory for the Eigen nLTE solver
  M_nlte_.resize(n_levels_,n_levels_);
  rates_.resize(n_levels_,n_levels_);
  b_nlte_.resize(n_levels_);
  x_nlte_.resize(n_levels_);

  //allocate memory for the rate matrix (this becomes M_nlte_)
  //rates_.reserve(n_levels_*n_levels_);

  //----------------------------------------------
  // allocate memory for arrays
  //----------------------------------------------
  // rates_ = new double*[n_levels_];
  // for (int i=0;i<n_levels_;++i) rates_[i] = new double[n_levels_];
  //
  // // matrix to solve
  // M_nlte_ = gsl_matrix_calloc(n_levels_,n_levels_);
  // gsl_matrix_set_zero(M_nlte_);
  //
  // // vector of level populations
  // x_nlte_ = gsl_vector_calloc(n_levels_);
  // gsl_vector_set_zero(x_nlte_);
  //
  // // right hand side vector
  // b_nlte_ = gsl_vector_calloc(n_levels_);
  // gsl_vector_set_zero(b_nlte_);
  //
  // // permuation vector, used internally for linear algebra solve
  // p_nlte_ = gsl_permutation_alloc(n_levels_);
  // gsl_permutation_init(p_nlte_);

  return 0;

}



void AtomicSpecies::print()
{

  cout << "--------------------- ions; n = " << n_ions_ << " ---------------------\n";
  cout << "# ion \t part \t frac \t chi (eV)\n";
  cout << "#---------------------------------------------------------------\n";


  for (int i=0;i<n_ions_;++i)
    cout << "#   "<<  adata_->ions_[i].stage << "\t" << ion_part_[i] << "\t" << ion_frac_[i]
	 << "\t" << adata_->ions_[i].chi << endl;


  cout << "\n";
  cout << "--------------------------------------------------------------------\n";
  cout << "--------------------levels; n = " << n_levels_ << "------------------------\n";
  cout << "# lev   ion     E_ex        g      pop          b_i       ion_to\n";
  cout << "#---------------------------------------------------------------\n";

  for (int i=0;i<n_levels_;++i)
  {
    printf("%5d %4d %12.3e %5d %12.3e %12.3e %5d\n",
	   i,adata_->levels_[i].ion,
	   adata_->levels_[i].E,adata_->levels_[i].g,lev_n_[i],
	   lev_n_[i]/lev_lte_[i],adata_->levels_[i].ic);
  }

  // printf("\n--- line data\n");
  //
  // for (int i=0;i<n_lines_;++i)
  // {
  //   printf("%8d %4d %4d %12.3e %12.3e %12.3e %12.3e %12.3e\n",
  //   	   i,lines_[i].ll,lines_[i].lu,lines_[i].nu,lines_[i].f_lu,
  //   	   lines_[i].A_ul,lines_[i].B_ul,lines_[i].B_lu);
  // }
  //
  // printf("\n--- line optical depths\n");
  //
  // for (int i=0;i<n_lines_;++i)
  // {
  //   int    ll = lines_[i].ll;
  //   double nl = levels_[ll].n;
  //
  //   printf("%8d %4d %4d %12.3e %12.3e %12.3e\n",
  //   	   i,lines_[i].ll,lines_[i].lu,lines_[i].nu,
	//    lines_[i].tau,nl);
  // }

}
