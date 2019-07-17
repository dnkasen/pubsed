#ifndef _ATOMIC_SPECIES_H
#define _ATOMiC_SPECIES_H 1

#include <string>
#include <vector>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>
#include "xy_array.h"
#include "locate_array.h"
#include "sedona.h"
#include "VoigtProfile.h"
#include "AtomicData.h"


class AtomicSpecies
{

private:

  // matrices, vectors to solve NLTE
  double     **rates_;
  gsl_matrix *M_nlte_;
  gsl_vector *b_nlte_;
  gsl_vector *x_nlte_;
  gsl_permutation *p_nlte_;


  // frequency bin array
  locate_array nu_grid_;

  // pointer to atomic data holder
  IndividualAtomData *adata_;

  // Voigt profile class
  VoigtProfile voigt_profile_;

  double blackbody_nu(double T, double nu);
  double Calculate_Milne(int lev, double temp);
  void   set_rates(double ne);

public:



  int atomic_number;      // Atomic number of atom

  double n_dens_;          // number density of this atom (cm^-3)
  double e_gamma_;         // radioactive energy deposited (ergs/sec/cm^3)
  double gas_temp_;        // temperature of gas

  double min_level_pop_;        // the minimum level population allowed
  double minimum_extinction_;   // minimum alpha = 1/mfp to calculate
  double line_beta_dop_;        // doppler width of lines = v/c
  int use_betas_;               // include escape probabilites in nlte
  int no_ground_recomb_;        // suppress recombinations to ground
  bool use_nlte_;               // treat this atom in nlte or not
  int use_collisions_nlte_;    // use collisional transitions in NLTE solve and heating/cooling

  // atomic state values
  int n_ions_;             // Number of ionic stages considered
  int n_levels_;           // number of energy levels
  int n_lines_;            // number of line transitions

  std::vector<double> ion_part_;    // ion partition function
  std::vector<double> ion_frac_;    // ion fraction
  std::vector<double> lev_n_;       // level number population
  std::vector<double> lev_lte_;     // level LTE number population
  std::vector<double> lev_Pic_;     // photoionization rate from level
  std::vector<double> lev_Rci_;     // recombination rate to level
  std::vector<double> line_J_;      // line radiation field

  // Constructor and Init
  AtomicSpecies();
  int initialize(int, AtomicData*);
  int set_use_nlte();
  int read_fuzzfile(std::string);

  // Destructor
  ~AtomicSpecies();

  // solve state
  void calculate_radiative_rates(std::vector<real> J_nu);
  int  solve_state(double ne);
  int  solve_lte (double ne);
  int  solve_nlte(double ne);
  void print();

  // sobolev
  double get_ion_frac();
  double compute_sobolev_tau(int i, double time);
  void   compute_sobolev_taus(double time);

  // opacities and heating/cooling rates
  void   bound_free_opacity (std::vector<double>&, std::vector<double>&, double);
  void   bound_free_opacity_general (std::vector<double>&, std::vector<double>&, double, double, int);
  void   bound_free_opacity_for_heating (std::vector<double>&, double,double);
  void   bound_free_opacity_for_cooling (std::vector<double>&, double,double);
  double collisional_net_cooling_rate(double, double);
  void   bound_bound_opacity(std::vector<double>&, std::vector<double>&);
  void   line_expansion_opacity(std::vector<double>&,double);
  void   fuzzline_expansion_opacity(std::vector<double>& opac, double time);

  // returns
  int get_n_fuzz_lines()
  {
    return adata_->get_n_fuzz_lines();
  }

  double partition(int ion)
  {
    for (int i = 0;i < n_ions_; i++)
      if (ion == adata_->ions_[i].stage) return ion_part_[i];
    return -1;
  }

  double ionization_fraction(int ion)
  {
    for (int i=0;i<n_ions_; i++)
      if (ion == adata_->ions_[i].stage) return ion_frac_[i];
    return 0;
  }

  double get_net_ion_fraction()
  {
    double x = 0;
    for (int i=0;i<n_levels_;++i)
      x += lev_n_[i]*adata_->get_lev_ion(i);
    return x;
  }

  double level_fraction(int lev)
  {
    if (lev >= n_levels_) return 0;
    return lev_n_[lev];
  }
  double level_depature(int lev)
  {
    if (lev >= n_levels_) return 0;
    return lev_n_[lev]/lev_lte_[lev];
  }


 };

#endif
