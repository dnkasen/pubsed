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


struct fuzz_line_structure
{
  int n_lines;
  std::vector<double> nu;
  std::vector<double> El;
  std::vector<double> gf;
  std::vector<int>   ion;
  std::vector<int>   bin;
};


struct AtomicIon
{
  int stage;          // ionization stage (0 = neutral, 1 = +, etc..)
  int ground;         // index of ground state level
  double chi;         // ionization energy, in eV
  double part;        // partition function
  double frac;        // fractional abundance
};

struct AtomicLine
{
  int lu,ll;           // index of upper/lower level
  double nu;           // rest frequency (Hz)
  double f_lu;         // oscillator strength
  double A_ul;         // Einstein A coeficient
  double B_ul;         // Einstein B coeficient
  double B_lu;         // Einstein B coeficient
  double J;            // radiation field in line
  double tau;          // sobolev optical depth
  double etau;         // exponential of soblev optical depth
  double beta;         // sobolev escape probability
  int    bin;          // index of the nu grid bin
};

struct AtomicLevel
{
  int   globalID;       // global id
  int  ion;             // ionization state (0 = neutral)
  int   ic;             // index of level we ionize to (-1 if none)
  int    g;             // statistical weight
  double E;             // excitation energy above ground (in eV)
  double E_ion;         // energy to ionize (in eV)
  double n;             // level population fraction
  double n_lte;         // lte level population
  double b;             // nlte departure coefficient

  //  ivector photo         // photoionization cross-section vector
  double  P_ic;            // rate of photoionization from this level
  double  R_ci;            // radiative recombination rate to this level

  // photoionization cross-section as a function of wavelength
  xy_array s_photo;
  // recombination coefficient as a function of temperature
  xy_array a_rec;

};


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

  // atomic data structs
  int n_ions_;             // Number of ionic stages considered
  int n_levels_;           // number of energy levels
  int n_lines_;            // number of line transitions
  AtomicLevel *levels_;       // array of level data
  AtomicLine  *lines_;        // array of line data
  AtomicIon    *ions_;        // array of ion data

  fuzz_line_structure fuzz_lines; // vector of fuzz lines

  // Constructor and Init
  AtomicSpecies();
  int initialize(std::string, int, locate_array, int&);
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
  void   bound_free_opacity_for_heating (std::vector<double>&, double,double);
  void   bound_bound_opacity(std::vector<double>&, std::vector<double>&);
  void   bound_free_opacity_for_cooling (std::vector<double>&, double,double);
  double collisional_net_cooling_rate(double, double);


  // returns
  int get_n_fuzz_lines()
  {
    return fuzz_lines.n_lines;
  }

  double partition(int ion)
  {
    for (int i = 0;i < n_ions_; i++)
      if (ion == ions_[i].stage) return ions_[i].part;
    return -1;
  }

  double ionization_fraction(int ion)
  {
    for (int i=0;i<n_ions_; i++)
      if (ion == ions_[i].stage) return ions_[i].frac;
    return 0;
 }

  double level_fraction(int lev)
  {
    if (lev >= n_levels_) return 0;
    return levels_[lev].n;
  }
  double level_depature(int lev)
  {
    if (lev >= n_levels_) return 0;
    return levels_[lev].b;
  }


 };

#endif
