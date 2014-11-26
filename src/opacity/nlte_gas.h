#ifndef NLTE_GAS_H
#define NLTE_GAS_H

#include "nlte_atom.h"
#include "locate_array.h"
#include <string>

class nlte_gas
{
 
 private:

  double ne_brent_method(double,double,double,int);
  double charge_conservation(double,int);

  locate_array nu_grid;
  int verbose;

  std::vector<std::vector<int> > nlteLineList_atom_;
  std::vector<std::vector<int> > nlteLineList_ID_;
  std::vector<int> nlteLineList_nlines_;
  

 public:

  std::vector<double>     mass_frac;  // vector of mass fractions
  std::vector<int>           elem_Z;  // vector of element atomic numbers
  std::vector<int>           elem_A;  // vector of element atomic weights
  std::vector<nlte_atom>      atoms;  // vector of atoms

  double dens;                   // total mass density (g cm^-3)
  double ne;                     // electron number density (cm^-3)
  double temp;                   // Temperature (K)
  double time;                   // Time since Explosion (days) 
  double e_gamma;                // gamma-ray deposited energy
  double A_mu;                   // mean atomic weight of the gas
  int no_ground_recomb;          // suppress ground recombinations

  // flags for what opacities to use
  int use_electron_scattering_opacity;   
  int use_line_expansion_opacity;
  int use_fuzz_expansion_opacity;
  int use_free_free_opacity;
  int use_bound_free_opacity;
  
  double grey_opacity_;
  double epsilon_;
  
  //***********************************************************
  // INITIALIZATION
  //***********************************************************

  //----------------------------------------------------------------
  // simple constructor
  //----------------------------------------------------------------
  nlte_gas();
  
  //----------------------------------------------------------------
  // initialize the gas by specifying the atoms that will
  // compose it, along with datafile and freq. array
  // inputs:
  // std::string atomfile: name of atom data file (in hdf5)
  // std::vector<int> e:  vector of atomic numbers
  // std::vector<int> A:  vector of atomic weights (in atomic units)
  // locate_array ng:  locate_array giving the freq. array
  //---------------------------------------------------------------
  void initialize
    (std::string, std::vector<int>, std::vector<int>, locate_array);

  
  //-----------------------------------------------------------
  // read fuzz lines from a file
  // input:
  // std::string fuzzfile: name of hdf5 file with fuzz data
  //-----------------------------------------------------------
  int  read_fuzzfile(std::string fuzzfile);
  
  //-----------------------------------------------------------------
  // Set mass fractions of each element in the gas
  // this function will enforce that the mass fractions are
  // normalized (i.e., add up to 1)
  //
  // input: 
  // std::vector<double> x: vector of mass fractions of each element
  //-----------------------------------------------------------------
  void set_mass_fractions(std::vector<double>);
  
  
  //***********************************************************
  // BASIC FUNCTIONALITY
  //***********************************************************

  //-----------------------------------------------------------
  // Solve for the gas state (excitation/ionization)
  // the level populations will be stored internally for
  // further calculations
  // input:
  // int lte: 1 = do it in LTE, 0 = do it in NLTE
  //-----------------------------------------------------------  
  void solve_state(int lte);

  //-----------------------------------------------------------
  // return the ionization state, i.e., electron density
  // over total ion number density
  //-----------------------------------------------------------
  double get_ionization_state();

  //***********************************************************
  // OPACITIES AND EMISSIVITIES
  //***********************************************************
  void computeOpacity(std::vector<double>&, std::vector<double>&, 
		      std::vector<double>&);
  double electron_scattering_opacity();
  void line_expansion_opacity(std::vector<double>&);
  void fuzz_expansion_opacity(std::vector<double>&);

  void print();
};


#endif
