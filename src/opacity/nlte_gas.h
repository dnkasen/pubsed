#ifndef NLTE_GAS_H
#define NLTE_GAS_H
#include <string>

#include "nlte_atom.h"
#include "locate_array.h"
#include "sedona.h"

class nlte_gas
{
 
 private:

  double ne_brent_method(double,double,double,std::vector<real>);
  double charge_conservation(double,std::vector<real>);

  locate_array nu_grid;
  int verbose;
  int solve_error_;
  std::string atomfile_;

  // global list of all levels of all atoms
  std::vector <int> globalLevelList_atom_;
  std::vector <int> globalLevelList_index_;


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
  int use_bound_bound_opacity;   
  int use_user_opacity_;
  double line_velocity_width_;

  std::vector <double> user_opacity_array_;
  
  // flags for nlte
  int use_nlte_;
  std::vector<int> atoms_in_nlte_;

  double grey_opacity_;
  double epsilon_;
  std::vector<int> atom_zero_epsilon_;



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
  // int: use_nlte: whether to use nlte
  //---------------------------------------------------------------
  void initialize
  (std::string, std::vector<int>, std::vector<int>, locate_array, int);

  //---------------------------------------------------------------
  // initialize: overload above with default is use_nlte = 0
  //---------------------------------------------------------------
  void initialize
  (std::string af, std::vector<int> e, std::vector<int> A, locate_array ng)
  {
    initialize(af,e,A,ng,0);
  }

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
  // input: a vector of the radiation field J_nu
  // output: a flag for status
  //    == 0 for no error
  //    == 1 root not bracketed in electron density solve
  //    == 2 maximum iterations reached in n_e solve
  //-----------------------------------------------------------  
  int solve_state(std::vector<real>);
  int solve_state();


  //***********************************************************
  // RETURN BASIC DATA
  //***********************************************************

  //-----------------------------------------------------------
  // get essential parameters 
  //-----------------------------------------------------------
  double get_density()             {return dens;} 
  double get_temperature()         {return temp;}
  double get_mean_atomic_weight()  {return A_mu;}
  double get_electron_density()    {return ne;}

  //-----------------------------------------------------------
  // return the the number of atoms, or 
  // number of ions in atom i
  //-----------------------------------------------------------
  int get_number_atoms()
  {
    return (int)(atoms.size());
  }
  int get_number_ions(int i)
  {
    return atoms[i].n_ions;
  }

  //-----------------------------------------------------------
  // return the ionization state, i.e., electron density
  // over total ion number density
  //-----------------------------------------------------------
  double get_ionization_state();

  //-----------------------------------------------------------
  // return the fraction of atoms of index i that are in
  // ionization state j.  
  //-----------------------------------------------------------
  double get_ionization_fraction(int i, int j);

  //-----------------------------------------------------------
  // return the fraction of atoms of index i that are in
  // ionization state j.  
  //-----------------------------------------------------------
  double get_partition_function(int i, int j)
  {
    return atoms[i].partition(j);
  }

  //-----------------------------------------------------------
  // return the fraction of atoms of index i that are in
  // level state j.  
  //-----------------------------------------------------------
  double get_level_fraction(int i, int j);
  double get_level_departure(int i, int j);


  //***********************************************************
  // OPACITIES AND EMISSIVITIES
  //***********************************************************
  void computeOpacity(std::vector<OpacityType>&, std::vector<OpacityType>&, 
		      std::vector<OpacityType>&);
  double electron_scattering_opacity();
  void free_free_opacity  (std::vector<double>&, std::vector<double>&);
  void bound_free_opacity (std::vector<double>&, std::vector<double>&);
  void bound_bound_opacity(std::vector<double>&, std::vector<double>&);
  void bound_bound_opacity(int, std::vector<double>&, std::vector<double>&);
  void line_expansion_opacity(std::vector<double>&);
  void fuzz_expansion_opacity(std::vector<double>&, std::vector<double>&);
  void get_line_opacities(std::vector<double>&);
 
  void get_user_defined_opacity
  (std::vector<double>&, std::vector<double>&,std::vector<double>&);

  void set_minimum_extinction(double d)
  {
    for (size_t i=0;i<atoms.size();++i) atoms[i].minimum_extinction_ = d;
  }

  void print_properties();
  void print();
};


#endif
