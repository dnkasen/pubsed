#ifndef GAS_STATE_H
#define GAS_STATE_H
#include <string>

#include "AtomicSpecies.h"
#include "locate_array.h"
#include "sedona.h"
#include "hdf5.h"
#include "hdf5_hl.h"

class GasState
{

 private:

  double ne_brent_method(double,double,double,std::vector<real>);
  double charge_conservation(double,std::vector<real>);

  locate_array nu_grid_;
  int verbose_;
  int solve_error_;

  // global list of all levels of all atoms
  std::vector <int> globalLevelList_atom_;
  std::vector <int> globalLevelList_index_;


 public:

  std::vector<double>     mass_frac;  // vector of mass fractions
  std::vector<int>           elem_Z;  // vector of element atomic numbers
  std::vector<int>           elem_A;  // vector of element atomic weights
  std::vector<AtomicSpecies>  atoms;  // vector of atoms

  // pointer to stored atomic data
  AtomicData *atomic_data_;
  std::string atomdata_file_;    // filename for data

  double dens_;                  // total mass density (g cm^-3)
  double n_elec_;                // electron number density (cm^-3)
  double temp_;                  // Temperature (K)
  double time_;                  // Time since Explosion (days)
  double e_gamma;                // gamma-ray deposited energy
  double mu_I;                   // // mean atomic/ionic mass (not including free electrons). Dimensionless; needs to be multiplied by amu (~ m_p) to get units of grams
  int no_ground_recomb;          // suppress ground recombinations

  // flags for what opacities to use
  int use_electron_scattering_opacity;
  int use_line_expansion_opacity;
  int use_fuzz_expansion_opacity;
  int use_free_free_opacity;
  int use_bound_free_opacity;
  int use_bound_bound_opacity;
  int use_user_opacity_;
  int use_zone_specific_grey_opacity_;
  double line_velocity_width_;

  // calculate means
  double get_planck_mean(const std::vector<OpacityType>& x);
  double get_rosseland_mean(const std::vector<OpacityType>& x);
  double get_planck_mean
  (const std::vector<OpacityType> &abs, const std::vector<OpacityType>& scat);
  double get_rosseland_mean
  (const std::vector<OpacityType>& abs, const std::vector<OpacityType>& scat);

  std::vector <double> user_opacity_array_;

  // flags for nlte
  int use_nlte_;
  int use_collisions_nlte_;

  double bulk_grey_opacity_;    // bulk component of the grey opacity (cm^2/g), which is the same in every zone
  double total_grey_opacity_;   // total grey opacity (cm^2/g), which is the sum of the bulk component and the zone-specific component
  double epsilon_;
  std::vector<int> atom_zero_epsilon_;



  //***********************************************************
  // INITIALIZATION
  //***********************************************************

  //----------------------------------------------------------------
  // simple constructor
  //----------------------------------------------------------------
  GasState();

  //----------------------------------------------------------------
  // initialize the gas by specifying the atoms that will
  // compose it and freq. array, and name of atomic data file
  // inputs:
  // std::string atomfile: name of atom data file (in hdf5)
  // std::vector<int> e:  vector of atomic numbers
  // std::vector<int> A:  vector of atomic weights (in atomic units)
  // locate_array ng:  locate_array giving the freq. array
  //---------------------------------------------------------------
  void initialize
  (std::string af, std::vector<int> e, std::vector<int> A, locate_array ng);

  // initialize with atomic data class already created and passed
  // directly as pointer
  void initialize
  (AtomicData*, std::vector<int>, std::vector<int>, locate_array);


  //-----------------------------------------------------------
  // set which atoms to treat in nlte
  //-----------------------------------------------------------
  void set_atoms_in_nlte
  (std::vector<int> useatoms);


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
  void set_mass_fractions(std::vector<double>&);


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
  int solve_state(std::vector<real>&);
  int solve_state();


  // basic setting of properties
  void set_density(double d) {
    dens_ = d;
  }
  void set_temperature(double T) {
    temp_ = T;
  }
  void set_time(double t) {
    time_ = t;
  }

  //***********************************************************
  // RETURN BASIC DATA
  //***********************************************************

  //-----------------------------------------------------------
  // get essential parameters
  //-----------------------------------------------------------
  double get_density()             {return dens_;}
  double get_temperature()         {return temp_;}
  double get_mean_atomic_weight()  {return mu_I;}
  double get_electron_density()    {return n_elec_;}

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
    return atoms[i].n_ions_;
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
  double free_free_heating_rate(double, std::vector<real>);
  double free_free_cooling_rate(double);
  void bound_free_opacity (std::vector<double>&, std::vector<double>&);
  double bound_free_heating_rate (double, std::vector<real>);
  double bound_free_cooling_rate(double);
  double collisional_net_cooling_rate(double);
  void bound_bound_opacity(std::vector<double>&, std::vector<double>&);
  void bound_bound_opacity(int, std::vector<double>&, std::vector<double>&);
  void line_expansion_opacity(std::vector<double>&,std::vector<double>&);
  void fuzz_expansion_opacity(std::vector<double>&,std::vector<double>&);
  void get_line_opacities(std::vector<double>&);

  void get_user_defined_opacity
  (std::vector<double>&, std::vector<double>&,std::vector<double>&);

  void set_minimum_extinction(double d)
  {
    for (size_t i=0;i<atoms.size();++i) atoms[i].minimum_extinction_ = d;
  }

  void print_properties();
  void print();
  void print_memory_footprint();
  void write_levels(int iz);

};


#endif
