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

 public:

  std::vector<double>     mass_frac;  // vector of mass fractions
  std::vector<int>           elem_Z;  // vector of element atomic numbers
  std::vector<int>           elem_A;  // vector of element atomic weights
  std::vector<nlte_atom>      atoms;  // vector of atoms

  nlte_gas();
  void init(std::string, std::vector<int>, std::vector<int>, locate_array);
  int  read_fuzzfile(std::string fuzzfile);

  double dens;                   // total mass density (g cm^-3)
  double ne;                     // electron number density (cm^-3)
  double temp;                   // Temperature (K)
  double time;                   // Time since Explosion (days) 
  double e_gamma;                // gamma-ray deposited energy
  double A_mu;                   // mean atomic weight of the gas
  
  int no_ground_recomb;          // suppress ground recombinations


  void set_mass_fractions(std::vector<double> x)
  { 
    A_mu = 0;
    for (int i=0;i<mass_frac.size();i++) {
      mass_frac[i] = x[i]; 
      A_mu += mass_frac[i]*elem_A[i]; }
  }

  void solve_state(int lte);
  void print();

  // opacities and emissivities
  double electron_scattering_opacity();
  void line_expansion_opacity(std::vector<double>&);
  void fuzz_expansion_opacity(std::vector<double>&);
};


#endif
