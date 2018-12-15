#include <iostream>
#include <fstream>
#include "nlte_gas.h"
#include "locate_array.h"
#include "physical_constants.h"
#include <mpi.h>

using namespace std;
namespace pc = physical_constants;


int main(int argc, char **argv)
{
  void LTE_Ionization_Test(std::string);
  void Expansion_Opacity_Test(std::string atomdata);

  std::string atomdata = "/Users/kasen/codes/sedona6/data/cmfgen_atomdata.hdf5";

  // initialize MPI parallelism
  MPI_Init( &argc, &argv );

  //Expansion_Opacity_Test(atomdata);
  LTE_Ionization_Test(atomdata);

}


void Expansion_Opacity_Test(std::string atomdata)
{
  // density to use
  double dens = 1e-13;
  double temp = 10000;
  double texp = 20.0*3600.0*24.0;
  
  // set up frequency grid
  double lam_start = 100;
  double lam_stop  = 25000;
  int    n_pts     = 10000;
  double nu_1 = pc::c/(lam_stop*pc::angs_to_cm);
  double nu_2 = pc::c/(lam_start*pc::angs_to_cm);
  double dnu = (nu_2 - nu_1)/(1.0*n_pts);
  locate_array nu_grid;
  nu_grid.init(nu_1,nu_2,dnu);

  // set up composition and abundance
  std::vector<int> Z;
  std::vector<int> A;
  std::vector<double> mass_frac;
  Z.push_back(26);
  A.push_back(56);

  nlte_gas gas;
  gas.initialize(atomdata,Z,A,nu_grid);
  int nl = gas.read_fuzzfile("/Users/kasen/codes/sedona6/data/kurucz_cd23_fuzz.hdf5");
  std::cout << "Read and Stored " << nl << " fuzz lines\n";
  mass_frac.push_back(1.0);

  gas.set_mass_fractions(mass_frac);
  gas.time = texp;
  gas.temp = temp;
  gas.dens = dens;
  gas.solve_state(1);

  std::vector<double> line_opac;
  std::vector<double> fuzz_opac;
  fuzz_opac.resize(n_pts);
  line_opac.resize(n_pts);
  gas.line_expansion_opacity(line_opac);
  gas.fuzz_expansion_opacity(fuzz_opac);

  ofstream outfile("Fe_exp_opac.txt");
  for (int i=0;i<n_pts;i++)
  {
    double lam = pc::c/(nu_grid[i]*pc::angs_to_cm);
    outfile << lam << "\t" << line_opac[i]/dens << " ";
    outfile << fuzz_opac[i]/dens << "\n";
  }
  outfile.close();
  
}


void LTE_Ionization_Test(std::string atomdata)
{
  // density to use
  double dens = 1e-13;

  // set up frequency grid
  double lam_start = 100;
  double lam_stop  = 20000;
  int    n_pts     = 2000;
  double nu_1 = pc::c/(lam_stop*pc::angs_to_cm);
  double nu_2 = pc::c/(lam_start*pc::angs_to_cm);
  double dnu = (nu_2 - nu_1)/(1.0*n_pts);
  locate_array nu_grid;
  nu_grid.init(nu_1,nu_2,dnu);

  // set up composition and abundance
  std::vector<int> Z;
  std::vector<int> A;
  std::vector<double> mass_frac;
  Z.push_back(1);
  A.push_back(1);
  mass_frac.push_back(1.0);

  // do hydrogen
  {
    // initialize the gas
    nlte_gas gas;
    gas.initialize(atomdata,Z,A,nu_grid);
    mass_frac.push_back(1.0);
    gas.set_mass_fractions(mass_frac);
    
    ofstream outfile("H_ionization.txt");
    for (double temp=1000;temp < 20000;temp+=200)
    {
      gas.temp = temp;
      gas.dens = dens;
      gas.solve_state(1);
      outfile  << temp << "\t" << gas.get_ionization_state() << "\n";
    }
    outfile.close();
  }
  
  // do iron
  {
    // initialize the gas
    nlte_gas gas;
    Z[0] = 26;
    A[0] = 56; 
    gas.initialize(atomdata,Z,A,nu_grid);
    mass_frac[0] = 1.0;
    gas.set_mass_fractions(mass_frac);
    
    ofstream outfile("Fe_ionization.txt");
    for (double temp=1000;temp < 100000;temp+=1000)
    {
      gas.temp = temp;
      gas.dens = dens;
      gas.solve_state(1);
      outfile  << temp << "\t" << gas.get_ionization_state() << "\n";
    }
    outfile.close();
  }

  

}
