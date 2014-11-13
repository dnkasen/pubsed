#include <iostream>
#include "nlte_gas.h"
#include "locate_array.h"
#include "physical_constants.h"
#include <mpi.h>

using std::cout;
namespace pc = physical_constants;


int main(int argc, char **argv)
{
  void LTE_Ionization_Test(std::string);

  std::string atomdata = "/Users/kasen/codes/sedona6/data/atomdata.hdf5";

  // initialize MPI parallelism
  MPI_Init( &argc, &argv );
 
  LTE_Ionization_Test(atomdata);

}

void LTE_Ionization_Test(std::string atomdata)
{
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
  Z.push_back(1);
  A.push_back(1.0);

  // initialize the gas
  nlte_gas gas;
  gas.init(atomdata,Z,A,nu_grid);
  gas.print();

}
