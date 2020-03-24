#include <iostream>
#include <fstream>
#include <stdio.h>
#include "GasState.h"
#include "locate_array.h"
#include "physical_constants.h"
#include <mpi.h>

using namespace std;
namespace pc = physical_constants;

int main(int argc, char **argv)
{
  MPI_Init( &argc, &argv );

  //-------------------------------------------------------
  // parameters
  //-------------------------------------------------------
  std::vector<int>      elems_Z = {14,16,20};
  std::vector<int>      elems_A = {28,32,40};
  std::vector<double> massfrac  = {0.55,0.35,0.1};
  std::vector<int>  atoms_in_nlte = {14,16,20};
  std::string atomdata = "../../../data/cmfgen_levelcap100.hdf5";

  double time = 20*3600.0*24.0;
  double temp = 1e4;
  double dens = 1e-13;

  double nu_start = 1e13;
  double nu_stop  =1e17;
  double nu_delta = 0.001;   // logarthimc spacing
  double line_velocity_width = 3e8;

  //-------------------------------------------------------
  // Set frequency grid and radiation field
  //-------------------------------------------------------
  locate_array nu_grid;
  nu_grid.log_init(nu_start,nu_stop,nu_delta);
  int n_nu = nu_grid.size();

  std::vector<double> J_nu(n_nu);
  for (int i=0;i<n_nu;i++)
  {
    double W = 1.0;
    double nu = nu_grid.center(i);
    double zeta = pc::h*nu/pc::k/temp;
    double B_nu = W*2*pc::h*nu*nu*nu/pc::c/pc::c/(exp(zeta) - 1);
    J_nu[i] = B_nu;
  }

  //-------------------------------------------------------
  // Set up and run the gas class
  //-------------------------------------------------------

  // check for atomic data and load
  std::ifstream afile(atomdata);
  if (!afile) {
    std::cerr << "Can't open atom datafile " << atomdata << "; exiting" << std::endl;
    exit(1);   }
  afile.close();
  AtomicData* atomic_data = new AtomicData;
  atomic_data->initialize(atomdata,nu_grid);

  // initialize gas class
  GasState gas;
  gas.initialize(atomic_data,elems_Z,elems_A,nu_grid);
  gas.set_mass_fractions(massfrac);
  gas.set_atoms_in_nlte(atoms_in_nlte);

  // set gas parameters
  gas.use_nlte_ = 1;
  gas.line_velocity_width_ = line_velocity_width;
  gas.time_ = time;
  gas.temp_ = temp;
  gas.dens_ = dens;

  // solve the NLTE
  gas.solve_state(J_nu);

  // print out stuff
  gas.print();


}
