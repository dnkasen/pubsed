#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <vector>
#include "transport.h"
#include "sedona.h"
#include "ParameterReader.h"
#include "particle.h"
#include <iostream>
#include "hdf5.h"
#include "hdf5_hl.h"
#include "h5utils.h"
#include "grid_general.h"
#include "grid_1D_sphere.h"
#include "grid_2D_cyln.h"
#include "grid_3D_cart.h"


// little test C++ code to test doing parallel mpi write to binary file
// Is an approach for writing checkpoint files, since particle data
// must (or probably should) be written in parallel


    // A possible structure of checkpoint files
	//
	// !!! grid data (nz = # of zones; nelem = # of elems)
    // time
	// Z[nz]   # atomic numbers
	// A[nz]   # atomic weights
    // density[nz]
    // temp[nz]
    // comp[nz*nelem]
	// !!! particle data  (np = # of particles)
    // x[np]
    // y[np]
    // z[np]
    // Dx[np]
    // Dy[np]
    // Dz[np]
    // e[np]
    // t[np]
    // nu[np]
    // type[np] 


int main(int argc, char **argv){
    MPI_Init(NULL, NULL);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    const int verbose = (rank == 0);

    std::string param_file = "param.lua";
    if( argc > 1 ) param_file = std::string( argv[ 1 ] );
    ParameterReader params(param_file,verbose);

    std::string test_fname = "test_check.h5";
    createFile(test_fname);
    grid_general *grid;

    // read the grid type
    string grid_type = params.getScalar<string>("grid_type");

    // create a grid of the appropriate type
    if      (grid_type == "grid_1D_sphere") grid = new grid_1D_sphere;
    else if (grid_type == "grid_2D_cyln"  ) grid = new grid_2D_cyln;
    else if (grid_type == "grid_3D_cart"  ) grid = new grid_3D_cart;
    else  {
      if(verbose) std::cerr << "# ERROR: the grid type is not implemented" << std::endl;
      exit(3);   }

    // initialize the grid (including reading the model file)
    std::cerr << "initializing grid" << std::endl;
    grid->init(&params);

    std::cerr << "initializing mcarlo" << std::endl;
    transport mcarlo;
    string transport_type = params.getScalar<string>("transport_module");
    int use_transport = 0;
    if (transport_type != "") use_transport = 1;
    if (use_transport) mcarlo.init(&params, grid);

    std::cerr << "testing particles" << std::endl;
    mcarlo.testCheckpointParticles(test_fname);

    for (int i = 0; i < grid->elems_A.size();i++) {
      std::cerr << grid->elems_A[i] << " " << grid->elems_Z[i] << std::endl;
    }

    std::cerr << "testing grid" << std::endl;
    grid->testCheckpointGrid(test_fname);
    std::cerr << "testing zones" << std::endl;
    grid->testCheckpointZones(test_fname);

    mcarlo.testCheckpointSpectrum(test_fname);
    delete grid;

    MPI_Finalize();

}
