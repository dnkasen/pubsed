#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <vector>
#include "particle.h"
#include <iostream>
#include "hdf5.h"
#include "hdf5_hl.h"
#include "h5utils.h"


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


int main(){
    MPI_Init(NULL, NULL);


    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    int n = 4;
    char *filename = "out.dat";

    // array to write to file
    double *x = new double[n];
    for (int i=0;i<n;i++) x[i] = rank;


    std::cerr << H5T_NATIVE_DOUBLE << std::endl;
    std::cerr << "set x" << std::endl;
    int chunk = floor(n/nprocs);
    int offset = rank*chunk;

    double *y = new double[chunk];
    for (int i = 0; i < chunk; i++) y[i] = rank * 2.2;

    std::cout << chunk << " " << offset << "\n";

//    int targetFile = 0;//
//    MPI_Comm fileComm;
//    MPI_Comm_split(MPI_COMM_WORLD, targetFile, rank, &fileComm);

    if (rank == 0) {
      hsize_t dims_g[1]={n};
      createFile("test.h5");
      createGroup("test.h5", "agroup");
      createDataset("test.h5", "agroup", "adataset", 1, dims_g, H5T_NATIVE_DOUBLE);
      writeSimple("test.h5", "agroup", "adataset", x, H5T_NATIVE_DOUBLE);
      createFile("test_par.h5");
      createGroup("test_par.h5", "agroup");
      createDataset("test_par.h5", "agroup", "adataset", 1, dims_g, H5T_NATIVE_DOUBLE);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    for (int i = 0; i < nprocs; i++) {
      if (i == rank) {
        std::cerr << "rank " << rank << std::endl;
        std::cerr << "offset " << offset << std::endl;
        std::cerr << "chunk " << chunk << std::endl;
        std::cerr << "glochunk " << nprocs * chunk << std::endl;
        std::cerr << "y " << y[offset] << std::endl;
        int starts[1]={offset};
        int loc_size[1]={chunk};
        int glo_size[1]={nprocs * chunk};
        writePatch("test_par.h5", "agroup", "adataset", y, H5T_NATIVE_DOUBLE, 1, starts, loc_size, glo_size);
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }

    double* xnew = new double[n];
    double* ynew = new double[chunk];
    // Now read files back in
    for (int i = 0; i < nprocs; i++) {
      if (i == rank) {
        int starts[1]={offset};
        int loc_size[1]={chunk};
        int glo_size[1]={nprocs * chunk};
        readPatch("test_par.h5", "agroup", "adataset", ynew, H5T_NATIVE_DOUBLE, 1, starts, loc_size, glo_size);
        for (int j=0; j < chunk; j++) {
          std::cerr << "y " << j << " " << ynew[j] << std::endl;
        }
        readSimple("test.h5", "agroup", "adataset", xnew, H5T_NATIVE_DOUBLE);
        for (int j=0; j < n; j++) {
          std::cerr << "x " << j << " " << xnew[j] << std::endl;
        }
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
        

    MPI_Finalize();

}
