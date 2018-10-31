#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <vector>
#include "particle.h"
#include <iostream>



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

    int n = 100;
    char *filename = "out.dat";

    // array to write to file
    double *x = new double[n];
    for (int i=0;i<n;i++) x[i] = rank;

    int chunk = floor(n/nprocs);
    int offset = rank*chunk;

    std::cout << chunk << " " << offset << "\n";

//    int targetFile = 0;//
//    MPI_Comm fileComm;
//    MPI_Comm_split(MPI_COMM_WORLD, targetFile, rank, &fileComm);

    MPI_File outFile;
    MPI_File_open(
        MPI_COMM_WORLD, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY,
        MPI_INFO_NULL, &outFile);
    MPI_File_set_view( outFile, 0, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL ) ;

    MPI_File_write_at_all(
        outFile, offset,
        &x[offset], chunk, MPI_DOUBLE, MPI_STATUS_IGNORE);

    MPI_File_close(&outFile);
    MPI_Finalize();
}
