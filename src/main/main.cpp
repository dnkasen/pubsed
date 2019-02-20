#include "SedonaClass.h"

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

int main(int argc, char **argv)
{
  // initialize MPI parallelism
#ifdef MPI_PARALLEL
    MPI_Init( &argc, &argv );
#endif

  // This runs the entire calculation
  // see functions in SedonaClass.cpp
  SedonaClass sedona;
  sedona.run("param.lua");

#ifdef MPI_PARALLEL
    MPI_Finalize();
#endif

}
