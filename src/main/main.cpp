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

  // get name of parameter file -- default is defined in sedona.h
  std::string param_file = DEFAULT_PARAM_FILE_NAME;
  if( argc > 1 ) param_file = std::string( argv[ 1 ] );

  // This runs the entire sedona calculation using the
  // passed parameter file name
  // See the functions in SedonaClass.cpp for how
  // it is all done
  SedonaClass sedona;
  sedona.run(param_file);

#ifdef MPI_PARALLEL
    MPI_Finalize();
#endif

}
