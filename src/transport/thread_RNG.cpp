#pragma warning disable 161

#include "sedona.h"
#include "thread_RNG.h"
#include <ctime>

#include <mpi.h>

#ifdef _OPENMP
#include <omp.h>
#endif

//-----------------------------------------------------------------
// initialize the RNG system
//-----------------------------------------------------------------
// ASSUMES the number of threads remains constant so it only has to be initialized once
void thread_RNG::init()
{
  int my_mpiID = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_mpiID);

  // set up the stuff that creates the random number generators
  const gsl_rng_type* TypeR = gsl_rng_default;
  gsl_rng_env_setup();
  int nthreads;
  
  #pragma omp parallel
  #pragma omp single
  {
    #ifdef _OPENMP
    nthreads = omp_get_num_threads();
    #else
    nthreads = 1;
    #endif
  }
  if(my_mpiID==0) printf("# Using %d threads on each MPI rank.\n", nthreads);
    
  // assign a unique RNG to each thread
  generators.resize(nthreads);
  for(int i=0; i<nthreads; i++){
    if(i==0)
      gsl_rng_default_seed = (unsigned int)time(NULL);
    else
      gsl_rng_default_seed = gsl_rng_get(generators[0]);
    generators[i] = gsl_rng_alloc (TypeR);
  }
}



//-----------------------------------------------------------------
// return a uniformily distributed random number (thread safe)
//-----------------------------------------------------------------
double thread_RNG::uniform(){
  #ifdef _OPENMP
  const int my_ompID = omp_get_thread_num();
  #else 
  const int my_ompID = 0;
  #endif

  return gsl_rng_uniform(generators[my_ompID]);
}
