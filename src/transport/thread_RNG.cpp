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
  int my_mpiID, mpi_nranks;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_mpiID);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_nranks);

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
  generators.resize(nthreads);

  // set generator for rank 0 thread 0
  unsigned long int seed = (unsigned long int)time(NULL);
  generators[0] = gsl_rng_alloc(TypeR);
  gsl_rng_set(generators[0], seed);

  // assign unique RNG to each rank.
  for(int i=0; i<my_mpiID; i++) seed = gsl_rng_get(generators[0]);
  gsl_rng_set(generators[0], seed);

  // assign a unique RNG to each thread
  for(int thread=1; thread<nthreads; thread++){
    seed = gsl_rng_get(generators[0]);
    generators[thread] = gsl_rng_alloc(TypeR);
    gsl_rng_set(generators[thread], seed);
  }
}



//-----------------------------------------------------------------
// return a uniformily distributed random number (thread safe)
//-----------------------------------------------------------------
double thread_RNG::uniform()
{
#ifdef _OPENMP
  const int my_ompID = omp_get_thread_num();
#else
  const int my_ompID = 0;
#endif

  return gsl_rng_uniform(generators[my_ompID]);
}

thread_RNG::~thread_RNG() {
  for (auto i = generators.begin(); i != generators.end(); ++i) {
    gsl_rng_free(*i);
  }
}
