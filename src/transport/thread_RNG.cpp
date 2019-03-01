#include "hdf5.h"
#include "h5utils.h"
#include "sedona.h"
#include "thread_RNG.h"
#include <ctime>
#include <cstring>

#include <mpi.h>

#ifdef _OPENMP
#include <omp.h>
#endif

//-----------------------------------------------------------------
// initialize the RNG system
//-----------------------------------------------------------------
// ASSUMES the number of threads remains constant so it only has to be initialized once
void thread_RNG::init(bool fix_seed, unsigned long int fixed_seed_val)
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
  if (not fixed_seed)
    unsigned long int seed = (unsigned long int)time(NULL);
  else
    unsigned long int seed = fixed_seed_val;
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

void thread_RNG::writeCheckpointRNG(std::string fname) {
  int my_mpiID, mpi_nranks;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_mpiID);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_nranks);
  int ndim1 = 1;
  hsize_t one[1] = {1};
  size_t rng_size = gsl_rng_size(generators[0]);
  if (my_mpiID == 0) {
    createGroup(fname, "RNG");
    createDataset(fname, "RNG", "RNG_size", ndim1, one, H5T_NATIVE_INT);
    writeSimple(fname, "RNG", "RNG_size", &rng_size, H5T_NATIVE_INT);
  }
  // Create new opaque data type to hold RNG states
  // This does NOT take into account that GSL RNGs can have different types.
  // We're not storing those, and just assume that you want the defaults, and that
  // they're set the same in the version you're checkpointing from as they are
  // in the version you're restarting in. If you want to do something fancy,
  // there are GSL environment variables you can set.
  hid_t H5_rng = H5Tcreate(H5T_OPAQUE, rng_size);
  MPI_Barrier(MPI_COMM_WORLD);
  for (int i = 0; i < mpi_nranks; i++) {
    if (i == my_mpiID) {
      void* buffer;
      size_t rng_size = gsl_rng_size(generators[0]);
      hsize_t n_rngs = generators.size();
      buffer = malloc(n_rngs * rng_size);
      for (int i_gen = 0; i_gen < n_rngs; i_gen++) {
        void* rng_state = gsl_rng_state(generators[i_gen]);
        std::memcpy(buffer + i_gen * rng_size, rng_state, rng_size);
      }
      createDataset(fname, std::string("RNG"), std::to_string(my_mpiID), ndim1, &n_rngs, H5_rng);
      writeSimple(fname, std::string("RNG"), std::to_string(my_mpiID), buffer, H5_rng);
      free(buffer);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
}
