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

  unsigned long int seed;
  // set generator for rank 0 thread 0
  if (not fix_seed)
    seed = (unsigned long int)time(NULL);
  else
    seed = fixed_seed_val;
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
    createDataset(fname, "RNG", "n_ranks", ndim1, one, H5T_NATIVE_INT);
    writeSimple(fname, "RNG", "n_ranks", &mpi_nranks, H5T_NATIVE_INT);
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
      uint8_t* buffer;
      size_t rng_size = gsl_rng_size(generators[0]);
      hsize_t n_rngs = generators.size();
      buffer = new uint8_t[n_rngs * rng_size];
      for (int i_gen = 0; i_gen < n_rngs; i_gen++) {
        uint8_t* rng_state = (uint8_t*) gsl_rng_state(generators[i_gen]);
        std::memcpy(&buffer[i_gen * rng_size], rng_state, rng_size);
      }
      createDataset(fname, std::string("RNG"), std::to_string(my_mpiID), ndim1, &n_rngs, H5_rng);
      writeSimple(fname, std::string("RNG"), std::to_string(my_mpiID), buffer, H5_rng);
      delete [] buffer;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
}

// Returns a status for whether or not this was successful
int thread_RNG::readCheckpointRNG(std::string fname) {
  int my_mpiID, mpi_nranks;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_mpiID);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_nranks);
  int n_ranks_chk;
  if (my_mpiID == 0) {
    readSimple(fname, "RNG", "n_ranks", &n_ranks_chk, H5T_NATIVE_INT);
  }
  MPI_Bcast(&n_ranks_chk, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (n_ranks_chk > mpi_nranks) {
    if (my_mpiID == 0) {
      std::cerr << "WARNING: More ranks checkpointed RNG state than there are currently running. " <<
        mpi_nranks << " now, were " << n_ranks_chk << ". Only first " << mpi_nranks <<
        " of RNGs are being read in." << std::endl;
    }
  }
  else if (n_ranks_chk < mpi_nranks) {
    if (my_mpiID == 0) {
      std::cerr << "Fewer ranks checkpointed RNG state than there are currently running. " <<
        mpi_nranks << " now, were " << n_ranks_chk << ". Generating new seeds/states." << std::endl;
    }
    return 1;
  }

  for (int i = 0; i < mpi_nranks; i++) {
    if (i == my_mpiID) {
      int rng_size;
      hsize_t num_rngs;
      readSimple(fname, std::string("RNG"), "RNG_size", &rng_size, H5T_NATIVE_INT);
      hid_t H5_rng = H5Tcreate(H5T_OPAQUE, rng_size);
      getH5dims(fname, std::string("RNG"), std::to_string(my_mpiID), &num_rngs);
      uint8_t* buffer = new uint8_t[num_rngs * rng_size];
      readSimple(fname, std::string("RNG"), std::to_string(my_mpiID), buffer, H5_rng);

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
      if (nthreads != num_rngs) {
        std::cerr << "The number of RNGs in the checkpoint file for rank " << my_mpiID <<
          " does not equal the number of threads on the rank. " << num_rngs << " vs " << 
          nthreads <<std::endl << ". Generating new seeds/states." << std::endl;
        return 1;
      }
      generators.resize(nthreads);
      for (int j = 0; j < nthreads; j++) {
        generators[j] = gsl_rng_alloc(TypeR);
        std::memcpy(gsl_rng_state(generators[j]), &buffer[j * rng_size], rng_size);
      }
      delete [] buffer;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  return 0;
}
