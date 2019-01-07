#include "hdf5.h"
#include "hdf5_hl.h"

#include <string.h>
#include <iostream>
#include <sstream>

#include "transport.h"
#include "physical_constants.h"

void transport::testCheckpointParticles() {
  for (int rank = 0; rank < MPI_nprocs; rank++) {
    if (rank == MPI_myID) {
      std::cerr << "rank " << rank << std::endl;
      for (int i = 0; i < particles.size(); i++) {
        particles[i].print();
      }
      std::cout << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  writeCheckpointParticles("test_particles.h5");

  MPI_Barrier(MPI_COMM_WORLD);

  readCheckpointParticles("test_particles.h5", true);

  for (int rank = 0; rank < MPI_nprocs; rank++) {
    if (rank == MPI_myID) {
      std::cerr << "rank " << rank << std::endl;
      for (int i = 0; i < particles_new.size(); i++) {
        particles_new[i].print();
        if (particles_new[i].type != particles[i].type) {
          std::cerr << "New particle type is different." << std::endl;
          exit(1);
        }
        if (particles_new[i].x[0] != particles[i].x[0]) {
          std::cerr << "New particle x0 is different." << std::endl;
          exit(1);
        }
        if (particles_new[i].x[1] != particles[i].x[1]) {
          std::cerr << "New particle x1 is different." << std::endl;
          exit(1);
        }
        if (particles_new[i].x[2] != particles[i].x[2]) {
          std::cerr << "New particle x2 is different." << std::endl;
          exit(1);
        }
        if (particles_new[i].D[0] != particles[i].D[0]) {
          std::cerr << "New particle D0 is different." << std::endl;
          exit(1);
        }
        if (particles_new[i].D[1] != particles[i].D[1]) {
          std::cerr << "New particle D1 is different." << std::endl;
          exit(1);
        }
        if (particles_new[i].D[2] != particles[i].D[2]) {
          std::cerr << "New particle D2 is different." << std::endl;
          exit(1);
        }
        if (particles_new[i].ind != particles[i].ind) {
          std::cerr << "New particle ind is different." << std::endl;
          exit(1);
        }
        if (particles_new[i].t != particles[i].t) {
          std::cerr << "New particle t is different." << std::endl;
          exit(1);
        }
        if (particles_new[i].e != particles[i].e) {
          std::cerr << "New particle e is different." << std::endl;
          exit(1);
        }
        if (particles_new[i].nu != particles[i].nu) {
          std::cerr << "New particle nu is different." << std::endl;
          exit(1);
        }
        if (particles_new[i].gamma != particles[i].gamma) {
          std::cerr << "New particle gamma is different." << std::endl;
          exit(1);
        }
        if (particles_new[i].dshift != particles[i].dshift) {
          std::cerr << "New particle dshift is different." << std::endl;
          exit(1);
        }
        if (particles_new[i].dvds != particles[i].dvds) {
          std::cerr << "New particle dvds is different." << std::endl;
          exit(1);
        }
        if (particles_new[i].fate != particles[i].fate) {
          std::cerr << "New particle fate is different." << std::endl;
          exit(1);
        }
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  std::cerr << "Partices written out and back in correctly on rank " << MPI_myID << std::endl;
}

void transport::testCheckpointSpectrum() {
  optical_spectrum.writeCheckpointSpectrum("spectrum_test.h5", "optical_spectrum");
  
  MPI_Barrier(MPI_COMM_WORLD);

  optical_spectrum_new.readCheckpointSpectrum("spectrum_test.h5", "optical_spectrum");

  bool complain = true;
  if (not optical_spectrum.is_equal(optical_spectrum_new, complain)) {
    std::cerr << "optical spectra not equal" << std::endl;
    exit(3);
  }
}
