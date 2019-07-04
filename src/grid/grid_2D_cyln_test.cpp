#include <cstdlib>
#include <math.h>

#include "grid_2D_cyln.h"

void grid_2D_cyln::testCheckpointGrid(std::string fname) {
  testCheckpointGeneralGrid(fname);
  for (int rank = 0; rank < nproc; rank++) {
    if (rank == my_rank) {
      bool fail = false;
      if (nx_ != nx_new_) {
        std::cerr << "issue at nx on rank " << rank << std::endl;
        fail = true;
      }
      if (nz_ != nz_new_) {
        std::cerr << "issue at nz on rank " << rank << std::endl;
        fail = true;
      }
      if (dx_ != dx_new_) {
        std::cerr << "issue at dx on rank " << rank << std::endl;
        fail = true;
      }
      if (dz_ != dz_new_) {
        std::cerr << "issue at dz on rank " << rank << std::endl;
        fail = true;
      }
      if (zcen_ != zcen_new_) {
        std::cerr << "issue at zcen on rank " << rank << std::endl;
        fail = true;
      }
      if (index_x_ != index_x_new_) {
        std::cerr << "issue at index_x on rank " << rank << std::endl;
        fail = true;
      }
      if (index_z_ != index_z_new_) {
        std::cerr << "issue at index_z on rank " << rank << std::endl;
        fail = true;
      }
      for (int i = 0; i < vol_.size(); i++) {
        if (vol_[i] != vol_new_[i]) {
          std::cerr << "issue at vol elem " << i <<std::endl;
        }
      }
      if (vol_.size() != vol_new_.size()) {
        std::cerr << "issue at vol size" << std::endl;
      }
      if (vol_ != vol_new_) {
        std::cerr << "issue at vol on rank " << rank << std::endl;
        fail = true;
      }

      if (fail) {
        std::cerr << "2D grid restart failed on rank " << my_rank << std::endl;
        exit(4);
      }
      else std::cerr << "2D grid restart succeeded on rank " << my_rank << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
}
