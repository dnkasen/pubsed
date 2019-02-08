#include <cstdlib>
#include <math.h>

#include "grid_1D_sphere.h"

void grid_1D_sphere::testCheckpointGrid(std::string fname) {
  testCheckpointGeneralGrid(fname);
  
  for (int rank = 0; rank < nproc; rank++) {
    if (rank == my_rank) {
      bool fail = false;
      if (v_inner_ != v_inner_new) {
        std::cerr << "issue at v_inner on rank " << rank << std::endl;
        fail = true;
      }
      if (vol.size() != vol_new.size()) {
        std::cerr << "vol arrays not same length on rank " << rank << std::endl;
        fail = true;
      }
      for (int i = 0; i < vol.size(); i++) {
        if (vol[i] != vol_new[i]) {
          std::cerr << "issue at vol, entry number " << i << " on rank " << rank << std::endl;
          fail = true;
        }
      }
      bool complain_about_locate_array = true;
      if (not r_out.is_equal(r_out_new, complain_about_locate_array)) {
        std::cerr << "issue at r_out on rank " << rank << std::endl;
        fail = true;
      }
      if (fail) {
        std::cerr << "1D grid restart failed on rank " << my_rank << std::endl;
        exit(3);
      }
      else std::cerr << "1D grid restart succeeded on rank " << my_rank << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
}

