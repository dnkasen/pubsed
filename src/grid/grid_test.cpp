#include <cstdlib>
#include <math.h>

#include "grid_general.h"
#include "physical_constants.h"
namespace pc = physical_constants;

void grid_general::testCheckpointZones() {
  std::cerr << "starting to test zones" << std::endl;

  writeCheckpointZones("test_zones.h5");

  std::cerr << "wrote checkpoint zones" << std::endl;

  MPI_Barrier(MPI_COMM_WORLD);

  readCheckpointZones("test_zones.h5");

  std::cerr << "read checkpoint zones" << std::endl;

  for (int rank = 0; rank < nproc; rank++) {
    std::cerr << "rank " << my_rank << std::endl;
    if (rank == my_rank) {
      bool fail = false;
      for (int i = 0; i < n_zones; i++) {
        if ((z_new[i].v[0] != z[i].v[0]) || (z_new[i].v[1] != z[i].v[1]) || (z_new[i].v[2] != z[i].v[2])) {
          std::cerr << "issue at v on zone " << i << " on rank " << rank << std::endl;
          fail = true;
        }
        if (z_new[i].rho != z[i].rho) {
          std::cerr << "issue at rho on zone " << i << " on rank " << rank << std::endl;
          fail = true;
        }
        if (z_new[i].cs != z[i].cs) {
          std::cerr << "issue at cs on zone " << i << " on rank " << rank << std::endl;
          fail = true;
        }
        if (z_new[i].p_gas != z[i].p_gas) {
          std::cerr << "issue at p_gas on zone " << i << " on rank " << rank << std::endl;
          fail = true;
        }
        if (z_new[i].T_gas != z[i].T_gas) {
          std::cerr << "issue at T_gas on zone " << i << " on rank " << rank << std::endl;
          fail = true;
        }
        if (z_new[i].mu != z[i].mu) {
          std::cerr << "issue at mu on zone " << i << " on rank " << rank << std::endl;
          fail = true;
        }
        if (z_new[i].e_rad != z[i].e_rad) {
          std::cerr << "issue at e_rad on zone " << i << " on rank " << rank << std::endl;
          std::cerr << z_new[i].e_rad << " " << z[i].e_rad << std::endl;
          fail = true;
        }
        if (z_new[i].e_abs != z[i].e_abs) {
          std::cerr << "issue at e_abs on zone " << i << " on rank " << rank << std::endl;
          fail = true;
        }
        if (z_new[i].fx_rad != z[i].fx_rad) {
          std::cerr << "issue at fx on zone " << i << " on rank " << rank << std::endl;
          fail = true;
        }
        if (z_new[i].fy_rad != z[i].fy_rad) {
          std::cerr << "issue at fy on zone " << i << " on rank " << rank << std::endl;
          fail = true;
        }
        if (z_new[i].fz_rad != z[i].fz_rad) {
          std::cerr << "issue at fz on zone " << i << " on rank " << rank << std::endl;
          fail = true;
        }
        if (z_new[i].eps_imc != z[i].eps_imc) {
          std::cerr << "issue at eps_imc on zone " << i << " on rank " << rank << std::endl;
          fail = true;
        }
        if (z_new[i].L_thermal != z[i].L_thermal) {
          std::cerr << "issue at L_thermal on zone " << i << " on rank " << rank << std::endl;
          fail = true;
        }
        if (z_new[i].L_radio_emit != z[i].L_radio_emit) {
          std::cerr << "issue at L_radio_emit on zone " << i << " on rank " << rank << std::endl;
          fail = true;
        }
        if (z_new[i].L_radio_dep != z[i].L_radio_dep) {
          std::cerr << "issue at L_radio_dep on zone " << i << " on rank " << rank << std::endl;
          fail = true;
        }
        if (z_new[i].G1 != z[i].G1) {
          std::cerr << "issue at g1 on zone " << i << " on rank " << rank << std::endl;
          fail = true;
        }
        if (z_new[i].G2 != z[i].G2) {
          std::cerr << "issue at g2 on zone " << i << " on rank " << rank << std::endl;
          fail = true;
        }
        if (z_new[i].G3 != z[i].G3) {
          std::cerr << "issue at g3 on zone " << i << " on rank " << rank << std::endl;
          fail = true;
        }
        if (z_new[i].P11 != z[i].P11) {
          std::cerr << "issue at p11 on zone " << i << " on rank " << rank << std::endl;
          fail = true;
        }
        if (z_new[i].P12 != z[i].P12) {
          std::cerr << "issue at p12 on zone " << i << " on rank " << rank << std::endl;
          fail = true;
        }
        if (z_new[i].P13 != z[i].P13) {
          std::cerr << "issue at p13 on zone " << i << " on rank " << rank << std::endl;
          fail = true;
        }
        if (z_new[i].P22 != z[i].P22) {
          std::cerr << "issue at p22 on zone " << i << " on rank " << rank << std::endl;
          fail = true;
        }
        if (z_new[i].P23 != z[i].P23) {
          std::cerr << "issue at p23 on zone " << i << " on rank " << rank << std::endl;
          fail = true;
        }
        if (z_new[i].P33 != z[i].P33) {
          std::cerr << "issue at P33 zone " << i << " on rank " << rank << std::endl;
          fail = true;
        }
        for (int j = 0; j < n_elems; j++) {
          if (z_new[i].X_gas[j] != z[i].X_gas[j]) {
            std::cerr << "issue at zone " << i << " on rank " << rank << " with elem " << j << std::endl;
            fail = true;
          }
        }
      }
      if (fail) {
        std::cerr << "zone restart failed on rank " << my_rank << std::endl;
        exit(3);
      }
      else std::cerr << "zone restart succeeded on rank " << my_rank << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
}


void grid_general::testCheckpointGeneralGrid() {
  writeCheckpointGrid("grid.h5");

  MPI_Barrier(MPI_COMM_WORLD);

  readCheckpointGrid("grid.h5");

  for (int rank = 0; rank < nproc; rank++) {
    std::cerr << "rank " << my_rank << std::endl;
    if (rank == my_rank) {
      bool fail = false;
      if (t_now != t_now_new) {
        std::cerr << "issue at t_now on rank " << rank << std::endl;
        fail = true;
      }
      if (n_zones != n_zones_new) {
        std::cerr << "issue at n_zones on rank " << rank << std::endl;
        fail = true;
      }
      if (n_elems != n_elems_new) {
        std::cerr << "issue at n_elems on rank " << rank << std::endl;
        fail = true;
      }
      for (int i = 0; i < n_elems_new; i++) {
        if (elems_Z[i] != elems_Z_new[i]) {
          std::cerr << "issue at elems_Z, element number " << i << " on rank " << rank << std::endl;
          std::cerr << elems_Z[i] << " " << elems_Z_new[i] << std::endl;
          fail = true;
        }
        if (elems_A[i] != elems_A_new[i]) {
          std::cerr << "issue at elems_A, element number " << i << " on rank " << rank << std::endl;
          std::cerr << elems_A[i] << " " << elems_A_new[i] << std::endl;
          fail = true;
        }
      }
      if (fail) {
        std::cerr << "grid restart failed on rank " << my_rank << std::endl;
        exit(3);
      }
      else std::cerr << "general grid restart succeeded on rank " << my_rank << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
}
