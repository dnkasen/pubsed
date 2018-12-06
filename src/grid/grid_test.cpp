#include <cstdlib>
#include <math.h>

#include "grid_general.h"
#include "physical_constants.h"
namespace pc = physical_constants;

void grid_general::testCheckpointZones() {
  int my_rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
  int nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  std::cerr << "starting to test zones" << std::endl;

  writeCheckpointZones("test_zones.h5");

  std::cerr << "wrote checkpoint zones" << std::endl;

  MPI_Barrier(MPI_COMM_WORLD);

  readCheckpointZones("test_zones.h5");

  std::cerr << "read checkpoint zones" << std::endl;

  for (int rank = 0; rank < nprocs; rank++) {
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




