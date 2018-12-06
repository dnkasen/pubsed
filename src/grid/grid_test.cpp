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

  writeCheckpointZones("test_zones.h5");

  MPI_Barrier(MPI_COMM_WORLD);

  readCheckpointZones("test_zones.h5");

  for (int rank = 0; rank < nprocs; rank++) {
    if (rank == my_rank) {
      bool fail = false;
      for (int i = 0; i < n_zones; i++) {
        if ((z_new[i].v[0] != z[i].v[0]) || (z_new[i].v[1] != z[i].v[1]) || (z_new[i].v[2] != z[i].v[2])
            || (z_new[i].rho != z[i].rho)
            || (z_new[i].cs != z[i].cs)
            || (z_new[i].p_gas != z[i].p_gas)
            || (z_new[i].T_gas != z[i].T_gas)
            || (z_new[i].mu != z[i].mu)
            || (z_new[i].e_rad != z[i].e_rad)
            || (z_new[i].e_abs != z[i].e_abs)
            || (z_new[i].fx_rad != z[i].fx_rad)
            || (z_new[i].fy_rad != z[i].fy_rad)
            || (z_new[i].fz_rad != z[i].fz_rad)
            || (z_new[i].eps_imc != z[i].eps_imc)
            || (z_new[i].L_thermal != z[i].L_thermal)
            || (z_new[i].L_radio_emit != z[i].L_radio_emit)
            || (z_new[i].L_radio_dep != z[i].L_radio_dep)
            || (z_new[i].G1 != z[i].G1)
            || (z_new[i].G2 != z[i].G2)
            || (z_new[i].G3 != z[i].G3)
            || (z_new[i].P11 != z[i].P11)
            || (z_new[i].P12 != z[i].P12)
            || (z_new[i].P13 != z[i].P13)
            || (z_new[i].P22 != z[i].P22)
            || (z_new[i].P23 != z[i].P23)
            || (z_new[i].P33 != z[i].P33))
          fail = true;
        for (int j = 0; j < n_elems; j++) {
          if (z_new[i].X_gas[j] != z[i].X_gas[j]) fail = true;
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




