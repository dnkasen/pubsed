#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include "transport.h"
#include "ParameterReader.h"
#include "spectrum_array.h"
#include "grid_general.h"
#include "grid_1D_sphere.h"

int main(int argc, char **argv){
    MPI_Init(NULL, NULL);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    const int verbose = (rank == 0);

    std::string param_file = "param.lua";
    if( argc > 1 ) param_file = std::string( argv[ 1 ] );
    ParameterReader params(param_file,verbose);

    std::vector<std::string> fnames = params.getVector<std::string>("spectrum_calc_particle_files");
    std::vector<int> file_rank_id = params.getVector<int>("spectrum_calc_file_to_rank");
    std::string chk_spectrum_fname = params.getScalar<std::string>("spectrum_calc_chk_file");
    std::string out_spectrum_fname = params.getScalar<std::string>("spectrum_calc_out_file");
    double rank_factor = params.getScalar<double>("spectrum_calc_n_ranks_old") / nprocs;

    if (fnames.size() != file_rank_id.size()) {
      std::cerr << "File list and rank correspondence not same size. Exiting." << std::endl;
      exit(1);
    }
    // See which files each rank will take.
    int n_files = fnames.size();

    std::vector<std::string> my_fnames;

    for (int i = 0; i < n_files; i ++) {
      if (file_rank_id[i] == rank) {
        my_fnames.push_back(fnames[i]);
      }
    }
    // read in and wipe reference spectrum to get the grids right
    spectrum_array spectrum;
    spectrum.readCheckpointSpectrum(chk_spectrum_fname, "optical spectrum");
    spectrum.wipe();
    spectrum.set_name(out_spectrum_fname);
    
    // minimum grid initialization for transport class
    grid_general* grid = new grid_1D_sphere;
    grid->n_zones = 1;
    grid->elems_A.push_back(1);
    grid->elems_Z.push_back(1);

    transport* transport_dummy = new transport;
    transport_dummy->init(&params, grid);
    for (auto i_fname = my_fnames.begin(); i_fname != my_fnames.end(); i_fname++) {
      std::vector<particle> particle_list;
      std::string fname = *i_fname;
      transport_dummy->readCheckpointParticles(particle_list, fname, "particles_escaped", false, true);
      std::cerr << particle_list.size() << std::endl;
      for (auto i_part = particle_list.begin(); i_part != particle_list.end(); i_part++) {
        spectrum.count(i_part->t, i_part->nu, i_part->e / rank_factor, i_part->D);
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    spectrum.MPI_average();
    if (verbose) {
      spectrum.print(1);
    }

    delete grid;
    delete transport_dummy;

    MPI_Finalize();
}
