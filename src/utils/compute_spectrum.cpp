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
#include "h5utils.h"
#include "physical_constants.h"

namespace pc = physical_constants;

int between(double val, std::vector<double> bounds) {
  if (bounds.size() == 2)
    return ((val > bounds[0]) and (val <= bounds[1]));
  else if (bounds.size() == 0)
    return 1;
  else {
    std::cerr << "Incorrectly sized bounds vector with " << bounds.size() << " elements. " <<
      "Should be 0 or 2." << std::endl;
    exit(4);
  }
}

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
    int save_particles = params.getScalar<int>("spectrum_calc_save_particles");
    std::string save_particles_fname = params.getScalar<std::string>("spectrum_calc_save_particles_file");

    std::vector<double> time_filter = params.getVector<double>("spectrum_calc_time_filter");
    std::vector<double> time_phys_filter = params.getVector<double>("spectrum_calc_time_phys_filter");
    std::vector<double> energy_filter = params.getVector<double>("spectrum_calc_energy_filter");
    std::vector<double> nu_filter = params.getVector<double>("spectrum_calc_nu_filter");
    std::vector<double> mu_filter = params.getVector<double>("spectrum_calc_mu_filter");
    std::vector<double> vel_filter = params.getVector<double>("spectrum_calc_velocity_filter");


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
    transport_dummy->setup_MPI();
    std::vector<particle> saved_particles;
    for (auto i_fname = my_fnames.begin(); i_fname != my_fnames.end(); i_fname++) {
      std::vector<particle> particle_list;
      std::string fname = *i_fname;
      transport_dummy->readCheckpointParticles(particle_list, fname, "particles_escaped", false, true);
      // Divide particle energy by number of processes to undo multiplication that happens on
      // reading in a particle list from checkpoint
      for (auto i_part = particle_list.begin(); i_part != particle_list.end(); i_part++) {
        i_part->e = i_part->e / nproc;
      }
      for (auto i_part = particle_list.begin(); i_part != particle_list.end(); i_part++) {
        double time_phys = i_part->t + i_part->x_dot_d() / pc::c;
        double x_inter_x_sep[3] = {i_part->x_interact[0] - i_part->x[0],
          i_part->x_interact[1] - i_part->x[1], i_part->x_interact[2] - i_part->x[2]};
        double x_interact_dist = sqrt(x_inter_x_sep[0] * x_inter_x_sep[0] +
              x_inter_x_sep[1] * x_inter_x_sep[1] + x_inter_x_sep[2] * x_inter_x_sep[2]);
        double time_interact = time_phys - x_interact_dist / pc::c;
        int time_filt_flag = between(i_part->t, time_filter);
        int time_phys_filt_flag = between(time_phys, time_phys_filter);
        int energy_filt_flag = between(i_part->e, energy_filter);
        int mu_filt_flag = between(i_part->D[2], mu_filter);
        int nu_filt_flag = between(i_part->nu, nu_filter);
        int vel_filt_flag = between(i_part->r_interact() / time_interact / pc::c, vel_filter);
        int filt_flag = time_filt_flag * time_phys_filt_flag * energy_filt_flag *
            mu_filt_flag * nu_filt_flag * vel_filt_flag;

        if (filt_flag) {
          spectrum.count(i_part->t, i_part->nu, i_part->e, i_part->D);
          if (save_particles)
            saved_particles.push_back(*i_part);
        } 
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    spectrum.MPI_average();
    if (verbose) {
      spectrum.print(1);
    }

    if (save_particles) {
      for (auto i_part = saved_particles.begin(); i_part != saved_particles.end(); i_part++) {
        i_part->e = i_part->e;
      }
      if (verbose)
        createFile(save_particles_fname);
      transport_dummy->writeCheckpointParticles(saved_particles, save_particles_fname, "particles_filtered");
    }

    delete grid;
    delete transport_dummy;

    MPI_Finalize();
}
