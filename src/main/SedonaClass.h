#ifndef _SEDONA_CLASS_H
#define _SEDONA_CLASS_H

#include "sedona.h"
#include "physical_constants.h"
#include "ParameterReader.h"
#include "grid_general.h"
#include "grid_1D_sphere.h"
#include "grid_2D_cyln.h"
#include "grid_3D_cart.h"
#include "grid_3D_sphere.h"
#include "hydro_general.h"
#include "hydro_homologous.h"
#include "hydro_1D_lagrangian.h"
#include "transport.h"

class SedonaClass
{

private:

  // The core classes used in the code
  grid_general   *grid_;
  hydro_general  *hydro_;
  transport      *transport_;
  ParameterReader params_;

  // some core control parameters
  int verbose_;
  int my_rank; // MPI info 
  int n_procs;
  int do_restart_;
  int do_checkpoint_;
  int use_transport_;
  int use_hydro_;
  std::string checkpoint_name_base_;

  int last_chk_timestep_;
  double last_chk_walltime_;
  double last_chk_simtime_;
  double start_step_wt_;

  int chk_timestep_interval_;
  double chk_walltime_interval_;
  double chk_simtime_interval_;
  double chk_walltime_max_buffer_;
  double chk_walltime_max_;

  int it_;
  double t_;
  double dt_;
  int i_chk_;

#ifdef MPI_PARALLEL
  double timer_start_;
#else
  clock_t timer_start_;
#endif

  // calling run() will execute these
  void setup();
  void evolve_to_start();
  void evolve_system();
  void finish();

  int do_checkpoint_now(int chk_force = 0);
  int do_checkpoint_iteration();
  int do_checkpoint_walltime();
  int do_checkpoint_before_end();
  int do_checkpoint_simulation_time();
  int do_checkpoint_triggered();

  // write a checkpoint file with numbered name
  void write_checkpoint(int i_chk);
  // write a checkpoint file with the passed name
  void write_checkpoint(std::string checkpoint_file_full);
  // read and write meta data (e.g., git version of code)
  void write_checkpoint_meta(std::string);
  void read_checkpoint_meta(std::string,std::string&,std::string&);
  // control the timer
  void start_timer();
  double get_timer();
  double get_and_reset_timer();

public:

  // this does everything, based on the passed param file
  int run(std::string param_file);

};

#endif
