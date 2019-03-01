#ifndef _SEDONA_CLASS_H
#define _SEDONA_CLASS_H

#include "sedona.h"
#include "physical_constants.h"
#include "ParameterReader.h"
#include "grid_general.h"
#include "grid_1D_sphere.h"
#include "grid_2D_cyln.h"
#include "grid_3D_cart.h"
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
  int do_restart_;
  int do_checkpoint_;
  int use_transport_;
  int use_hydro_;
  std::string checkpoint_name_base_;

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
