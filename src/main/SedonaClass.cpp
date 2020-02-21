#include <time.h>
#include <iostream>
#include <iomanip>      // std::setprecision
#include <vector>
#include <ctime>

#include "sedona.h"
#include "SedonaClass.h"
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

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

namespace pc = physical_constants;
using std::string;
using std::cout;
using std::cerr;
using std::endl;

//--------------------------------------------------------
// ***** Run everything ******
//--------------------------------------------------------
int SedonaClass::run(std::string param_file)
{
#ifdef MPI_PARALLEL
  MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
  MPI_Comm_size( MPI_COMM_WORLD, &n_procs);
#else
  my_rank = 0;
  n_procs = 1;
#endif
  verbose_ = (my_rank == 0);

int n_threads = 1;
#ifdef _OPENMP
  n_threads = omp_get_max_threads();
#endif

  // initial printout
  if (verbose_)
  {
    cout << "##################################" << endl;
    cout << "############  sedona  ############" << endl;
    cout << "##################################" << endl;
    cout << "#" << endl;
    cout << "# MPI tasks = " << n_procs << endl;
    cout << "# threads per task = " << n_threads << endl;
    cout << "#" << endl;
  }

  // open the parameter reader and setup everything
  start_timer();
  params_.initialize(param_file,verbose_);
  setup();
  if (verbose_)
    cout << "# Setup Done: took " << get_timer() << " seconds" << endl;

  last_chk_timestep_ = 0;
  last_chk_walltime_ = 0;
  last_chk_simtime_ = 0;

  // evolve the entire system in time (or iteration)
  evolve_system();

  // printout completion
  if (verbose_)
  {
    double time_wasted = get_timer();
    cout << "#" << endl;
    cout << "# CALCULATION took " << time_wasted << " seconds ";
    cout << " = " << time_wasted/60.0 << " minutes";
    cout << " = " << time_wasted/60.0/60 << " hours" << endl;
  }

  // Clean up
  delete grid_;
  if (use_hydro_) delete hydro_;
  if (use_transport_) delete transport_;
  return 0;
}

//--------------------------------------------------------
// Setup everything
//--------------------------------------------------------
void SedonaClass::setup()
{
  std::string git_version = "";
#ifdef SEDONA_GIT_VERSION
  git_version = std::string(SEDONA_GIT_VERSION);
  if (verbose_)
    std::cout << "# git version " << git_version << std::endl;
#endif

  std::string compile_time = "";
#ifdef COMPILE_DATETIME
  compile_time = std::string(COMPILE_DATETIME);
  if (verbose_)
    std::cout << "# compiled on " << compile_time << std::endl;
#endif
  cout << "#" << endl;

  // Handle restart bookkeeping
  do_restart_ = params_.getScalar<int>("run_do_restart");
  do_checkpoint_ = params_.getScalar<int>("run_do_checkpoint");
  chk_timestep_interval_ = params_.getScalar<int>("run_chk_timestep_interval");
  chk_walltime_interval_ = params_.getScalar<double>("run_chk_walltime_interval");
  chk_simtime_interval_ = params_.getScalar<double>("run_chk_simtime_interval");
  chk_walltime_max_buffer_ = params_.getScalar<double>("run_chk_walltime_max_buffer");
  chk_walltime_max_ = params_.getScalar<double>("run_chk_walltime_max");
  int do_checkpoint_test = params_.getScalar<int>("run_do_checkpoint_test");
  i_chk_ = params_.getScalar<int>("run_chk_number_start");

  std::string restart_file;
  if (do_restart_)
  {
    restart_file = params_.getScalar<string>("run_restart_file");
    std::string compile_time_chk, git_version_chk;
    read_checkpoint_meta(restart_file, compile_time_chk, git_version_chk);
    if (verbose_)
    {
      cout << "# Restarting from " << restart_file << endl;
      if (git_version != git_version_chk)
      {
        cerr << "# WARNING: restarting from file generated from git version ";
        cerr << git_version_chk << endl;
        cerr << "# Different from current git version " << git_version << endl;
      }
      else if (compile_time != compile_time_chk)
      {
        cerr << "# WARNING: restarting from file generated from binary";
        cerr << " compiled at " << compile_time_chk << endl;
        cerr << "# Different from current binary compile time ";
        cerr << compile_time << endl;
      }
      else
      {
        cout << "# Restarting from checkpoint file generated from this binary";
        cout << endl;
      }
    }
  }
  if (verbose_)
    cout << "##################################" << endl << "#" << endl;

  if (do_checkpoint_)
  checkpoint_name_base_ = params_.getScalar<string>("run_checkpoint_name_base");
  // TODO: complain if old and new code versions aren't the same.

  //---------------------------------------------------------------------
  // Setup the grid
  //--------------------------------------------------------------------_
  // read grid type and set up appropripately
  grid_ = NULL;
  string grid_type = params_.getScalar<string>("grid_type");

  // create a grid of the appropriate type
  if      (grid_type == "grid_1D_sphere") grid_ = new grid_1D_sphere;
  else if (grid_type == "grid_2D_cyln"  ) grid_ = new grid_2D_cyln;
  else if (grid_type == "grid_3D_cart"  ) grid_ = new grid_3D_cart;
  else
  {
    if(verbose_) cerr << "# ERROR: the grid type is not implemented" << endl;
    exit(3);
  }
  // initialize the grid (including reading the model file)
  grid_->init(&params_);

  //---------------------------------------------------------------------
  // Set up the hydro module
  //---------------------------------------------------------------------
  hydro_ = NULL;
  string hydro_type = params_.getScalar<string>("hydro_module");

  if (hydro_type == "homologous")
    hydro_ = new hydro_homologous;
  else if (hydro_type == "none")
    hydro_ = NULL;
  else if (hydro_type == "1D_lagrangian")
    hydro_ = new hydro_1D_lagrangian;
  else
  {
    if (verbose_) cerr << "# ERROR: the hydro type is not implemented" << endl;
    exit(3);
  }
  use_hydro_ = (hydro_ != NULL);
  if (use_hydro_) hydro_->init(&params_, grid_);

  // At this point, check if we should evolve system to user start time
  evolve_to_start();

  //---------------------------------------------------------------------
  // Set up the transport module
  //---------------------------------------------------------------------
  string transport_type = params_.getScalar<string>("transport_module");
  use_transport_ = 0;
  // only have 1 transport module now = monte carlo
  if (transport_type != "")
  {
    use_transport_ = 1;
    transport_ = new transport;
    transport_->init(&params_, grid_);
  }

  // you can immediately write out a checkpoint to test correctness
  if ((do_restart_) && (do_checkpoint_test))
  {
    std::string checkpoint_file_init = checkpoint_name_base_ + "_init.h5";
    write_checkpoint(checkpoint_file_init);
  }



}

//-----------------------------------------------------------
// For homologous models, evolve from the model time
// to the user parameter specified time
//-----------------------------------------------------------
void SedonaClass::evolve_to_start()
{
  // Evolve to start time if homologous
  // by adiabatically expanding or compress rho and T
  // But this should only happen if there's not a restart
  double t_start    = params_.getScalar<double>("tstep_time_start");
  string hydro_type = params_.getScalar<string>("hydro_module");
  if ((hydro_type == "homologous")&&(t_start > 0) && (not do_restart_))
  {
    // printout what we are going to do to the model
    if (verbose_)
    {
      cout << "# t_start = " << t_start << ", t_now = " << grid_->t_now;
      if (t_start < grid_->t_now)
      {
        cout << "#" << endl;
        cout << "# Adiabatically compressing rho and T from input model" << endl;
        cout << "# Assumes rad pressure dominated; no radioactive heating" << endl;
      }
      else
      {
        cout << "#" << endl;
        cout << "# Adiabatically expanding rho and T from input model" << endl;
        cout << "# Assumes rad pressure dominated; includes radioactive heating" << endl;
      }
    }
    int force_rproc  = params_.getScalar<int>("force_rprocess_heating");
    hydro_->evolve_to_start(t_start, force_rproc);
    grid_->t_now = t_start;
  }
}


//-----------------------------------------------------------
// Do main time (or iteration loop)
//-----------------------------------------------------------
void SedonaClass::evolve_system()
{
  // check for steady state iterative calculation
  // or a time dependent calculation
  int steady_iterate  = params_.getScalar<int>("transport_steady_iterate");
  int n_steps   = steady_iterate;
  if (steady_iterate) use_hydro_ = 0;

  // read in time stepping parameters
  double t_stop = 0;
  if (!steady_iterate)
  {
    n_steps  = params_.getScalar<int>("tstep_max_steps");
    t_stop   = params_.getScalar<double>("tstep_time_stop");
  }
  double dt_max  = params_.getScalar<double>("tstep_max_dt");
  double dt_min  = params_.getScalar<double>("tstep_min_dt");
  double dt_del  = params_.getScalar<double>("tstep_max_delta");
  // TODO: read in last time step information if restart. Given current time-
  // stepping scheme, not strictly necessary

  // parameters for writing data to file
  int write_levels      = params_.getScalar<int>("output_write_atomic_levels");
  int write_radiation   = params_.getScalar<int>("output_write_radiation");
  int write_mass_fracs  = params_.getScalar<int>("output_write_mass_fractions");
  double write_out_step = params_.getScalar<double>("output_write_plt_file_time");
  double write_out_log  = params_.getScalar<double>("output_write_plt_log_space");
  int    i_write = 0;
  double next_write_out = grid_->t_now;

  // print out initial state
  if (verbose_)
  {
    cout << "# writing initial plot file" << endl;
    grid_->write_plotfile(0,grid_->t_now,write_mass_fracs);
  }

  std::cout << std::scientific;
  std::cout << std::setprecision(2);
  // loop over time/iterations. dt at iteration i doesn't depend on dt at
  // iteration i + 1. Not currently checkpointed, but may need to be later.
  t_ = grid_->t_now;
  it_ = 1;
  while (it_ <= n_steps)
  //for(int it=1; it<=n_steps; it++,t+=dt_)
  {
    start_step_wt_ = get_timer();
    // get this time step
    if (!steady_iterate)
    {
      dt_ = dt_max;
      if (use_hydro_)
      {
        double dt_hydro = hydro_->get_time_step();
        if (dt_hydro < dt_) dt_ = dt_hydro;
      }
      if ((dt_del > 0)&&(t_ > 0)) if (dt_ > t_*dt_del) dt_ = t_*dt_del;
      if (dt_ < dt_min) dt_ =  dt_min;
    }
    else
    {
      dt_ = 0;
      if (it_ == n_steps) transport_->set_last_iteration_flag();
    }

    // printout basic time step information
    if (verbose_)
    {
      cout << "#-------------------------------------------" << endl;
      if (steady_iterate) cout << "# ITERATION: " << it_ << ";  t = " << t_ << "\t";
      else cout << "# TSTEP #" << it_ << " ; t = " << t_ << " sec (" << t_/3600/24.0 << " days); dt = " << dt_;
      cout << endl;
      if (use_transport_) cout << "# particles on grid = " << transport_->n_particles() << endl;
    }

    // do hydro step
    if (use_hydro_) hydro_->step(dt_);

    if (steady_iterate)	transport_->wipe_spectra();

    // do transport step
    if (use_transport_)
    {
      transport_->write_levels = 0;
      if(((t_>=next_write_out)||(steady_iterate))&&write_levels)
      transport_->write_levels = 1;

      transport_->step(dt_);
      // print out spectrum if an iterative calc
      if (steady_iterate) transport_->output_spectrum(it_);
    }

    // writeout output files when appropriate
    if ((t_ >= next_write_out)||(steady_iterate))
    {
      double t_write = t_ + dt_;
      if (steady_iterate) t_write = t_;

      // write out what we want
      if (verbose_)
      {
        cout << "# writing plot file " << i_write + 1;
        cout << " at time " << t_write << endl;
        grid_->write_plotfile(i_write+1,t_write,write_mass_fracs);
        if ((use_transport_)&&(write_radiation)) {
          transport_->write_radiation_file(i_write+1);
          if(write_levels) transport_->write_levels_to_plotfile(i_write+1);
        }
      }

      //write spectrum
      if (use_transport_)
        transport_->output_spectrum(i_write+1);


      // determine next write out
      if ((write_out_log > 0)&&(i_write > 0))
      next_write_out = next_write_out*(1.0 + write_out_log);
      else
      next_write_out = next_write_out + write_out_step;
      i_write++;
    }

    int chk_now;
#ifdef MPI_PARALLEL
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    if (verbose_) {
      chk_now = do_checkpoint_now();
    }
#ifdef MPI_PARALLEL
    MPI_Bcast(&chk_now, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
    if ((do_checkpoint_) && (chk_now))
    {
      write_checkpoint(i_chk_);
      i_chk_++;
    }

    // check for end
    if ((!steady_iterate)&&(t_ > t_stop)) break;
    it_++;
    t_ += dt_;
  }

  // print out final spectrum
  if ((use_transport_)&&(!steady_iterate))  transport_->output_spectrum(-1);
  if (do_checkpoint_)
  write_checkpoint(checkpoint_name_base_ + "_final.h5");
}

int SedonaClass::do_checkpoint_now(int chk_force) {
  int chk_timestep = do_checkpoint_iteration();
  int chk_walltime = do_checkpoint_walltime();
  int chk_end = do_checkpoint_before_end();
  int chk_simtime = do_checkpoint_simulation_time();
  int chk_total = chk_timestep + chk_walltime + chk_end + chk_simtime + chk_force;
  return chk_total;
}

int SedonaClass::do_checkpoint_iteration() {
  if (chk_timestep_interval_ <= 0) {
    return 0;
  }
  else
  {
    if (it_ >= last_chk_timestep_ + chk_timestep_interval_) {
      return 1;
    }
  }
  return 0;
}

int SedonaClass::do_checkpoint_walltime() {
  if (chk_walltime_interval_ <= 0) {
    return 0;
  }
  else
  {
    double curr_time = get_timer();
    if (curr_time >= last_chk_walltime_ + chk_walltime_interval_) {
      return 1;
    }
  }
  return 0;
}

int SedonaClass::do_checkpoint_before_end() {
  if ((chk_walltime_max_ <= 0) || (chk_walltime_max_buffer_ == 0)) {
    return 0;
  }
  else
  {
    double curr_time = get_timer();
    double last_step = curr_time - start_step_wt_;
    if ((curr_time + chk_walltime_max_buffer_ * last_step) > chk_walltime_max_)
    {
      return 1;
    }
  }
  return 0;
}

int SedonaClass::do_checkpoint_simulation_time() {
  if (chk_simtime_interval_ <= 0) {
    return 0;
  }
  else
  {
    if (t_ >= last_chk_simtime_ + chk_simtime_interval_) {
      return 1;
    }
  }
  return 0;
}

//-----------------------------------------------------------
// Write a checkpoint file with numbered name
//-----------------------------------------------------------
void SedonaClass::write_checkpoint(int i_chk_)
{
  string i_chk_str = std::to_string(i_chk_);
  i_chk_str.insert(i_chk_str.begin(), 5 - i_chk_str.length(), '0');
  string checkpoint_file_full = checkpoint_name_base_ + "_" + i_chk_str + ".h5";
  write_checkpoint(checkpoint_file_full);
}

//-----------------------------------------------------------
// Write a checkpoint file with passed name
//-----------------------------------------------------------
void SedonaClass::write_checkpoint(std::string checkpoint_file_full)
{
  if (verbose_) {
    cout << "# writing checkpoint file " << checkpoint_file_full;
    cout << " at time " << t_ << endl;
    createFile(checkpoint_file_full);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  write_checkpoint_meta(checkpoint_file_full);
  transport_->writeCheckpointParticles(checkpoint_file_full);
  grid_->writeCheckpointZones(checkpoint_file_full);
  transport_->writeCheckpointSpectra(checkpoint_file_full);
  transport_->writeCheckpointRNG(checkpoint_file_full);
  grid_->writeCheckpointGrid(checkpoint_file_full);
  last_chk_timestep_ = it_;
  last_chk_walltime_ = get_timer();
  last_chk_simtime_ = t_;
}

//-----------------------------------------------------------
// Write meta data re git version, compile time
//-----------------------------------------------------------
void SedonaClass::write_checkpoint_meta
(std::string checkpoint_file)
{
  if (verbose_)
  {
    createGroup(checkpoint_file, "meta");
#ifdef COMPILE_DATETIME
    std::string datetime_string = std::string(COMPILE_DATETIME);
#else
    std::string datetime_string = "";
#endif
    writeString(checkpoint_file, "meta", "compile_time", datetime_string);

#ifdef SEDONA_GIT_VERSION
    std::string version_string = std::string(SEDONA_GIT_VERSION);
#else
    std::string version_string = "";
#endif
    writeString(checkpoint_file, "meta", "git_version", version_string);
  }
}


//-----------------------------------------------------------
// check if meta data from checkpoint consistent
//-----------------------------------------------------------
void SedonaClass::read_checkpoint_meta
(std::string checkpoint_file, std::string& datetime_string_new,
std::string& version_string_new)
{

  for (int rank = 0; rank < n_procs; rank++)
  {
    if (rank == my_rank)
    {
      readString(checkpoint_file, "meta", "compile_time", datetime_string_new);
      readString(checkpoint_file, "meta", "git_version", version_string_new);
    }
    #ifdef MPI_PARALLEL
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
  }
}



//-----------------------------------------------------------
// Functions to start and stop the timer
//-----------------------------------------------------------
void SedonaClass::start_timer()
{
  // start timer
#ifdef MPI_PARALLEL
  timer_start_ = MPI_Wtime();
#else
  timer_start_ = clock();
#endif
}

double SedonaClass::get_timer()
{
  // stop timer
  double time_wasted;
#ifdef MPI_PARALLEL
  double timer_stop = MPI_Wtime();
  time_wasted = timer_stop - timer_start_;
#else
  clock_t timer_stop = clock();
  time_wasted = double(timer_stop - timer_start_) / CLOCKS_PER_SEC;
#endif
  return time_wasted;
}

double SedonaClass::get_and_reset_timer()
{
  double time_wasted = get_timer();
  start_timer();
  return time_wasted;
}
