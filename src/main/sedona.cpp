#include <time.h>
#include <iostream>
#include <iomanip>      // std::setprecision
#include <vector>
#include <ctime>

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

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

namespace pc = physical_constants;
using std::string;
using std::cout;
using std::cerr;
using std::endl;



namespace
{
  void writeCheckpointMeta(int verbose, std::string checkpoint_file) {
    if (verbose) {
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

  void readCheckpointMeta(std::string checkpoint_file, std::string& datetime_string_new,
      std::string& version_string_new) {
    int my_rank, n_procs;
#ifdef MPI_PARALLEL
    MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
    MPI_Comm_size( MPI_COMM_WORLD, &n_procs);
#else
    my_rank = 0;
    n_procs = 1;
#endif

    for (int rank = 0; rank < n_procs; rank++) {
      if (rank == my_rank) {
        readString(checkpoint_file, "meta", "compile_time", datetime_string_new);
        readString(checkpoint_file, "meta", "git_version", version_string_new);
      }
#ifdef MPI_PARALLEL
      MPI_Barrier(MPI_COMM_WORLD);
#endif
    }
  }

  void writeCheckpoint(int verbose, std::string checkpoint_file, 
      string chk_label, transport &mcarlo, grid_general* grid) {
    string checkpoint_file_full = checkpoint_file + "_" + chk_label + ".h5";
    if (verbose)
    {
      createFile(checkpoint_file_full);
    }
    writeCheckpointMeta(verbose, checkpoint_file_full);
    mcarlo.writeCheckpointParticles(checkpoint_file_full);
    grid->writeCheckpointZones(checkpoint_file_full);
    mcarlo.writeCheckpointSpectra(checkpoint_file_full);
    grid->writeCheckpointGrid(checkpoint_file_full);
  }

  void writeCheckpoint(int verbose, std::string checkpoint_file, int i_chk, transport &mcarlo, grid_general* grid)
  {
    string i_chk_str = std::to_string(i_chk);
    i_chk_str.insert(i_chk_str.begin(), 5 - i_chk_str.length(), '0');
    writeCheckpoint(verbose, checkpoint_file, i_chk_str, mcarlo, grid);
  }
}


//--------------------------------------------------------
// The main code
//--------------------------------------------------------
int main(int argc, char **argv)
{
 // initialize MPI parallelism
  int my_rank,n_procs;

#ifdef MPI_PARALLEL
  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
  MPI_Comm_size( MPI_COMM_WORLD, &n_procs);
#else
  my_rank = 0;
  n_procs = 1;
#endif

  // verbocity
  const int verbose = (my_rank == 0);
  if (verbose)
  {
    cout << "##################################" << endl;
    cout << "############  sedona  ############" << endl;
    cout << "##################################" << endl;
    cout << "# " << endl;
    cout << "# MPI tasks = " << n_procs << endl << "#" << endl;
  }

// start timer
#ifdef MPI_PARALLEL
  double proc_time_start = MPI_Wtime();
#else
  clock_t time_start = clock();
#endif

  std::string git_version = "";
  std::string compile_time = "";
#ifdef SEDONA_GIT_VERSION
  git_version = std::string(SEDONA_GIT_VERSION);
  if (verbose)
    std::cout << "# git version " << git_version << std::endl;
#endif

#ifdef COMPILE_DATETIME
  compile_time = std::string(COMPILE_DATETIME);
  if (verbose)
    std::cout << "# compiled on " << compile_time << std::endl;
#endif

  //---------------------------------------------------------------------
  // BEGIN SETTING UP
  //---------------------------------------------------------------------

  // open up the parameter reader
  std::string param_file = "param.lua";
  if( argc > 1 ) param_file = std::string( argv[ 1 ] );
  ParameterReader params(param_file,verbose);

  // Handle restart bookkeeping
  int do_restart = params.getScalar<int>("run_do_restart");
  int do_checkpoint = params.getScalar<int>("run_do_checkpoint");
  int do_checkpoint_test = params.getScalar<int>("run_do_checkpoint_test");
  std::string restart_file;
  std::string checkpoint_name_base;
  if (do_restart) {
    restart_file = params.getScalar<string>("run_restart_file");
    cout << "# Restarting from " << restart_file << endl;
    std::string compile_time_chk;
    std::string git_version_chk;
    readCheckpointMeta(restart_file, compile_time_chk, git_version_chk);
    if (verbose) {
      if (git_version != git_version_chk) {
        cerr << "# WARNING: restarting from file generated from git version " << git_version_chk << endl;
        cerr << "# Different from current git version " << git_version << endl;
      }
      else if (compile_time != compile_time_chk) {
        cerr << "# WARNING: restarting from file generated from binary compiled at " << compile_time_chk << endl;
        cerr << "# Different from current binary compile time " << compile_time << endl;
      }
      else {
        cout << "# Restarting from checkpoint file generated from this binary" << endl;
      }
    }
  }
  if (do_checkpoint)
    checkpoint_name_base = params.getScalar<string>("run_checkpoint_name_base");
// TODO: complain if old and new code versions aren't the same.

  //---------------------------------------------------------------------
  // SET UP THE GRID
  //---------------------------------------------------------------------
  grid_general *grid = NULL;

  // read the grid type
  string grid_type = params.getScalar<string>("grid_type");

  // create a grid of the appropriate type
  if      (grid_type == "grid_1D_sphere") grid = new grid_1D_sphere;
  else if (grid_type == "grid_2D_cyln"  ) grid = new grid_2D_cyln;
  else if (grid_type == "grid_3D_cart"  ) grid = new grid_3D_cart;
  else  {
    if(verbose) cerr << "# ERROR: the grid type is not implemented" << endl;
    exit(3);   }

  // initialize the grid (including reading the model file)
  grid->init(&params);



  //---------------------------------------------------------------------
  // SET UP the hydro module
  //---------------------------------------------------------------------
  hydro_general *hydro = NULL;
  string hydro_type = params.getScalar<string>("hydro_module");

  // create a hydro module of the appropriate type
  if (hydro_type == "homologous")
    hydro = new hydro_homologous;
  else if (hydro_type == "none")
    hydro = NULL;
  else if (hydro_type == "1D_lagrangian")
    hydro = new hydro_1D_lagrangian;
  else {
    if (verbose) cerr << "# ERROR: the hydro type is not implemented" << endl;
    exit(3);
  }
  int use_hydro = (hydro != NULL);
  if (use_hydro) hydro->init(&params, grid);

  // Evolve to start time if homologous
  // by adiabatically expanding or compress rho and T
  // But this should only happen if there's not a restart
  double t_start = params.getScalar<double>("tstep_time_start");
  if ((hydro_type == "homologous")&&(t_start > 0) && (not do_restart))
  {
    if (verbose)
    {
      cout << "# t_start = " << t_start << ", t_now = " << grid->t_now;
      if (t_start < grid->t_now)
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
    int force_rproc  = params.getScalar<int>("force_rprocess_heating");
    hydro->evolve_to_start(t_start, force_rproc);
    grid->t_now = t_start;
  }

  //---------------------------------------------------------------------
  // SET UP the transport module
  //---------------------------------------------------------------------
  transport mcarlo;
  string transport_type = params.getScalar<string>("transport_module");
  int use_transport = 0;
  if (transport_type != "") use_transport = 1;
  if (use_transport) mcarlo.init(&params, grid);

  //---------------------------------------------------------------------
  // DO TIME/ITERATION LOOP
  //---------------------------------------------------------------------

  // read in time stepping parameters
  int steady_iterate  = params.getScalar<int>("transport_steady_iterate");
  if (steady_iterate) use_hydro = 0;
  double dt_max  = params.getScalar<double>("tstep_max_dt");
  double dt_min  = params.getScalar<double>("tstep_min_dt");
  double dt_del  = params.getScalar<double>("tstep_max_delta");
  // TODO: read in last time step information if restart. Given current time-
  // stepping scheme, not strictly necessary Might need a sedona class to
  // do that in a non-annoying way, but that's for another time/branch


  // check for steady state iterative calculation
  // or a time dependent calculation
  int n_steps   = steady_iterate;
  double t_stop = 0;
  if (!steady_iterate)
  {
    n_steps  = params.getScalar<int>("tstep_max_steps");
    t_stop   = params.getScalar<double>("tstep_time_stop");
  }

  // parameters for writing data to file
  int write_levels     = params.getScalar<int>("output_write_atomic_levels");
  int write_radiation  = params.getScalar<int>("output_write_radiation");
  int write_mass_fractions  = params.getScalar<int>("output_write_mass_fractions");

  double write_out_step = params.getScalar<double>("output_write_plt_file_time");
  double write_out_log  = params.getScalar<double>("output_write_plt_log_space");
  int    i_write = 0;
  double next_write_out = grid->t_now;

  // print out initial state
  if (verbose)
  {
    cout << "# writing initial plot file" << endl;
    grid->write_plotfile(0,grid->t_now,write_mass_fractions);
  }

  if ((do_restart) && (do_checkpoint_test))
  {
    // Only one rank needs to create the file.
    // Using the verbose variable as a proxy for this when MPI is used
    writeCheckpoint(verbose, checkpoint_name_base, "init", mcarlo, grid);
  }

  std::cout << std::scientific;
  std::cout << std::setprecision(2);

  // loop over time/iterations. dt at iteration i doesn't depend on dt at
  // iteration i + 1. Not currently checkpointed, but may need to be later.
  int i_chk = 0;
  double dt, t = grid->t_now;
  for(int it=1; it<=n_steps; it++,t+=dt)
  {
    // get this time step
    if (!steady_iterate)
    {
      dt = dt_max;
      if (use_hydro)
      {
        double dt_hydro = hydro->get_time_step();
        if (dt_hydro < dt) dt = dt_hydro;
      }
      if ((dt_del > 0)&&(t > 0)) if (dt > t*dt_del) dt = t*dt_del;
      if (dt < dt_min) dt =  dt_min;
    }
    else
    {
      dt = 0;
      if (it == n_steps) mcarlo.set_last_iteration_flag();
    }

    // printout basic time step information
    if (verbose)
    {
      std::cout << "#-------------------------------------------" << std::endl;
      if (steady_iterate) cout << "# ITERATION: " << it << ";  t = " << t << "\t";
      else cout << "# TSTEP #" << it << " ; t = " << t << " sec (" << t/3600/24.0 << " days); dt = " << dt;
      cout << std::endl;
      cout << "# particles on grid = " << mcarlo.n_particles() << std::endl;
    }

    // do hydro step
    if (use_hydro) hydro->step(dt);

    if (steady_iterate)	mcarlo.wipe_spectra();

    // do transport step
    if (use_transport)
    {
      mcarlo.write_levels = 0;
      if(((t>=next_write_out)||(steady_iterate))&&write_levels)
        mcarlo.write_levels = 1;

      mcarlo.step(dt);
      // print out spectrum if an iterative calc
      if (steady_iterate) mcarlo.output_spectrum(it);
    }

    // writeout output files when appropriate
    if ((t >= next_write_out)||(steady_iterate))
    {
      double t_write = t + dt;
      if (steady_iterate) t_write = t;

      // write out what we want
      if (verbose)
      {
        cout << "# writing plot file " << i_write + 1;
        cout << " at time " << t_write << endl;
        grid->write_plotfile(i_write+1,t_write,write_mass_fractions);
        if ((use_transport)&&(write_radiation)){
          mcarlo.write_radiation_file(i_write+1);
          if(write_levels) mcarlo.write_levels_to_plotfile(i_write+1);
        }
      }

      //write spectrum
      if (use_transport)
        mcarlo.output_spectrum(i_write+1);

      if (do_checkpoint) {
        writeCheckpoint(verbose, checkpoint_name_base, i_chk, mcarlo, grid);
        i_chk++;
      }

      // determine next write out
      if ((write_out_log > 0)&&(i_write > 0))
        next_write_out = next_write_out*(1.0 + write_out_log);
      else
        next_write_out = next_write_out + write_out_step;
      i_write++;
    }

    // check for end
    if ((!steady_iterate)&&(t > t_stop)) break;
  }

  // print out final spectrum
  if ((use_transport)&&(!steady_iterate))  mcarlo.output_spectrum(-1);
  if (do_checkpoint)
    writeCheckpoint(verbose, checkpoint_name_base, "final", mcarlo, grid);

  //---------------------------------------------------------------------
  // CALCULATION DONE; WRITE OUT AND FINISH
  //---------------------------------------------------------------------

  // stop timer
double time_wasted;
#ifdef MPI_PARALLEL
  double proc_time_end = MPI_Wtime();
  time_wasted = proc_time_end - proc_time_start;
#else
  clock_t time_stop = clock();
  time_wasted = double(time_stop - time_start) / CLOCKS_PER_SEC;
#endif

  if (verbose)
  {
      cout << "#" << endl;
      cout << "# CALCULATION took " << time_wasted << " seconds ";
      cout << " = " << time_wasted/60.0 << " minutes";
      cout << " = " << time_wasted/60.0/60 << " hours" << endl;
  }

  // Clean up
  delete grid;
  delete hydro;
#ifdef MPI_PARALLEL
  // finish up mpi
  MPI_Finalize();
#endif

  return 0;
}
