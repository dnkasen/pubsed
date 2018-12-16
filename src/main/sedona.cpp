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

#ifdef SEDONA_GIT_VERSION
  if (verbose)
    std::cout << "# git version " << std::string(SEDONA_GIT_VERSION) << std::endl;
#endif

#ifdef COMPILE_DATETIME
  if (verbose)
    std::cout << "# compiled on " << std::string(COMPILE_DATETIME) << std::endl;
#endif

  //---------------------------------------------------------------------
  // BEGIN SETTING UP
  //---------------------------------------------------------------------

  // open up the parameter reader
  std::string param_file = "param.lua";
  if( argc > 1 ) param_file = std::string( argv[ 1 ] );
  ParameterReader params(param_file,verbose);


  //---------------------------------------------------------------------
  // SET UP THE GRID
  //---------------------------------------------------------------------
  grid_general *grid;

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
  double t_start = params.getScalar<double>("tstep_time_start");
  if ((hydro_type == "homologous")&&(t_start > 0))
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

  std::cout << std::scientific;
  std::cout << std::setprecision(2);

  // loop over time/iterations
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

#ifdef MPI_PARALLEL
  // finish up mpi
  MPI_Finalize();
#endif

  return 0;
}
