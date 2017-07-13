#include <mpi.h>
#include <time.h>
#include <iostream>
#include <vector>

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

namespace pc = physical_constants;
using std::string;
using std::cout;
using std::endl;


//--------------------------------------------------------
// The main code
//--------------------------------------------------------
int main(int argc, char **argv)
{
 // initialize MPI parallelism
  int my_rank,n_procs;
  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
  MPI_Comm_size( MPI_COMM_WORLD, &n_procs);

  // verbocity
  const int verbose = (my_rank == 0);
  if (verbose) 
  {
    cout << "##################################\n";
    cout << "############  sedona  ############\n";
    cout << "##################################\n";
    cout << "#\n# MPI tasks = " << n_procs << endl << "#" << endl;

  }
  
  // start timer 
  double proc_time_start = MPI_Wtime();
 
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
    if(verbose) cout << "# ERROR: the grid type is not implemented\n";
    exit(3);   }
  
  // initialize the grid (including reading the model file)
  grid->init(&params);


  //---------------------------------------------------------------------
  // SET UP the transport module
  //---------------------------------------------------------------------
  transport mcarlo;
  string transport_type = params.getScalar<string>("transport_module");
  int use_transport = 0;
  if (transport_type != "") use_transport = 1;
  if (use_transport) mcarlo.init(&params, grid);
  

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
    if (verbose) cout << "# ERROR: the hydro type is not implemented\n";
    exit(3); 
  }
  int use_hydro = (hydro != NULL);
  if (use_hydro) hydro->init(&params, grid);

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
  int write_levels = params.getScalar<int>("output_write_levels");
  int write_grid   = params.getScalar<int>("output_write_grid");
  double write_out_step = params.getScalar<double>("output_write_times");
  double write_out_log  = params.getScalar<double>("output_write_log_times");
  int    i_write = 0;
  // number iteration output from 1
  if (steady_iterate) i_write = 1;
  double next_write_out = grid->t_now;


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
    else dt = 0;

    // printout time step
    if (verbose) 
    {
      if (steady_iterate) cout << "# ITERATION: " << it << "\t";
      else cout << "# TSTEP #" << it << " ; t = " << t << " sec (" << t/3600/24.0 << " days); dt = " << dt;
      cout << "; particles on grid = " << mcarlo.n_particles() << "\n";
    }

    // do hydro step
    if (use_hydro) hydro->step(dt);

    // do transport step
    if (use_transport) 
    {
      mcarlo.step(dt);
      // print out spectrum if an iterative calc
      if (steady_iterate) mcarlo.output_spectrum(it);
    }

    // writeout zone state when appropriate 
    if ((verbose)&&((t >= next_write_out)||(steady_iterate)))
    {
      double t_write = t + dt;
      if (steady_iterate) t_write = t;
      printf("# writing zone file %d at time %e\n",i_write+1, t_write);
      grid->write_out(i_write,t_write);

      if ((write_out_log > 0)&&(i_write > 0))
        next_write_out = next_write_out*(1.0 + write_out_log);
      else
        next_write_out = next_write_out + write_out_step;

      if (use_transport)
      {
        if (write_grid)   mcarlo.write_opacities(i_write);
        if (write_levels) mcarlo.write_levels(i_write);
      }
      i_write++;
    }

    // check for end
    if ((!steady_iterate)&&(t > t_stop)) break;
  }

  // print out final spectrum
  if ((use_transport)&&(!steady_iterate))  mcarlo.output_spectrum();

  //---------------------------------------------------------------------
  // CALCULATION DONE; WRITE OUT AND FINISH
  //---------------------------------------------------------------------
  
  // calculate the elapsed time 
  double proc_time_end = MPI_Wtime();
  double time_wasted = proc_time_end - proc_time_start;
  
  if (verbose)
    printf("#\n# CALCULATION took %.3e seconds or %.3f mins or %.3f hours\n",
	   time_wasted,time_wasted/60.0,time_wasted/60.0/60.0);
  
  // finish up mpi
  MPI_Finalize();
  
  return 0;
}
