#include <mpi.h>
#include <time.h>
#include <iostream>
#include <vector>

#include "physical_constants.h"
#include "ParameterReader.h"
#include "grid_general.h"
#include "grid_1D_sphere.h"
#include "hydro_general.h"
#include "hydro_homologous.h"
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
    cout << "################################\n";
    cout << "##########  SEDONA  ############\n";
    cout << "################################\n";
    cout << "# MPI cores used = " << n_procs << endl;
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
  if     (grid_type == "grid_1D_sphere") grid = new grid_1D_sphere;
  //else if(grid_type == "grid_3D_cart"  ) grid = new grid_3D_cart;
  else  {
    if(verbose) cout << "# ERROR: the grid type is not implemented\n";
    exit(3);   }
  
  // initialize the grid (including reading the model file)
  grid->init(&params);

  
  //---------------------------------------------------------------------
  // SET UP the transport module
  //---------------------------------------------------------------------
  transport mcarlo;
  mcarlo.init(&params, grid);
  
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
  
  // check for steady state iterative calculation
  // or a time dependent calculation
  int n_steps   = steady_iterate;
  double t_stop = 0; 
  if (!steady_iterate)  {
    n_steps  = params.getScalar<int>("tstep_max_steps");
    t_stop   = params.getScalar<double>("tstep_time_stop"); }

  // parameters for writing data to file
  double write_out   = params.getScalar<double>("grid_write_out");
  int    iw = 0;

  // loop over time/iterations
  double dt, t = grid->t_now;
  for(int it=1; it<=n_steps; it++,t+=dt)
  {
    // get this time step
    if (!steady_iterate)
    {
      double dt_max  = params.getScalar<double>("tstep_max_dt");
      double dt_min  = params.getScalar<double>("tstep_min_dt");
      double dt_del  = params.getScalar<double>("tstep_max_delta");
      dt = dt_max;
      if ((dt_del > 0)&&(t > 0)) if (dt > t*dt_del) dt = t*dt_del;
      if (dt < dt_min) dt =  dt_min;
    }
    else dt = 0;

    // printout time step
    if (verbose) 
    {
      if (steady_iterate) cout << "# ITERATION: " << it << "\t";
      else cout << "# TIME STEP: " << it << "\t" << t << "\t" << dt << "\t";
      cout << mcarlo.n_particles() << "\n";
    }

    // integrate mass 
     double mass =0;
     for (int i=0;i<grid->n_zones;i++)
       mass += grid->z[i].rho*grid->zone_volume(i);
     //std::cout << "mass =" << mass/pc::m_sun << "\n";

    // do hydro step
    if (use_hydro) hydro->step(dt);

    // do transport step
    mcarlo.step(dt);

    // print out spectrum if an iterative calc
    if (steady_iterate) mcarlo.output_spectrum(it);

    // writeout zone state when appropriate 
    if ((verbose)&&((t >= write_out*iw)||(steady_iterate)))
    {
      double t_write = t + dt;
      if (steady_iterate) t_write = t;
      printf("# writing zone file %d at time %e\n",iw+1, t_write);
      grid->write_out(iw+1,t_write);
      iw++;
    }

    // check for end
    if ((!steady_iterate)&&(t > t_stop)) break;
  }

  // print out final spectrum
  if (!steady_iterate)  mcarlo.output_spectrum();

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
