#include <limits>

#include "hydro_homologous.h"
#include "grid_1D_sphere.h"
#include "radioactive.h"
#include <stdlib.h>
#include "physical_constants.h"

namespace pc = physical_constants;

using std::cerr;


void hydro_homologous::init(ParameterReader *params, grid_general *g)
{
  this->grid = g;
}


double hydro_homologous::get_time_step()
{
  return std::numeric_limits<double>::infinity();
}



void hydro_homologous::step(double dt)
{
  double e = (grid->t_now+dt)/grid->t_now;

  for (int i=0;i<grid->n_zones;i++) 
    grid->z[i].rho = grid->z[i].rho/e/e/e;
  grid->expand(e);
  
  grid->t_now += dt;
}



// Adiabatically expands or compresses temperature from t_now to t_start
// Accounts for radioactive energy deposition if t_start > t_now
void hydro_homologous::evolve_to_start(double t_start, int force_rproc)
{
  radioactive radio;
  double gfrac;
  double e;
  double dt;
  double eps_nuc;
  double u_old;
  double time_tol = 1.0e-4;
  int incr_max = 9999;

  // Compress and return if t_start < current time
  if (t_start < grid->t_now)
  {
    e = t_start/grid->t_now;

    // Compress rho and T_gas
    // Set T_rad to T_gas
    for (int i=0; i<grid->n_zones; i++)
    {
      grid->z[i].rho = grid->z[i].rho/e/e/e;
      grid->z[i].T_gas = grid->z[i].T_gas/e;
      grid->z[i].e_rad = pc::a*pow(grid->z[i].T_gas,4);
    }

    // Update grid and time
    grid->t_now = t_start;
    grid->expand(e);

    return;
  }

  // If we get to here, t_start > current time, so we expand
  // Just forward Euler for now
  for (int incr=0; incr <= incr_max; incr++)
  {
    if (incr == incr_max)
    {
      cerr << "Reached " << incr_max << " iterations evolving to t_start\n";
      exit(1);
    }

    // Set dt to change lnT by at most 0.1 over all zones
    dt = 1e99;
    for (int i=0;i<grid->n_zones;i++)
    {
      eps_nuc = 
        radio.decay(grid->elems_Z,grid->elems_A,grid->z[i].X_gas,grid->t_now,&gfrac,force_rproc);
      u_old = pc::a*pow(grid->z[i].T_gas,4)/grid->z[i].rho;

      dt = std::min(dt, fabs( 0.1/(eps_nuc/(4.0*u_old) - 1.0/grid->t_now)) );
    }

    // If t_now is within dt of t_start, adjust
    if (grid->t_now + dt > t_start) dt = t_start - grid->t_now;

    e = (grid->t_now + dt)/grid->t_now;

    // Find new T and rho
    for (int i=0; i<grid->n_zones; i++)
    {
      eps_nuc =
        radio.decay(grid->elems_Z,grid->elems_A,grid->z[i].X_gas,grid->t_now,&gfrac,force_rproc);
      u_old = pc::a*pow(grid->z[i].T_gas,4)/grid->z[i].rho;
      grid->z[i].T_gas +=
        grid->z[i].T_gas*dt*(eps_nuc/(4.0*u_old) - 1.0/grid->t_now);
      grid->z[i].rho = grid->z[i].rho/e/e/e;
    }

    // Update grid and time
    grid->t_now += dt;
    grid->expand(e);

    // If we made it to t_start, exit the loop
    if ( fabs(grid->t_now-t_start)/t_start < time_tol ) 
    {
      grid->t_now = t_start; // Make it exact
      break;
    }

  }

  // Set T_rad to T_gas
  for (int i=0; i<grid->n_zones; i++)
  {
    grid->z[i].e_rad = pc::a*pow(grid->z[i].T_gas,4);
  }

}
