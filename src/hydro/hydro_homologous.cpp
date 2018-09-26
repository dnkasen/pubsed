#include <limits>

#include "hydro_homologous.h"
#include "grid_1D_sphere.h"
#include "radioactive.h"
#include <stdlib.h>
#include "physical_constants.h"

namespace pc = physical_constants;


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



// Adiabatically expands temperature from input time to start time
// Accounts for radioactive energy deposition
void hydro_homologous::evolve_to_start(double t_start, int force_rproc)
{
  radioactive radio;
  double gfrac;
  double e;
  double dt = 0.1*grid->t_now; // Timestep should be << relevant timescales
  double dt_max = 100.0;  
  double eps_nuc;
  double u_old;
  double time_tol = 1.0e-4;

  // Just forward Euler for now
  for ( ; ; )
  {
    e = (grid->t_now + dt)/grid->t_now;

    for (int i=0;i<grid->n_zones;i++)
    {
      eps_nuc = 
        radio.decay(grid->elems_Z,grid->elems_A,grid->z[i].X_gas,grid->t_now,&gfrac,force_rproc);
      u_old = pc::a*pow(grid->z[i].T_gas,4)/grid->z[i].rho;
      grid->z[i].T_gas +=
        grid->z[i].T_gas*dt*(eps_nuc/(4.0*u_old) - 1.0/grid->t_now);

      grid->z[i].rho = grid->z[i].rho/e/e/e;
    }

    grid->expand(e);

    // Update time
    grid->t_now += dt;

    // If we made it to t_start, exit
    if ( fabs(grid->t_now-t_start)/t_start < time_tol ) 
    {
      grid->t_now = t_start; // Make it exact
      break;
    }

    // Adjust dt; don't let it go above some set value (currently = 100s)
    dt = std::min(0.1*grid->t_now, dt_max);

    // If t_now is within dt of t_start, adjust dt
    if (grid->t_now + dt > t_start) dt = t_start - grid->t_now;
  }

  // Set T_rad to T_gas
  for (int i=0; i<grid->n_zones; i++)
  {
    grid->z[i].e_rad = pc::a*pow(grid->z[i].T_gas,4);
  }

}
