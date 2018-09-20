#include <limits>

#include "hydro_homologous.h"
#include "grid_1D_sphere.h"

void hydro_homologous::init(ParameterReader *params, grid_general *g)
{
  this->grid = g;
  this->t_start = params->getScalar<double>("tstep_time_start");
}


double hydro_homologous::get_time_step()
{
  return std::numeric_limits<double>::infinity();
}



void hydro_homologous::step(double dt)
{
  double e = (grid->t_now+dt)/grid->t_now;

  for (int i=0;i<grid->n_zones;i++) {
    grid->z[i].rho = grid->z[i].rho/e/e/e;
    // If t_now < t_start, adiabatically expand T as well
    // assuming radiation pressure dominated
    if (grid->t_now < t_start) grid->z[i].T_gas = grid->z[i].T_gas/e;
  }
  grid->expand(e);
  
  grid->t_now += dt;
}
