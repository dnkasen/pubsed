#include <limits>

#include "hydro_homologous.h"
#include "grid_1D_sphere.h"

void hydro_homologous::init(ParameterReader *params, grid_general *g)
{
  this->grid = g;
  this->t_now = g->t_now;
}


double hydro_homologous::get_time_step()
{
  return std::numeric_limits<double>::infinity();
}



void hydro_homologous::step(double dt)
{
  double e = (t_now+dt)/t_now;

  for (int i=0;i<grid->n_zones;i++) 
    grid->z[i].rho = grid->z[i].rho/e/e/e;
  grid->expand(e);
  
  t_now += dt;
}
