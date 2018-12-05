#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <limits>

#include "hydro_1D_movingmesh.h"
#include "physical_constants.h"
namespace pc = physical_constants;


void hydro_1D_movingmesh::init(ParameterReader *params, grid_general *g)
{
  grid = g;

  /*
  gamfac_      = params->getScalar<double>("hydro_gamma_index");
  cfl_         = params->getScalar<double>("hydro_cfl");
  C_q_         = params->getScalar<double>("hydro_viscosity_parameter");
  M_center_    = params->getScalar<double>("hydro_central_point_mass");
  use_gravity_ = params->getScalar<int>("hydro_use_gravity");
  r_accrete_   = params->getScalar<double>("hydro_accrete_radius");

  grid->get_radial_edges(r_out_,r_min_,v_out_,v_min_);
  v_min_ = params->getScalar<double>("hydro_v_piston");


  // thermal bomb parameters
  double E_bomb = params->getScalar<double>("hydro_bomb_energy");
  double r_bomb = params->getScalar<double>("hydro_bomb_radius");
  // default if r_bomb undefined, put  in inner 5 zones
  if (r_bomb == 0)
  {
    if (nz_ <= 10)
      r_bomb = r_out_[1];
    else
      r_bomb = r_out_[10];
  }
  double vol_bomb = 4.0*pc::pi/3.0*r_bomb*r_bomb*r_bomb;
  double p_bomb = (E_bomb/vol_bomb/3.0);

  // make sure things are calculated
  for (int i = 0;i < nz_;i++)
  {
    double vol = grid->zone_volume(i);
    // radiation pressure EOS
    //	grid->z[i].p_gas = pc::a*pow(grid->z[i].T_gas,4)/3.0;
    // gas pressure EOS
    grid->z[i].p_gas = pc::k*grid->z[i].rho/pc::m_p*grid->z[i].T_gas;

    // add in bomb
    if (r_out_[i] <= 3*r_bomb)
    {
      grid->z[i].p_gas += p_bomb*exp(-r_out_[i]*r_out_[i]/r_bomb/r_bomb);
      grid->z[i].T_gas = pow(3.0*grid->z[i].p_gas/pc::a,0.25);
    }

  }

  // time to write out
  time_write_ = 0;
  // clear output file
  std::ofstream output;
  output.open("hydro_data.dat");
  output.close();

  */

}


double hydro_1D_movingmesh::get_time_step()
{
  double tstep = 0;
  return tstep;
}



void hydro_1D_movingmesh::step(double dt)
{


}
