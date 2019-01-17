#ifndef _HYDRO_1D_LAGRANGIAN_H
#define _HYDRO_1D_LAGRANGIAN_H

#include <vector>
#include <fstream>
#include <iostream>
#include "zone.h"
#include "hydro_general.h"



class hydro_1D_lagrangian : public hydro_general
{


 public:

  void   init(ParameterReader*, grid_general*);
  double get_time_step();
  void   step(double dt);
  
  // zone properties
  int nz_;                     // number of zones
  int z_start_;                // start zone to do hydro on
  std::vector<double> r_out_;  // radius of outer zone edge
  std::vector<double> v_out_;  // velocity at outer zone edge
  std::vector<double> mass_;   // mass of zone
  std::vector<double> eden_;   // gas energy density per unit mass
  double r_min_,v_min_;        // radius and velocity at innermost edge
  double r_accrete_;           // inner boundary for things to get accreted
  double time_;
  double time_write_;


  // artificial viscosity
  std::vector<double> visq_;
  double C_q_;

  // gravity params
  int use_gravity_;
  double M_center_;

  // boundary conditions
  int boundary_outflow_;
  int boundary_rigid_outer_wall_;

  double get_dr(int i)
  {
    if (i == 0) return r_out_[0] - r_min_;
    else return r_out_[i] - r_out_[i-1];
  }

  double compute_artificial_viscosity(int i);
  
  //int update_temp;
 // double p_ext;
 // double alpha_imc;
//  void Update_Temperature(int);
//  void Remap_Zones();
//  void Solve_EOS_E(int);
//  void Solve_EOS_T(int);

};

#endif
