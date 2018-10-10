#ifndef _HOMOLOGOUS_HYDRO_H
#define _HOMOLOGOUS_HYDRO_H

#include "hydro_general.h"
#include "zone.h"

//*******************************************
// 1-Dimensional Spherical geometry
//*******************************************
class hydro_homologous: public hydro_general
{

 public:

  void   init(ParameterReader*, grid_general*);
  double get_time_step();
  void   step(double dt);
  void   evolve_to_start(double t_start, int force_rproc);


};

#endif
