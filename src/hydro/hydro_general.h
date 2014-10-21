#ifndef _HYDRO_GENERAL_H
#define _HYDRO_GENERAL_H

#include "grid_general.h"
#include "Lua.h"

class hydro_general
{

 private:

 public:

  grid_general *grid;
  double t_now;
  
  virtual void init(Lua*, grid_general*) = 0;
  virtual double get_time_step() = 0;
  virtual void step(double)      = 0;


};

#endif
