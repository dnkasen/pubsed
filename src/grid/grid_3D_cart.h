#ifndef _GRID_1D_SPHERE_H
#define _GRID_1D_SPHERE_H 1

#include <fstream>
#include <vector>
#include "grid_general.h"
#include "locate_array.h"


//*******************************************
// 1-Dimensional Spherical geometry
//*******************************************
class grid_3D_cart: public grid_general
{

private:

  // specifics to this geometry
 // specifics to this geometry

  int    nx, ny, nz; // number of zones in each dimension
  double dx, dy, dz; // length of each zone in each dimension
  double x0, y0, z0; // leftmost points
  double vol;        // volume of each zone = dx*dy*dz
  double min_ds;
  int *ix,*iy,*iz;
  int reflect_x, reflect_y, reflect_z;


public:

  virtual ~grid_3D_cart() {}

  void read_model_file(Lua *lua);


  // required functions
  int     get_zone(const double *) const;
  double  zone_volume(const int) const;
  double  zone_min_length(const int) const;
  void    sample_in_zone(int, std::vector<double>, double[3]);
  void    velocity_vector(const int i, const double[3], double[3]);
  void    write_out(int);
  void    expand(double);

};


#endif
