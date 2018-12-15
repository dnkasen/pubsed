#ifndef _GRID_3D_CART_H
#define _GRID_3D_CART_H 1

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

  int    nx_, ny_, nz_; // number of zones in each dimension
  double dx_, dy_, dz_; // length of each zone in each dimension
  double x0_, y0_, z0_; // leftmost points
  double vol_;        // volume of each zone = dx*dy*dz
  int *index_x_,*index_y_,*index_z_;

  int get_index(int i, int j, int k)
  {
    int ind =  i*ny_*nz_ + j*nz_ + k;
    return ind;
  }

public:

  virtual ~grid_3D_cart() {}

  // required functions
  void    read_model_file(ParameterReader*);
  void    write_plotfile(int,double,int);
  int     get_zone(const double *) const;
  double  zone_volume(const int) const;
  double  zone_min_length(const int) const;
  void    sample_in_zone(int, std::vector<double>, double[3]);
  void    get_velocity(int i, double[3], double[3], double[3], double*);
  void    expand(double);
  int     get_next_zone(const double *x, const double *D, int, double, double *dist) const;
  void    coordinates(int i,double r[3]);

};


#endif
