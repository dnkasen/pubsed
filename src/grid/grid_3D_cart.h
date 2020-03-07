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

  int nx_, ny_, nz_; // number of zones in each dimension

  // the right zone edge in each direction
  // the lengths of these arrays are nx_, ny_, and nz_, respectively
  // these arrays are indexed by the x-index, y-index, and z-index, respectively
  locate_array x_out_, y_out_, z_out_;

  // store precomputed zone widths in each direction. These arrays are indexed by the x-index, y-index, and z-index, respectively
  std::vector<double> dx_, dy_, dz_;

  // store precomputed zone volumes. These arrays are indexed by the index in the flattened 1D array of all zones
  std::vector<double> vol_;

  std::vector<int> index_x_; // map to x index from the index in the flattened 1D array of all zones
  std::vector<int> index_y_; // map to y index from the index in the flattened 1D array of all zones
  std::vector<int> index_z_; // map to z index from the index in the flattened 1D array of all zones

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
