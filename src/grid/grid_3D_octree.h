#ifndef _GRID_3D_OCTREE_H
#define _GRID_3D_OCTREE_H 1

#include <fstream>
#include <vector>
#include "grid_general.h"
#include "locate_array.h"


//*******************************************
// 3-Dimensional Octree AMR
//*******************************************
class grid_3D_octree: public grid_general
{

private:

  // specifics to this geometry

  long n_zones;

  double *width_; // width of zone
  double *xmin_, *ymin_, *zmin_; // coordinates of corner of cell
  long *parent_cell_id_;


public:

  // required functions
  void    read_model_file(ParameterReader*);
  void    write_plotfile(int,double);
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
