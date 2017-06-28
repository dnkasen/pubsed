#ifndef _GRID_2D_SPHERE_H
#define _GRID_2D_SPHERE_H 1

#include <fstream>
#include <vector>
#include "grid_general.h"
#include "locate_array.h"


//*******************************************
// 2-Dimensional Cylndrical geometry
//*******************************************
class grid_2D_cyln: public grid_general
{

private:


  int    nx_, nz_; // number of zones in each dimension
  double dx_, dz_; // length of each zone in each dimension
  double zcen_   ; // center coordinate of the z-axis

  std::vector<int> index_x_; // map to x index from 1D index
  std::vector<int> index_z_; // map to z index from 1D index

  // store precomputed zone volumes
  std::vector<double> vol_; 

public:

  virtual ~grid_2D_cyln() {}


  // required functions
  void    read_model_file(ParameterReader*);
  void    write_out(int,double);
  int     get_zone(const double *) const;
  double  zone_volume(const int)   const;
  void    sample_in_zone(int, std::vector<double>, double[3]);
  void    get_velocity(int i, double[3], double[3], double[3], double*);
  void    expand(double);
  int     get_next_zone(const double *x, const double *D, int, double, double *dist) const;
  void    coordinates(int i,double r[3]) {
    r[0] = 0; r[1] = 0; r[2] = 0;}

};


#endif
