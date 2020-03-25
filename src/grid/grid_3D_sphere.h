#ifndef _GRID_3D_SPHERE_H
#define _GRID_3D_SPHERE_H 1

#include <fstream>
#include <vector>
#include "grid_general.h"
#include "locate_array.h"


//*******************************************
// 3-Dimensional Spherical geometry
//*******************************************
class grid_3D_sphere: public grid_general
{

private:

  // Flags to remember how grid was read in
  herr_t status_dr_;
  herr_t status_rmin_;
  herr_t status_rthetaphi_;

  int nr_, ntheta_, nphi_; // number of zones in each dimension

  // the right zone edge in each direction
  // the lengths of these arrays are nx_, ny_, and nz_, respectively
  // these arrays are indexed by the x-index, y-index, and z-index, respectively
  locate_array r_out_, theta_out_, phi_out_;

  // store precomputed zone widths in each direction. These arrays are indexed by the x-index, y-index, and z-index, respectively
  std::vector<double> dr_, dtheta_, dphi_;

  // store precomputed zone volumes. These arrays are indexed by the index in the flattened 1D array of all zones
  std::vector<double> vol_;

  std::vector<int> index_r_; // map to x index from the index in the flattened 1D array of all zones
  std::vector<int> index_theta_; // map to y index from the index in the flattened 1D array of all zones
  std::vector<int> index_phi_; // map to z index from the index in the flattened 1D array of all zones

  int get_index(int i, int j, int k)
  {
    int ind =  i*ntheta_*nphi_ + j*nphi_ + k;
    return ind;
  }

public:

  virtual ~grid_3D_sphere() {}

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
  void    coordinates(int i,double r[3]) {r[0] = 0; r[1] = 0; r[2] = 0;}

};


#endif