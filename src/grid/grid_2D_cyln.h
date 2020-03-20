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

  // Flags to remember how grid was read in
  herr_t status_dr_;
  herr_t status_rmin_;
  herr_t status_xz_;

  int nx_, nz_; // number of zones in each dimension

  // the right zone edge in each direction
  // the lengths of these arrays are nx_ and nz_, respectively
  // these arrays are indexed by the x-index and z-index, respectively
  locate_array x_out_, z_out_;

  // store precomputed zone widths in each direction. These arrays are indexed by the x-index and z-index, respectively
  std::vector<double> dx_, dz_;

  // store precomputed zone volumes. These arrays are indexed by the index in the flattened 1D array of all zones
  std::vector<double> vol_;

  std::vector<int> index_x_; // map to x index from 1D index
  std::vector<int> index_z_; // map to z index from 1D index

  /* For testing */
  int    nx_new_, nz_new_;

  locate_array x_out_new_, z_out_new_;

  std::vector<double> dx_new_, dz_new_;

  std::vector<int> index_x_new_;
  std::vector<int> index_z_new_;

  std::vector<double> vol_new_; 
public:

  // required functions
  void    read_model_file(ParameterReader*);
  void    write_plotfile(int,double,int);
  int     get_zone(const double *) const;
  double  zone_volume(const int)   const;
  void    sample_in_zone(int, std::vector<double>, double[3]);
  void    get_velocity(int i, double[3], double[3], double[3], double*);
  void    expand(double);
  int     get_next_zone(const double *x, const double *D, int, double, double *dist) const;
  void    coordinates(int i,double r[3]) {r[0] = 0; r[1] = 0; r[2] = 0;}

  void writeCheckpointGrid(std::string fname);
  void readCheckpointGrid(std::string fname, bool test=false);
  void testCheckpointGrid(std::string fname);

  void restartGrid(ParameterReader* params);
};


#endif
