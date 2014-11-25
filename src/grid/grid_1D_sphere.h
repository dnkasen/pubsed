#ifndef _GRID_1D_SPHERE_H
#define _GRID_1D_SPHERE_H 1

#include <fstream>
#include <vector>
#include "grid_general.h"
#include "locate_array.h"


//*******************************************
// 1-Dimensional Spherical geometry
//*******************************************
class grid_1D_sphere: public grid_general
{

private:

  // specifics to this geometry
  double r_inner;

  // store location of the outer edge of the zone.
  locate_array r_out;

  // store volumes explicitly
  std::vector<double> vol;
  
  // functions for reading in files
  void read_SNR_file(std::ifstream &, int);


public:

  virtual ~grid_1D_sphere() {}

  void read_model_file(ParameterReader*);


  // required functions
  int     get_zone(const double *) const;
  double  zone_volume(const int) const;
  double  zone_min_length(const int) const;
  void    sample_in_zone(int, std::vector<double>, double[3]);
  void    velocity_vector(int i, double[3], double[3]);
  void    write_out(int);
  void    expand(double);

  void    coordinates(int i,double r[3]) {
    r[0] = r_out[i]; r[1] = 0; r[2] = 0;}

};


#endif
