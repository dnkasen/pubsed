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

  // store location of the outer edge of the zone.
  locate_array r_out;
  // velocity at inner boundary
  double v_inner_; 

  // store volumes explicitly
  std::vector<double> vol;
  
  // functions for reading in files
  void read_SNR_file(std::ifstream &, int, int, double);


public:

  virtual ~grid_1D_sphere() {}

  void read_model_file(ParameterReader*);

  virtual void get_radial_edges
     (std::vector<double>&, double&, std::vector<double>&, double&) const;
  virtual void set_radial_edges
    (const std::vector<double>, const double, const std::vector<double>, const double);


  // required functions
  int     get_zone(const double *) const;
  double  zone_volume(const int)   const;
  void    sample_in_zone(int, std::vector<double>, double[3]);
  void    get_velocity(int i, double[3], double[3], double[3], double*);
  void    write_out(int,double);
  void    expand(double);

  int get_next_zone(const double *x, const double *D, int, double, double *dist) const;

  void    coordinates(int i,double r[3]) {
    r[0] = r_out[i]; r[1] = 0; r[2] = 0;}

};


#endif
