//------------------------------------------------------------------
//*****************************************************************
//*************************  GRID ********************************
//*****************************************************************
// The grid class is a construct whose main purpose is to handle
// the geometry of the model.  It does two main things: (1) Reads
// in the input density,temperature,composition files (that of
// course must have a specific geometry). (2) Given a set of
// 3-d coordinates, it will give the corosponding zone index (or
// note that the coords are off the grid).
//
// The grid holds an array of zones, where key fluid data is stored
//
// The grid class is an abstract class that will be used to
// create subclasses (e.g. grid_1D_sphere, grid_3D_cart, etc...
//*****************************************************************

#ifndef _GRID_GENERAL_H
#define _GRID_GENERAL_H 1

#include <string>
#include <iostream>


#include "sedona.h"
#include "zone.h"
#include "ParameterReader.h"

#include "hdf5.h"
#include "hdf5_hl.h"

class grid_general
{

 protected:

  // fill the grid with data from a model file
  virtual void read_model_file(ParameterReader*) = 0;

 public:

  // set everything up
  void init(ParameterReader* params);

  // descrption of grid type
  std::string grid_type;

  // vector of zones
  std::vector<zone> z;
  int n_zones;

  double t_now;

  // vector of isotopes to use
  std::vector<int> elems_Z;
  std::vector<int> elems_A;
  int n_elems;

  // mpi reduce quantities
  void reduce_radiation();
  void reduce_radiation_block(int, int);
  void reduce_T_gas();


  // output functions
  void write_integrated_quantities(int iw, double tt);
  void write_hdf5_plotfile_zones(hid_t file_id, hsize_t *dims_g, int ndims, double);

  //****** virtual functions (geometry specific)

  // get zone index from x,y,z position
  virtual int get_zone(const double *) const   = 0;

  // get zone index from x,y,z position
  virtual int get_next_zone(const double *, const double *, int, double, double *) const = 0;

  // return volume of zone i
  virtual double zone_volume(const int i) const         = 0;

  // randomly sample a position within the zone i
  virtual void sample_in_zone(int,std::vector<double>,double[3]) = 0;

  // give the velocity vector at this point in zone i
  virtual void get_velocity(int i, double[3], double[3], double[3], double*) = 0;

  // get the coordinates at the center of the zone i
  virtual void coordinates(int i,double r[3]) = 0;

  // write out the grid state
  virtual void write_plotfile(int,double,int) = 0;

  // expand the grid by a factor of e
  virtual void expand(double) = 0;

  //****** empty functions to overide
  virtual void get_radial_edges
  (std::vector<double>&, double&, std::vector<double>&, double&) const
  { }

  virtual void set_radial_edges
  (const std::vector<double>, double, const std::vector<double>, double)
  { }

  virtual void get_zone_size(int, double*)
  {
  }


};


#endif
