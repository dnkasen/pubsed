#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cassert>
#include "mpi.h"
#include "Lua.h"
#include "grid_3D_cart.h"
#include "physical_constants.h"

//------------------------------------------------------------
// Read in a cartesian model file
//------------------------------------------------------------
void grid_3D_cart::read_model_file(Lua* lua)
{

}


//------------------------------------------------------------
// Overly simple search to find zone
//------------------------------------------------------------
int grid_3D_cart::get_zone(const double *x) const
{
  int i = floor((x[0]-x0)/dx);
  int j = floor((x[1]-y0)/dy);
  int k = floor((x[2]-z0)/dz);

  // check for off grid
  if ((i < 0)||(i > nx-1)) return -2;
  if ((j < 0)||(j > ny-1)) return -2;
  if ((k < 0)||(k > nz-1)) return -2;
  
  int ind =  i*ny*nz + j*nz + k;
  return ind;
}


//------------------------------------------------------------
// return volume of zone (precomputed)
//------------------------------------------------------------
double grid_3D_cart::zone_volume(const int i) const
{
  return vol;
}

//------------------------------------------------------------
// sample a random position within the cubical cell
//------------------------------------------------------------
void grid_3D_cart::sample_in_zone
(const int i, const std::vector<double> ran,double r[3])
{
  r[0] = x0 + (ix[i] + ran[0])*dx;
  r[1] = y0 + (iy[i] + ran[1])*dy;
  r[2] = z0 + (iz[i] + ran[2])*dz;
}



/------------------------------------------------------------
// return length of zone
//------------------------------------------------------------
double  grid_3D_cart::zone_min_length(const int i) const
{
  return min_ds;
}



//------------------------------------------------------------
// get the velocity vector 
//------------------------------------------------------------
void grid_3D_cart::velocity_vector(const int i, const double x[3], double v[3]) const
{
  // may want to interpolate here
  v[0] = z[i].v[0];
  v[1] = z[i].v[1];
  v[2] = z[i].v[2];
}

//------------------------------------------------------------
// cell-centered coordinates of zone i
//------------------------------------------------------------
void grid_3D_cart::coordinates(const int i,double r[3]) const
{
  r[0] = x0 + (ix[i]+0.5)*dx;
  r[1] = y0 + (iy[i]+0.5)*dy;
  r[2] = z0 + (iz[i]+0.5)*dz;
}

