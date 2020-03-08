#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "xy_array.h"

using namespace std;

//---------------------------------------------------------
// Just allocation the memory for this
//---------------------------------------------------------
void xy_array::init(int n) 
{
  x.resize(n);
  y.resize(n);
  for (int i=0;i<x.size();i++)  {
    x[i] = 0;
    y[i] = 0; }
      
}

//---------------------------------------------------------
// Initialize with start, stop and delta
//---------------------------------------------------------
void xy_array::init(double start, double stop, double delta)
{
  int n = (int)((stop-start)/delta) + 1;
  if (n <= 1) {n = 1; delta = stop-start;}
  if (delta < 0) delta = 0;
  init(n);

  for (int i=0;i<x.size();i++) 
    x[i]    = start + i*delta;

}

//---------------------------------------------------------
// locate
// If off the boundaries of the array, return the
// boundary value
//---------------------------------------------------------
int xy_array::locate(double z)
{
  // the degenerate case always returns 0
  if (x.size() == 1) return 0;
  
  // a form of locate from numerical recipes
  int bm;                             // mid point
  int bl = 0;                         // lower bound
  int bu = x.size()-1;               // upper bound

  // check if we are off the ends of the array
  if (z >= x[bu]) return (int)x.size();
  if (z <= x[bl]) return 0;
    
  // search the array for this index
  while (bu-bl > 1)
  {
    bm = (bu+bl)/2;
    if (x[bm] <= z) bl = bm;
    else bu = bm;
  }
  return bl;
} 


//---------------------------------------------------------
// Linear Interpolation of a passed array
//---------------------------------------------------------
double xy_array::value_at(double z)
{
  int ind = locate(z);
  
  double v,slope;
  if (ind < x.size()-1)
  {
    int i1    = ind;
    int i2    = ind + 1;
    slope = (y[i2]-y[i1])/(x[i2] - x[i1]);
    v     = y[ind] + slope*(z - x[i1]);
  }
  else
  {
    int i2    = ind;
    int i1    = ind - 1;
    slope = (y[i2]-y[i1])/(x[i2] - x[i1]);
    v     = y[ind] + slope*(z - x[i1]);
  }

  return v;
}



//---------------------------------------------------------
// Linear Interpolation of a passed array
// return zero if off array
//---------------------------------------------------------
double xy_array::value_at_with_zero_edges(double z)
{
  if (z < x[0]) return 0;
  if (z > x[x.size()-1]) return 0;
  else return value_at(z);
}


void xy_array::print()
{
  printf("# Print Locate Array; n_elements = %lu\n",x.size());
  for (int i=0;i<x.size();i++)
    printf("%4d %12.4e %12.4e\n",i,x[i],y[i]);
}
  
