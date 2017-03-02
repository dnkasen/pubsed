#include <algorithm>
#include <stdio.h>
#include "cdf_array.h"


//------------------------------------------------------
// return the actual y value, not the integrated
//------------------------------------------------------
double cdf_array::get_value(const int i) const
{
  if (i==0) return y[0];
  else return (y[i] - y[i-1]);  
}

//------------------------------------------------------
// set the actual y value, not the integrated
// must be called in order
//------------------------------------------------------
void cdf_array::set_value(const int i, const double f)   
{
  if (i==0) y[0] = f;
  else y[i] = y[i-1] + f;
}

//------------------------------------------------------
// Normalize such that the last entry is 1.0
//------------------------------------------------------
void cdf_array::normalize() 
{
  // check for zero array, set to all constant
  if (y.back() == 0) y.assign(y.size(),1.0);

  // normalize to end = 1.0
  N = y.back();
  for (int i=0;i<y.size();i++)   y[i] /= N;
}


//---------------------------------------------------------
// Sample the probability distribution using binary search.
// Pass a random number betwen 0 and 1.  
// Returns the index of the first value larger than yval
// if larger than largest element, returns size
//---------------------------------------------------------
int cdf_array::sample(const double yval) const
{
  if (y.size() == 1) return 0;
  return upper_bound(y.begin(), y.end(), yval) - y.begin();
}


//------------------------------------------------------
// Simple printout
//------------------------------------------------------
void cdf_array::print() const{
  for (int i=0;i<y.size();i++) 
    printf("%5d %10.4e %10.4e\n",i,get_value(i),y[i]);
}
  
//------------------------------------------------------
// Clear the arrays
//------------------------------------------------------
void cdf_array::wipe()
{
  y.assign(y.size(), 0.0);
}
  
//------------------------------------------------------------
// just returning the size of the array
//------------------------------------------------------------
int cdf_array::size() const
{
  return y.size();
}
