#ifndef _CDF_H
#define _CDF_H 1

#include <vector>
#include <algorithm>


//**********************************************************
// CDF == Comulative Distribution Function
//
// This simple class just holds a vector which should be
// monitonically increasing and reaches unity
// We can sample from it using a binary search.
// the CDF value at locate_array's "min" is assumed to be 0
//**********************************************************

template < class T> class cdf_array
{

private:
  
  std::vector<T> y;
  
public:

  void resize(const int n)  {y.resize(n); }


  //------------------------------------------------------
  // return the CDF value at index i
  //------------------------------------------------------
  T get(const int i) const 
    {return y[i];}   
  
  //------------------------------------------------------
  // set the CDF value at index i
  //------------------------------------------------------
  void   set(const int i, T f)  {y[i] = f;}      

  //------------------------------------------------------
  // return the actual y value, not the integrated CDF
  //------------------------------------------------------
  T get_value(const int i) const
  {
    if (i==0) return y[0];
    else return (y[i] - y[i-1]);  
  }

  //------------------------------------------------------
  // set the actual y value, not the integrated
  // must be called in order
  //------------------------------------------------------
  void set_value(const int i, T f)   
  {
  if (i==0) y[0] = f;
  else y[i] = y[i-1] + f;
  }

//------------------------------------------------------
// Normalize such that the last entry is 1.0
//------------------------------------------------------
void normalize() 
{

  // check for nan
  for (int i=0;i<y.size();i++) if (std::isnan(y[i])) y[i] = 0;

  // check for zero array, set to all constant
  if (y.back() == 0) y.assign(y.size(),1.0);

  // normalize to end = 1.0
  double N = y.back();
  for (int i=0;i<y.size();i++)   y[i] /= N;

}


//---------------------------------------------------------
// Sample the probability distribution using binary search.
// Pass a random number betwen 0 and 1.  
// Returns the index of the first value larger than yval
// if larger than largest element, returns size
//---------------------------------------------------------
int sample(const double yval) const
{
  if (y.size() == 1) return 0;
  int v = upper_bound(y.begin(), y.end(), yval) - y.begin();
  return v;
}


//------------------------------------------------------
// Simple printout
//------------------------------------------------------
void print() const{
  for (int i=0;i<y.size();i++) 
    printf("%5d %10.4e %10.4e\n",i,get_value(i),y[i]);
}
  
//------------------------------------------------------
// Clear the arrays
//------------------------------------------------------
void wipe()
{
  y.assign(y.size(), 0.0);
}
  
//------------------------------------------------------------
// just returning the size of the array
//------------------------------------------------------------
int size() const
{
  return y.size();
}


};

#endif
