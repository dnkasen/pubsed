#ifndef _LOCATE_ARRAY_H
#define _LOCATE_ARRAY_H

#include <vector>
#include <math.h>
#include <hdf5.h>
#include <limits>

#include "h5utils.h"

class locate_array {

public:

  // where applicable, these values are the right bin wall (i.e. not the left)
  std::vector<double> x;
  double min;

  // other parameters
  int do_log_interpolate;

  // constructors
  locate_array()  {}
  locate_array(int n) {init(n);}

  // Return size of array (also, # of bins)
  int size() const {return (int)x.size();}

  void init(const int);
  void init(const double,const double,const double);
  void log_init(const double,const double,const double);
  void init(const double,const double,const int);
  void init(const std::vector<double>, const double minval);
  void copy(locate_array l);
  void swap(locate_array new_array);

  // operators for easy access
  double  operator[] (const int i) const {return x[i];};
  double& operator[] (const int i)       {return x[i];};
  void resize(int i) {x.resize(i);};

  // equality
  bool is_equal(locate_array l, bool complain);

  // center of the bin left of nu_i
  double center(const int i) const{
    if (i == 0) return 0.5*(min    + x[0]);
    else        return 0.5*(x[i-1] + x[i]);
  }

  // left side of bin i
  double left(const int i) const{
    if (i == 0) return min;
    else return x[i-1];}

  // right side of bin i
  double right(const int i) const{
    return x[i];  }

  // width of the bin left of nu_i
  double delta(const int i) const{
    if (i == 0) return x[0] - min;
    else        return x[i] - x[i-1];
  }


  double maxval() {return x[x.size() - 1]; }
  double minval() {return min; }

  int    locate(const double) const;
  int    locate_within_bounds(const double xval) const;

  double sample(const int, const double) const;
  void   print() const;

  void writeCheckpoint(std::string fname, std::string gname, std::string dset);
  void readCheckpoint(std::string fname, std::string gname, std::string dset);

  //---------------------------------------------------------
// Linear Interpolation of a passed array
//---------------------------------------------------------
template<typename T>
T interpolate_between(const double xval, const int i1, const int i2, const std::vector<T>& y) const
{
  if (x.size() == 1) return y[0];
  double slope = (y[i2]-y[i1]) / (x[i2]-x[i1]);
  double yval = y[i1] + slope*(xval - x[i1]);
  return yval;
}


//---------------------------------------------------------
// Log-Log Interpolation of a passed array
//---------------------------------------------------------
template<typename T>
T log_interpolate_between(const double xval, const int i1, const int i2, const std::vector<T>& y) const
{
  if (x.size() == 1) return y[0];

  // safeguard against equal opacities
  if(y[i1]==y[i2]) return y[i1];

  // safeguard against nonsensical values
  if(y[i1]<=0 || y[i2]<=0) return interpolate_between(xval, i1, i2, y);

  // do logarithmic interpolation
  double slope = log(y[i2]/y[i1]) / log(x[i2]/x[i1]);
  double logyval = log(y[i1]) + slope*log(xval/x[i1]);
  return exp(logyval);
}

//---------------------------------------------------------
// find the value of y at the locate_array's value of xval
// assumes 1-1 correspondence between y and locate_array
// will extrapolate to regions off either end of the
// array grid
//---------------------------------------------------------
template<typename T>
T value_at_extrapolate(const double xval, const std::vector<T>& y) const{

  int ind = locate(xval);
  return y[ind];
  /*
  int i1, i2;

  if (x.size() == 1) return y[0];

  if(ind == 0){                // If off left side of grid
    i1 = 0;
    i2 = 1;
  }
  else if(ind < x.size()){    // If within expected region of grid
    i1 = ind - 1;
    i2 = ind;
  }
  else{ //if(ind == x.size()) // If off the right side of the grid
    i1 = x.size() - 2;
    i2 = x.size() - 1;
  }

  if(do_log_interpolate) return log_interpolate_between(xval, i1, i2, y);
  else                   return     interpolate_between(xval, i1, i2, y);
  */
}


//---------------------------------------------------------
// find the value of y at the locate_array's value of xval
// assumes 1-1 correspondence between y and locate_array
// will not extrapolate off of array -- if passed xval
// is out of bounds, will just return the end values
//
// overloaded so that doesn't need to do the locate if
// already known and passed.
//---------------------------------------------------------
template<typename T>
T value_at(const double xval, const std::vector<T>& y) const
{
  int ind = locate(xval);
  return value_at(xval,y,ind);
}

template<typename T>
T value_at(const double xval, const std::vector<T>& y,int ind) const
{
  return y[ind];
  /*
  int i1, i2;

  if (x.size() == 1) return y[0];

  if(ind == 0){                // If off left side of grid
    return y[0];
  }
  else if(ind < x.size()){    // If within expected region of grid
    i1 = ind - 1;
    i2 = ind;
  }
  else{ //if(ind == x.size()) // If off the right side of the grid
    return y[x.size()-1];
  }

  if(do_log_interpolate) return log_interpolate_between(xval, i1, i2, y);
  else                   return     interpolate_between(xval, i1, i2, y);
  */
}

};
#endif
