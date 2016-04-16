#ifndef _LOCATE_ARRAY_H
#define _LOCATE_ARRAY_H 

#include <vector>

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
  void init(const double,const double,const int);
  void init(const std::vector<double>, const double minval);
  void copy(locate_array l);
  void swap(locate_array new_array);

  // center of the bin left of nu_i
  double center(const int i) const{
    if (i == 0) return 0.5*(min    + x[0]);
    else        return 0.5*(x[i-1] + x[i]);
  }

  // width of the bin left of nu_i
  double delta(const int i) const{
    if (i == 0) return x[0] - min;
    else        return x[i] - x[i-1];
  }


  double maxval() {return x[x.size() - 1]; }
  double minval() {return min; }

  int    locate(const double) const;
  double interpolate_between(const double,const int,const int,const std::vector<double>&) const;
  double log_interpolate_between(const double,const int,const int,const std::vector<double>&) const;
  double sample(const int, const double) const;
  void   print() const;
  double value_at(const double nu, const std::vector<double>& array,int) const;
  double value_at(const double nu, const std::vector<double>& array) const;
  
  double value_at_extrapolate(const double nu, const std::vector<double>& array) const;

  // operators for easy access
  double  operator[] (const int i) const {return x[i];};
  double& operator[] (const int i)       {return x[i];};
  void resize(int i) {x.resize(i);};
};

#endif
