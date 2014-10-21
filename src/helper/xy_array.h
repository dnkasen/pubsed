#ifndef _XY_ARRAY_H
#define _XY_ARRAY_H 1

#include <vector>

class xy_array {
 


public:

  std::vector<double> x;
  std::vector<double> y;

  // constructors
  xy_array()  {}
  xy_array(int n) {init(n);}
  xy_array(double s1,double s2,double sd) {init(s1,s2,sd);}

  // Return
  int    size()        {return (int)x.size();}

  void init(int);
  void init(double,double,double);

  int    locate(double);
  double value_at(double);
  void   print();
  
};

#endif
