#ifndef _THREAD_RNG_H
#define _THREAD_RNG_H
#include <gsl/gsl_rng.h>
#include <vector>

class thread_RNG
{

protected:

  // OPTIMIZE - to prevent false sharing, make sure all the RNG's are on different cache lines
  // vector of generators
  std::vector<gsl_rng*> generators;
  
public:

  void   init();
  double uniform();
};

#endif
