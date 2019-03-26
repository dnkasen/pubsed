#ifndef _THREAD_RNG_H
#define _THREAD_RNG_H
#include <gsl/gsl_rng.h>
#include <vector>
#include <string>

class thread_RNG
{

protected:

  // OPTIMIZE - to prevent false sharing, make sure all the RNG's are on different cache lines
  // vector of generators
  std::vector<gsl_rng*> generators;
  std::vector<int> counters;
  
public:

  void   init(bool fix_seed = false, unsigned long int fixed_seed_val = 0);
  double uniform();
  ~thread_RNG();

  void writeCheckpointRNG(std::string fname);
  int readCheckpointRNG(std::string fname);
};

#endif
