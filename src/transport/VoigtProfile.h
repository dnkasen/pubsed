#ifndef _VOIGT_PROFILE_H
#define _VOIGT_PROFILE_H

#include <gsl/gsl_rng.h>


class VoigtProfile
{

 private:
  
  double u0_, eu0_;
  gsl_rng  *rangen_;    // random number generator

public:
  
  VoigtProfile();
  void   setU0(double);
  double getProfile(double,double);
  double sampleU(double,double);

};




#endif

