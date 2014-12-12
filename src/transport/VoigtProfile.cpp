#include <time.h>
#include <math.h>

#include "VoigtProfile.h"


VoigtProfile::VoigtProfile()
{
  u0_  = 0;
  eu0_ = 1;

  // setup and seed random number generator
  const gsl_rng_type * TypeR;
  gsl_rng_env_setup();
  gsl_rng_default_seed = (unsigned int)time(NULL);
  TypeR = gsl_rng_default;
  rangen_ = gsl_rng_alloc (TypeR);
}


void VoigtProfile::setU0(double upass)
{
  u0_     = upass;
  eu0_    = exp(-u0_*u0_);
}


double VoigtProfile::getProfile(double x, double a)
{
  double sqrt_pi = 1.77245385091;
  double pi = 3.14159265359;

  double xsq = x*x;
  double c   = (xsq - 0.855)/(xsq + 3.42);

  double q;
  if (c < 0) q = 0;
  else 
  {
    double pic = 5.674*c*c*c*c -9.207*c*c*c + 4.421*c*c + 0.1117*c;
    q = (1 + 21/xsq)*a/pi/(xsq + 1)*pic; 
  }
 
  double H = q*sqrt_pi + exp(-xsq);
  return H;
}



double VoigtProfile::sampleU(double x, double a)
{
  double pi = 3.14159265359;
  double u,th;

  double sgn = 1;
  if (x < 0) {sgn = -1; x = -1*x;}

  double theta0 = atan((u0_ - x)/a);
  double denom  = (1 - eu0_)*theta0 + (1 + eu0_)*pi/2;
  double p0 = (theta0 + pi/2)/denom;
  
  int stop = 0;
  while (!stop) 
  {
    double r1 = gsl_rng_uniform(rangen_);
    double r2 = gsl_rng_uniform(rangen_);
    double r3 = gsl_rng_uniform(rangen_);

    if (r1 < p0) 
    {
      th = -1*pi/2.0 + r2*(theta0 + pi/2);
      u = a*tan(th) + x; 
      if (r3 < exp(-u*u)) stop = 1;
    }
    else 
    {
      th = theta0 + r2*(pi/2 - theta0);
      u = a*tan(th) + x; 
      if (r3 < exp(-u*u)/eu0_) stop = 1;
    }
  }

  return sgn*u;
}
