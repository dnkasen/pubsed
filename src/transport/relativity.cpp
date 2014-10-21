#include <cassert>
#include "transport.h"
#include "particle.h"
#include "physical_constants.h"

namespace pc = physical_constants;

//---------------------------------------------------------
// lorentz factor ("gamma")
// v_rel = v_newframe - v_oldframe
//---------------------------------------------------------
double lorentz_factor(const double (&v_rel)[3])
{
  double beta2 = (v_rel[0]*v_rel[0] + v_rel[1]*v_rel[1] + v_rel[2]*v_rel[2]) / (pc::c*pc::c);
  return 1.0 / sqrt(1.0 - beta2);
}

//------------------------------------------------------------
// dot product of v_rel and relativistic particle direction
// v_rel = v_newframe - v_oldframe
// D = direction vector of relativistic particle in old frame
//------------------------------------------------------------
double v_dot_d(const double (&v_rel)[3], const double (&D)[3])
{
  return v_rel[0]*D[0] + v_rel[1]*D[1] + v_rel[2]*D[2];
}

//------------------------------------------------------------
// v_rel = v_newframe - v_oldframe
// v_dot_v is the dot product of the relative velocity 
// and the relativistic particle's direction
//------------------------------------------------------------
double doppler_shift(const double gamma, const double vdd)
{
  return gamma * (1.0 - vdd/pc::c);
}

//------------------------------------------------------------
// apply a lorentz transform to the particle
// v_rel = v_newframe - v_oldframe
//------------------------------------------------------------
void lorentz_transform(particle* p, const double (&v_rel)[3])
{
  // calculate the doppler shift, v dot D, and lorentz factors
  double gamma  = lorentz_factor(v_rel);
  double vdd    = v_dot_d(v_rel, p->D);
  double dshift = doppler_shift(gamma, vdd);

  // transform the 0th component (energy and frequency)
  p->e  *= dshift; 
  p->nu *= dshift;

  // transform the 1-3 components (direction)
  // See Mihalas & Mihalas eq 89.8
  p->D[0] = 1.0/dshift * (p->D[0] - gamma*v_rel[0]/pc::c * (1 - gamma*vdd/pc::c/(gamma+1)) );
  p->D[1] = 1.0/dshift * (p->D[1] - gamma*v_rel[1]/pc::c * (1 - gamma*vdd/pc::c/(gamma+1)) );
  p->D[2] = 1.0/dshift * (p->D[2] - gamma*v_rel[2]/pc::c * (1 - gamma*vdd/pc::c/(gamma+1)) );

  // for security, make sure it is properly normalized
  double norm = p->D[0]*p->D[0] + p->D[1]*p->D[1] + p->D[2]*p->D[2];
  p->D[0] /= norm;
  p->D[1] /= norm;
  p->D[2] /= norm;
}


//------------------------------------------------------------
// get the doppler shift when moving from frame_to_frame
// does not change any particle properties
//------------------------------------------------------------
double transport::dshift_comoving_to_lab(particle* p) const
{
  if(p->ind < 0){
    std::cout << p->r()-1e7<< std::endl;
    std::cout << p->ind<< std::endl;
  }
  assert(p->ind >= 0);
  double v[3];
  grid->velocity_vector(p->ind,p->x,v); // v_comoving - v_lab

  // new frame is lab frame. old frame is comoving frame.
  // v_rel = v_lab - v_comoving  --> v must flip sign.
  v[0] *= -1;
  v[1] *= -1;
  v[2] *= -1;

  double gamma = lorentz_factor(v);
  double vdd = v_dot_d(v, p->D);
  return doppler_shift(gamma, vdd);
}

double transport::dshift_lab_to_comoving(particle* p) const
{
  double v[3];
  grid->velocity_vector(p->ind,p->x,v); // v_comoving - v_lab

  // new frame is comoving frame. old frame is lab frame.
  // v_rel = v_comoving - v_lab  -->  v keeps its sign

  double gamma = lorentz_factor(v);
  double vdd = v_dot_d(v, p->D);
  return doppler_shift(gamma, vdd);
}


//------------------------------------------------------------
// do a lorentz transformation; modifies the energy, frequency
// and direction vector of the particle
//------------------------------------------------------------
void transport::transform_comoving_to_lab(particle* p) const
{
  double v[3];
  grid->velocity_vector(p->ind,p->x,v); // v_comoving - v_lab

  // new frame is lab frame. old frame is comoving frame.
  // v_rel = v_lab - v_comoving  --> v must flip sign.
  v[0] *= -1;
  v[1] *= -1;
  v[2] *= -1;

  lorentz_transform(p,v);
}

void transport::transform_lab_to_comoving(particle* p) const
{
  double v[3];
  grid->velocity_vector(p->ind,p->x,v); // v_comoving - v_lab

  // new frame is lab frame. old frame is comoving frame.
  // v_rel = v_comoving - v_lab  --> v keeps its sign.

  lorentz_transform(p,v);
}
