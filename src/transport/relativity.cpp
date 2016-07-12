#include <cassert>
#include "transport.h"
#include "particle.h"
#include "physical_constants.h"

namespace pc = physical_constants;

//------------------------------------------------------------
// get the doppler shift when moving from frame_to_frame
// values saved inside particle class
// tolab = 0 for lab_to_comoving
// tolab = 1 for comoving_to_lab (flips sign)
//------------------------------------------------------------
 double transport::do_dshift(particle* p, int tolab) 
 {
  assert(p->ind >= 0);

  // get velocity information here
  double v_rel[3], dvds;
  grid->get_velocity(p->ind,p->x,p->D,v_rel,&dvds); 

  // if new frame is lab frame. old frame is comoving frame.
  // v_rel = v_lab - v_comoving  --> v must flip sign.
  if (tolab)
  {
    v_rel[0] *= -1;
    v_rel[1] *= -1;
    v_rel[2] *= -1;
  }

  // get relativistic quantities
  double beta2 = (v_rel[0]*v_rel[0] + v_rel[1]*v_rel[1] + v_rel[2]*v_rel[2])/pc::c/pc::c;
  double vdd   = v_rel[0]*p->D[0] + v_rel[1]*p->D[1] + v_rel[2]*p->D[2];
  double gamma = 1.0/sqrt(1 - beta2);
  double dshift = gamma*(1 - vdd/pc::c);

  // store quantities
  p->gamma  = gamma;
  p->dshift = dshift;
  p->dvds   = dvds;

  return dshift;
 }

//------------------------------------------------------------
// doppler shift calls
//------------------------------------------------------------
double transport::dshift_lab_to_comoving(particle* p) 
{
  return do_dshift(p,0);
}

double transport::dshift_comoving_to_lab(particle* p) 
{
  return do_dshift(p,1);
}



//------------------------------------------------------------
// do a lorentz transformation; modifies the energy, frequency
// and direction vector of the particle
//------------------------------------------------------------
void transport::transform_lab_to_comoving(particle* p)
{
  assert(p->ind >= 0);
  lorentz_transform(p,0);
}

//------------------------------------------------------------
// do a lorentz transformation; modifies the energy, frequency
// and direction vector of the particle
//------------------------------------------------------------
void transport::transform_comoving_to_lab(particle* p)
{
  lorentz_transform(p,1);
}


//------------------------------------------------------------
// apply a lorentz transform to the particle
// v_rel = v_newframe - v_oldframe
//------------------------------------------------------------
void transport::lorentz_transform(particle* p, int tolab)
{
  assert(p->ind >= 0);

  // get velocity information here
  double v_rel[3],dvds;
  grid->get_velocity(p->ind,p->x,p->D,v_rel,&dvds); 

  // new frame is lab frame. old frame is comoving frame.
  // v_rel = v_lab - v_comoving  --> v must flip sign.
  if (tolab)
  {
    v_rel[0] *= -1;
    v_rel[1] *= -1;
    v_rel[2] *= -1;
  }

  // get relativistic quantities
  double beta2 = (v_rel[0]*v_rel[0] + v_rel[1]*v_rel[1] + v_rel[2]*v_rel[2])/pc::c/pc::c;
  double vdd   = v_rel[0]*p->D[0] + v_rel[1]*p->D[1] + v_rel[2]*p->D[2];
  double gamma = 1.0/sqrt(1 - beta2);
  double dshift = gamma*(1 - vdd/pc::c);

  // store quantities
  p->gamma  = gamma;
  p->dshift = dshift;
  p->dvds   = dvds;

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





