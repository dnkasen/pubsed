#ifndef _ZONE_H
#define _ZONE_H
#include <vector>
#include "sedona.h"

//-------------------------------------------------
// Class to store properties of one zone
//-------------------------------------------------
class zone
{

public:

  // fluid properties
  real v[3];            // velocity vector (cm/s)
  real rho;             // density (g/cm^3)
  real cs;              // sound speed (cm/s)
  real p_gas;           // gas pressure
  real e_gas;           // gas energy density per gram
  real E_gas;           // gas total energy
  real T_gas;           // gas temperature

  // composition of gas
  std::vector<real> X_gas;   // mass fractions of elements in zone
  real mu;                   // mean atomic mass

  // radiation quantities
  real e_rad;      // radiation energy density  (ergs/cm^3) in lab frame
  real e_abs;      // radiation energy deposition density rate (ergs/cm^3/s)
  real fx_rad;     // radiation x-force in lab frame
  real fy_rad;     // radiation y-force in lab frame
  real fz_rad;     // radiation z-force in lab frame
  real eps_imc;    // fleck factor effective absorption

  // radioactive quantities
  real L_radio_emit;
  real L_radio_dep;

  // four force vector in lab frame
  real G1, G2, G3;

  // radiation pessure tensor components (symmetric)
  real P11, P12, P13, P22, P23, P33;

};

#endif
