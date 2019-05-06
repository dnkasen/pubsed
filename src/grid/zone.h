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
  real e_gas;           // gas pressure
  real p_gas;           // gas pressure
  real T_gas;           // gas temperature
  real grey_opacity;    // grey opacity (cm^2/g) if the user defines a zone-dependent grey opacity

  // composition of gas
  std::vector<real> X_gas;   // mass fractions of elements in zone
  real mu;                   // mean atomic mass

  // radiation quantities
  real e_rad;      // radiation energy density  (ergs/cm^3) in lab frame
  real e_abs;      // radiation energy deposition density rate (ergs/cm^3/s)
  real e_abs_compton;      // radiation energy deposition density rate (ergs/cm^3/s)
  real fx_rad;     // radiation x-force in lab frame
  real fy_rad;     // radiation y-force in lab frame
  real fz_rad;     // radiation z-force in lab frame
  real fr_rad;     // radiation radial force in lab frame
  real eps_imc;    // fleck factor effective absorption
  real L_thermal;  // thermal luminosity 

  // radioactive quantities
  real L_radio_emit;
  real L_radio_dep;

  // four force vector in lab frame
  real G1, G2, G3;

  // radiation pessure tensor components (symmetric)
  real P11, P12, P13, P22, P23, P33;

};

#endif
