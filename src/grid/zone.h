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
  real n_elec;          // number of free electrons
  real bulk_grey_opacity;           // bulk component of the grey opacity (cm^2/g), which is the same in every zone
  real zone_specific_grey_opacity;  // zone-specific component of the grey opacity (cm^2/g), which varies from zone to zone
  real total_grey_opacity;          // total grey opacity (cm^2/g), which is the sum of the bulk component and the zone-specific component

  // composition of gas
  std::vector<real> X_gas;   // mass fractions of elements in zone
  real mu_I;                   // mean atomic/ionic mass (not including free electrons). Dimensionless; needs to be multiplied by amu (~ m_p) to get units of grams

  // radiation quantities
  real e_rad;      // radiation energy density  (ergs/cm^3) in lab frame
  real e_abs;      // radiation energy deposition density rate (ergs/cm^3/s)
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
