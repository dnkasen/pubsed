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
  SedonaReal v[3];            // velocity vector (cm/s)
  SedonaReal rho;             // density (g/cm^3)
  SedonaReal cs;              // sound speed (cm/s)
  SedonaReal e_gas;           // gas pressure
  SedonaReal p_gas;           // gas pressure
  SedonaReal T_gas;           // gas temperature
  SedonaReal n_elec;          // number of free electrons
  SedonaReal bulk_grey_opacity;           // bulk component of the grey opacity (cm^2/g), which is the same in every zone
  SedonaReal zone_specific_grey_opacity;  // zone-specific component of the grey opacity (cm^2/g), which varies from zone to zone
  SedonaReal total_grey_opacity;          // total grey opacity (cm^2/g), which is the sum of the bulk component and the zone-specific component

  // composition of gas
  std::vector<SedonaReal> X_gas;   // mass fractions of elements in zone
  SedonaReal mu_I;                   // mean atomic/ionic mass (not including free electrons). Dimensionless; needs to be multiplied by amu (~ m_p) to get units of grams

  // radiation quantities
  SedonaReal e_rad;      // radiation energy density  (ergs/cm^3) in lab frame
  SedonaReal e_abs;      // radiation energy deposition density rate (ergs/cm^3/s)
  SedonaReal fx_rad;     // radiation x-force in lab frame
  SedonaReal fy_rad;     // radiation y-force in lab frame
  SedonaReal fz_rad;     // radiation z-force in lab frame
  SedonaReal fr_rad;     // radiation radial force in lab frame
  SedonaReal eps_imc;    // fleck factor effective absorption
  SedonaReal L_thermal;  // thermal luminosity

  // radioactive quantities
  SedonaReal L_radio_emit;
  SedonaReal L_radio_dep;

  // four force vector in lab frame
  SedonaReal G1, G2, G3;

  // radiation pessure tensor components (symmetric)
  SedonaReal P11, P12, P13, P22, P23, P33;

};

#endif
