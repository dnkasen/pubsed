#ifndef PHYS_CONSTANTS
#define PHYS_CONSTANTS 1


namespace physical_constants
{

  // fundamental constants
  const double pi   =  3.14159;         // just pi
  const double c    =  2.99792458e10;   // speed of light (cm/s)
  const double h    =  6.6260755e-27;   // plancks constant (ergs s)
  const double k    =  1.380658e-16;    // boltz constatn (ergs/K)
  const double k_ev =  8.6173324e-5;    // boltzmann constant (ev/K)
  const double m_p  =  1.67262158e-24;  // proton mass (g)
  const double m_e  = 9.10938188e-28;   // mass of electron (g)
  const double sb   = 5.6704e-5;        // stefan-boltzmann (ergs cm^-2 s^-1 K^-4) 
  const double a    = 7.5657e-15;       // radiation constant (ergs cm-3 K-4)
  const double m_e_MeV =  0.510998910511;    // rest energy of electron in Mev
  const double thomson_cs = 0.66523e-24;       // Thomson cross-section cm-2 
  const double sigma_tot = 0.0265400193567;    // integrated line coefficent (cm^2 Hz)
  const double alpha_fs  = 0.007;              // fine structure constant
  const double rydberg   = 2.1798741e-11;     // rydberg constant in ergs


  // conversion factors
  const double ev_to_ergs   =  1.60217646e-12; 
  const double ergs_to_ev   =  624150975230.2815;
  const double Mev_to_ergs  =  1.60217646e-6;  
  const double cm_to_angs   =  1.0e8;
  const double angs_to_cm   =  1.0e-8;

  // astrophysical constants
  const double m_sun   = 1.98892e33;      // solar mass in grams

}


#endif
