#ifndef _RADIOACTIVE_H
#define _RADIOACTIVE_H 

#include <vector>

// radioactive decay constants

// Junde et al. (2011, NDS, 112, 1513)
// https://ui.adsabs.harvard.edu/#abs/2011NDS...112.1513J/abstract
// And ENSDF for average positron kinetic energy
#define TAU_NI  757241.8       // e-folding time = 8.764 days
#define TAU_CO  9627379.       // e-folding time = 111.428 days 
// average energy per decay in MeV
#define AVERAGE_NI_ENERGY 1.725 // 1.7206 (gamma+x-ray) + 0.0000 (pos KE) + 0.0043 (Auger+CE)
#define AVERAGE_CO_ENERGY 3.730 // 3.6052 (gamma+x-ray) + 0.1210 (pos KE) + 0.0037 (Auger+CE)
#define CO_POSITRON_FRACTION 0.0324

// Dong & Junde (2015, NDS, 128, 185)
// https://ui.adsabs.harvard.edu/#abs/2015NDS...128..185D/abstract
// And ENSDF for average positron kinetic energy
// Ignoring 1.8% branch of 52Mn metastable --> 52Mn ground
#define TAU_FE  42977.885        // e-folding time = 11.9383 hours
#define TAU_MN  1826.4519        // e-folding time = 30.4409 minute
// average energy per decay in MeV
#define AVERAGE_FE_ENERGY 0.924 // 0.7358 (gamma) + 0.1887 (pos KE) + ? (Auger+CE)
#define AVERAGE_MN_ENERGY 3.534 // 2.4016 (gamma+x-ray) + 1.1326 (pos KE) + 0.0001 (Auger+CE)

// Burrows (2006, NDS, 107, 1747)
// https://ui.adsabs.harvard.edu/#abs/2006NDS...107.1747B/abstract
// And ENSDF for average positron kinetic energy
#define TAU_CR  111976.        // e-folding time = 31.1045 hours  
#define TAU_VN  1.99107843e6   // e-folding time = 23.0449 days
// average energy per decay in MeV
#define AVERAGE_CR_ENERGY 0.442 // 0.4336 (gamma+x-ray) + 0.0015 (pos KE) + 0.0068 (Auger+CE)
#define AVERAGE_VN_ENERGY 3.066 // 2.9193 (gamma+x-ray) + 0.1450 (pos KE) + 0.0019 (Auger+CE)


class radioactive
{
 
public:

	double rprocess_heating_rate(double t, double *gfrac);

	void decay_composition
	(std::vector<int> elem_Z, std::vector<int> elem_A, std::vector<double>& X, double t);

	double decay(std::vector<int> elem_Z, std::vector<int> elem_A, 
	       std::vector<double> X, double t, double *gfrac, int force_rproc);

	double decay_energy_rate(int, int, double, double*);		
  //  double sample_particle_energy(int, int, double);
  //void   Fractions(int i, double *x, double t);

  //  static const double tau_56Ni =  760579.0;
  //static const double tau_56Co =  9.81599e06;   // decay time = 113.611 days 

  
};


#endif
