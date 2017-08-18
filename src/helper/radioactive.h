#ifndef _RADIOACTIVE_H
#define _RADIOACTIVE_H 

#include <vector>

// radioactive decay constonstants

// 56Ni chain data can be found
// http://adsabs.harvard.edu/abs/1994ApJS...92..527N
#define TAU_NI  757241.7         // decay time = 8.803   days
#define TAU_CO  9627378.       // decay time = 113.611 days 
// average energy per decay in MeV
#define AVERAGE_NI_ENERGY 1.728
#define AVERAGE_CO_ENERGY 3.73
#define CO_POSITRON_FRACTION 0.0321

#define TAU_FE  29790.0        // 8.275 hours
#define TAU_MN  1266.0         // 21.1 minute
#define TAU_CR  77616.0        // 21.56 hours  
#define TAU_VN  1.38007e6      // 15.973 days
#define AVERAGE_FE_ENERGY 0.86
#define AVERAGE_MN_ENERGY 3.415
#define AVERAGE_CR_ENERGY 0.42
#define AVERAGE_VN_ENERGY 2.874

#define INTEGRATION_DELTA 0.001 // integration step of 0.0001 minutes

class radioactive
{
 
public:

	double rprocess_heating_rate(double t, double *gfrac);

	void decay_composition
	(std::vector<int> elem_Z, std::vector<int> elem_A, std::vector<double>& X, double t);

	double decay(std::vector<int> elem_Z, std::vector<int> elem_A, 
	       std::vector<double> X, double t, double *gfrac);

	double decay_energy_rate(int, int, double, double*);		
  //  double sample_particle_energy(int, int, double);
  //void   Fractions(int i, double *x, double t);

  //  static const double tau_56Ni =  760579.0;
  //static const double tau_56Co =  9.81599e06;   // decay time = 113.611 days 

  
};


#endif
