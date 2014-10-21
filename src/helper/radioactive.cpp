#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "radioactive.h"
#include "physical_constants.h"

namespace pc = physical_constants;


// energy (Mev) an probability (# per decay) of gamma rays from a Co56 decay
int n_co56_decays   = 17;
double co56_prob[]  = { 0.360,   1.000,   0.015,   0.137,  0.022,  0.670,
                         0.043,   0.158,   0.031,   0.079,  0.166,  0.058,
                         0.030,   0.074,   0.018,   0.009,  0.000};
double co56_energy[] = { 0.511,   0.847,   0.980,   1.040,  1.180,  1.240,
                         1.360,   1.770,   2.015,   2.030,  2.600,  3.010,
                         3.200,   3.250,   3.270,   3.450,  4.000};
 
int n_ni56_decays    =  7;
double ni56_prob[]   =  {1.000, 0.360, 0.360, 0.500, 0.870, 0.140, 0.000};
double ni56_energy[] =  {0.158, 0.270 ,0.480, 0.750, 0.812, 1.562, 2.000};



int    n_R_proc_fit = 11;
double R_proc_fit[] =  {17.608179, -2.0442059, -0.42565322, 0.39830095,  -0.0059089906,
			-0.054805836, 0.014068697, -0.00086706160, -5.7056758e-05, 2.6401842e-06,
			3.7186979e-07};


//--------------------------------------------------------------
// returns the energy from radioactive decay, per time 
// integrated over all species
//--------------------------------------------------------------
double radioactive::decay
(std::vector<int> elem_Z, std::vector<int> elem_A, std::vector<double> X, double t, double *gfrac)
{
  double total  = 0;
  double gtotal = 0;
  double gf;
  for (int k=0;k<elem_Z.size();k++)
  {
    int el_Z = elem_Z[k];
    int el_A = elem_A[k];
    double val = decay_energy_rate(el_Z,el_A,t,&gf);
    val = val*X[k]/(el_A*pc::m_p);;
    total  += val;
    gtotal += val*gf;
  }
  *gfrac = gtotal/total;
  return total;
       
       // adjust the compositions here
}
  
//--------------------------------------------------------------
// returns the energy from radioactive decay, per time per volume
//--------------------------------------------------------------
double radioactive::decay_energy_rate(int Z, int A, double t, double *gfrac)
{

  double total  = 0;
  double gtotal = 0;

  // Do 56Ni decay
  if ((Z == 28)&&(A == 56))
  {
    // exponential factors to be used
    double e_ni = exp(-t/TAU_NI);
    double e_co = exp(-t/TAU_CO);
    // number divided by decay time)
    double ni56 = (e_ni/TAU_NI);
    double co56 = 1.0/(TAU_NI-TAU_CO)*(e_ni - e_co);
    // get the energy from decays in ergs/s, using unit conversions
    double ni_E = ni56*(AVERAGE_NI_ENERGY*pc::Mev_to_ergs);
    double co_E = co56*(AVERAGE_CO_ENERGY*pc::Mev_to_ergs);

    gtotal = ni_E + 0.98*co_E;
    total  = ni_E + co_E;
  }
    
  // // Do 52Fe decay
  // if (ZONE::radio[i] == 26) 
  //   {
  //     // amount of 52Fe at t=0
  //     double fe52_0 = z.rho*z.R_gas[i]/(52*M_PROTON);

  //     // exponential factors to be used
  //     double e_fe = exp(-time/TAU_FE);
  //     double e_mn = exp(-time/TAU_MN);
  //     // number divided by decay time)
  //     double fe52 = fe52_0*(e_fe/TAU_FE);
  //     double mn52 = fe52_0/(TAU_FE-TAU_MN)*(e_fe - e_mn);
  //     // get the energy from decays in ergs/s, using unit conversions
  //     double fe_E = fe52*(AVERAGE_FE_ENERGY*MEV_TO_ERGS);
  //     double mn_E = mn52*(AVERAGE_MN_ENERGY*MEV_TO_ERGS);

  //     gtotal += fe_E + mn_E;
  //     total  += fe_E + mn_E;
  //   }

  //   // Do 48Cr decay
  //   if (ZONE::radio[i] == 24) 
  //   {
  //     // amount of 48cr at t=0
  //     double cr48_0 = z.rho*z.R_gas[i]/(48*M_PROTON);

  //     // exponential factors to be used
  //     double e_cr = exp(-time/TAU_CR);
  //     double e_vn = exp(-time/TAU_VN);
  //     // number divided by decay time)
  //     double cr48 = cr48_0*(e_cr/TAU_CR);
  //     double vn48 = cr48_0/(TAU_CR-TAU_VN)*(e_cr - e_vn);
  //     // get the energy from decays in ergs/s, using unit conversions
  //     double cr_E = cr48*(AVERAGE_CR_ENERGY*MEV_TO_ERGS);
  //     double vn_E = vn48*(AVERAGE_VN_ENERGY*MEV_TO_ERGS);

  //     gtotal += cr_E + vn_E;
  //     total  += cr_E + vn_E;
  //   }


  //   // R-process decay
  //   if (ZONE::radio[i] > 50)
  //   {
  //     double at = log10(time);
  //     double rproc = 0;
  //     for (int j=0;j<n_R_proc_fit;j++) rproc += R_proc_fit[j]*pow(at,j);
  //     rproc = pow(10.0,rproc)*z.rho;

  //     // fission fragments versus beta decay
  //     double fission_E  = 0.1*rproc;
  //     double beta_E     = 0.9*rproc;
      
  //     // account for lost neutrino energy
  //     beta_E = 0.75*beta_E;
      
  //     total  += (fission_E + beta_E);
  //     // half of the beta energy comes out as gamma-rays
  //     gtotal += 0.5*beta_E;
  //   }
      

  *gfrac = gtotal/total;
  return total;
}



// void RADIOACTIVE::Fractions(int i, double *x, double t)
// {
//   if (i == 28) 
//   {
//     x[0] = exp(-t/TAU_NI);
//     x[1] = TAU_CO/(TAU_NI-TAU_CO)*(exp(-t/TAU_NI) - exp(-t/TAU_CO));
//     x[2] = 1 - x[0] - x[1];
//   }

//   if (i == 26) 
//   {
//     x[0] = exp(-t/TAU_FE);
//     x[1] = TAU_MN/(TAU_FE-TAU_MN)*(exp(-t/TAU_FE) - exp(-t/TAU_MN));
//     x[2] = 1 - x[0] - x[1];
//   }

//   if (i == 24) 
//   {
//     x[0] = exp(-t/TAU_CR);
//     x[1] = TAU_VN/(TAU_CR-TAU_VN)*(exp(-t/TAU_CR) - exp(-t/TAU_VN));
//     x[2] = 1 - x[0] - x[1];
//   }
// }



//double radioactive::sample_particle_energy(int Z, int A, double t, gsl_rng *rangen)
//{

  // // Ni56
  // if ((Z == 28)&&(A == 56))
  // {
  
  //   // The Ratio Of Energy Coming Out In Nickel
  //   double E_Ni = exp(-t/TAU_NI);
  //   double E_Co = exp(-t/TAU_CO);
  //   double Ni_E = AVERAGE_NI_ENERGY*(E_Ni/TAU_NI);
  //   double Co_E = AVERAGE_CO_ENERGY/(TAU_NI - TAU_CO)*(E_Ni - E_Co);
  //   double nico_ratio = Ni_E/(Ni_E + Co_E);
  
  //   // pick emission wavelength
  //   double x_val;
  //   int d;
  //   double z1 = gsl_rng_uniform(rangen);
  //   if (z1 < nico_ratio)
  //   {
  //     while (true)
  //     {
  // 	double z2 = gsl_rng_uniform(rangen);
  // 	double z3 = gsl_rng_uniform(rangen);
  // 	d = (int)(n_ni56_decays*z2);
  // 	if (z3 < ni56_prob[d]) break;
  //     }
  //     x_val = ni56_energy[d];
  //   }
  //   else
  //   {
  //   while (true)
  //   {
  //     double z2 = gsl_rng_uniform(rangen);
  //     double z3 = gsl_rng_uniform(rangen);
  //     d = (int)(n_co56_decays*z2);
  //     if (z3 < co56_prob[d]) break;
  //   }
  //   x_val = co56_energy[d];
  //   }
  // }
  
  // return x_val;
  //}




