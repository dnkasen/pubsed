#include <gsl/gsl_rng.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <iostream>

#include "hydro_1D_lagrangian.h"
#include "physical_constants.h"
namespace pc = physical_constants;

void hydro_1D_lagrangian::init(ParameterReader *params, grid_general *g)
{
  grid = g;
  t_now = g->t_now;

  nz_ = grid->n_zones;
  r_out_.resize(nz_);
  v_out_.resize(nz_);
  visq_.resize(nz_);
  eden_.resize(nz_);
  mass_.resize(nz_);

  gamfac_ = params->getScalar<double>("hydro_gamma_index");
  cfl_ = params->getScalar<double>("hydro_cfl");
  C_q_ = params->getScalar<double>("hydro_viscosity_parameter");

  grid->get_radial_edges(r_out_,r_min_,v_out_,v_min_);

  // make sure things are calculated
  for (int i = 0;i < nz_;i++)
  {
  	grid->z[i].p_gas = pc::k*grid->z[i].rho/pc::m_p*grid->z[i].T_gas;
  	grid->z[i].cs = sqrt(gamfac_ * grid->z[i].p_gas/grid->z[i].rho);
 	eden_[i] = grid->z[i].p_gas/(gamfac_- 1.0)/g->z[i].rho;
    visq_[i] = compute_artificial_viscosity(i);
    mass_[i] = grid->z[i].rho*grid->zone_volume(i);
  }

}


double hydro_1D_lagrangian::get_time_step()
{	
  double tstep = std::numeric_limits<double>::infinity();

  // courant condition
  for (int i = 0;i < nz_;i++)
  {
  	double dr = get_dr(i);
  	double v;
  	if (i == 0) v = fabs(v_min_) + grid->z[i].cs;
  	else   v = fabs(v_out_[i-1]) + grid->z[i].cs;
  	double dt = cfl_*dr/v;
  	if (dt < tstep) tstep = dt;
  }
  return tstep;
}



void hydro_1D_lagrangian::step(double dt)
{


 //  // -----------------------------------------------------
 //  // compute gravity etc..
 //  // -----------------------------------------------------
 //  double msum = M_gravity;
 // for (int i=0;i<n_zones;i++)
 // {
    // compute gravity
 // msum += grid->z[i].mass;
 //    if (use_gravity)
 //      grid->z[i].grav_up = -1*msum*NEWTON_G/grid->z[i].r_up/grid->z[i].r_up;
 //    else grid->z[i].grav_up = 0;
 
 //    // compute p stiff
 //    //double dstop = 5e-5;
 //    // double dsize = (dr/grid->z[i].r_dn);
 //    //double w  = 1.0/(1 + pow(dsize/dstop,20));
 //    //grid->z[i].p_stiff = w*RAD_CONST*pow(1e7,4);
 //  }
  

  // -----------------------------------------------------
  // update the velocities  
  // -----------------------------------------------------
  for (int i=0;i<nz_;i++)
  { 
    // get zones to do derivatives
    int z1, z2;
    if (i < nz_-1) { z1 = i; z2 = i+1;}
    else {z1 = i-1; z2 = i;}

    // gas pressure and viscosity gradiant (times -1)
    double dp = grid->z[z1].p_gas - grid->z[z2].p_gas;
    double dq  = visq_[z1] - visq_[z2];

    // radial element
    double dr1 = get_dr(z1);
    double dr2 = get_dr(z2);
    double rhomean = 0.5*(grid->z[z2].rho*dr2 + grid->z[z1].rho*dr1);
    
    // acceleration from pressure/visc gradients
    double accel = (dp + dq)/rhomean; 
    
 //    // acceleration from monte carlo flux
 //    if (use_transport) 
 //    {
 //      grid->z[i].v_up += grid->z[i].fx_rad/grid->z[i].rho*dt;
 //      // zone[z].f_rad = zone[z].f_rad/zone[z].rho/C_LIGHT; 
 //      // radiation pressure gradiant in diffusion regime
 //      //if (use_transport) //&&(zone[z].tau > tau_diffuse)) 
 //      // dp += (zone[z1].E_dif/3.0 - zone[z2].E_dif/3.0);
 //    }
  
 //    // acceleration from gravity
 //    if (use_gravity)
 //      grid->z[i].v_up += grid->z[i].grav_up *dt;

    // update velocities and boudnaries
    v_out_[i] += accel*dt;
    r_out_[i] += v_out_[i]*dt;
  }

 //  // inner boundary conditions
 //  grid->z[0].v_dn = v_inner;
 //  //if (v_inner == 0) zone[0].v_dn = zone[0].v_up;

 //  // boundaries
 //  //if (v_inner == 0) grid->z[0].r_dn  = 0; 
 //  grid->z[0].r_dn += grid->z[0].v_dn*dt;
 
 //  if (grid->z[0].r_dn < 0) grid->z[0].r_dn = 0;
 //  grid->z[n_zones-1].r_up += grid->z[n_zones-1].v_up*dt;

  // -----------------------------------------------------
  // get new densities, energies, pressures
  // -----------------------------------------------------
  for (int i=0;i<nz_;i++)
  {
    // new volume, density, artificial viscosity
    double r1 = r_out_[i];
    double r0 = r_min_;
    if (i != 0) r0 = r_out_[i-1];
    double new_vol = 4.0*pc::pi/3.0*(r1*r1*r1 - r0*r0*r0);
    double new_rho = mass_[i]/new_vol;
    double new_q   = compute_artificial_viscosity(i);

    // Calculate adiabatic and viscosity total energy input/output
    double dtau = (1.0/new_rho - 1.0/grid->z[i].rho);
    double qgam = -0.5*(new_q + visq_[i])*dtau;
 
    //if (qgam < 0) qgam = 0;
    //double Sgam = qgam - grid->z[i].p_gas*dtau;
    //grid->z[i].Sgam = (1 - grid->z[i].eps_imc)*Sgam*grid->z[i].mass;
    // determine total energy sources/sink
    // eps_imc should be by default forced to 1 if thermal_equilibrium is enforced
    //double dE = grid->z[i].eps_imc*Sgam;
    //if (use_transport) 
    //  dE += (grid->z[i].e_abs*dt - grid->z[i].e_emit)/new_rho;
    
    double dE = qgam - grid->z[i].p_gas*dtau;
    eden_[i] += dE;

	if (eden_[i] < 0)   
	 	printf("ERROR: e_gas < 0; %d %e %e %e %e\n",i,qgam,grid->z[i].p_gas,eden_[i],dE);
 
    // calculate other thermo properties, given the energy and density
    // Remember, in Sedona e_gas has units of ergs per gram, not ergs per cc
    double p_gas = (gamfac_ - 1)*eden_[i]*new_rho;
    grid->z[i].p_gas  = p_gas;  
    grid->z[i].cs     = sqrt(gamfac_*p_gas/new_rho);
    grid->z[i].T_gas  = p_gas/(pc::k*new_rho/grid->z[i].mu/pc::m_p);
    grid->z[i].rho    = new_rho;
    visq_[i] = new_q;
  }
  grid->set_radial_edges(r_out_,r_min_,v_out_,v_min_);

}


double hydro_1D_lagrangian::compute_artificial_viscosity(int i)
{
  double dr, dv;
  if (i == 0)
  { 
    dr = r_out_[0] - r_min_;
    dv = v_min_ - v_out_[0];
  }
  else
  {
    dr = r_out_[i] - r_out_[i-1];
    dv = v_out_[i-1] - v_out_[i];
  }

  // limit fluctuations
  //double dvdx = -1*dv/dr;
  //if (dvdx > -0.02) dv = 0;
  
  if (dv < 0) dv = 0;
  double q = C_q_*grid->z[i].rho*dv*dv;
    
  // linear term
  //  q += L_q*grid->z[i].rho*grid->z[i].cs*dv;
  return q;
}


// double HYDRO::Get_Total_Energy()
// {

//   double total_energy = 0;

//   int n_zones = grid->n_x;

//   double msum = M_gravity;

//   double grav_potential = 0.;
//   double kinetic = 0.;
//   double internal = 0.;

//   for (int i=0;i<n_zones;i++)
//   {

//     msum += grid->z[i].mass;
//     grav_potential = -1*msum*NEWTON_G/grid->z[i].r_up;

//     //    kinetic = 0.5 * (grid->z[i].mass + grid->z[i+1].mass) * 0.5 * pow(grid->z[i].v_up,2); // castor eq 3.7
//     kinetic = 0.5 * (grid->z[i].mass) * pow(grid->z[i].v[0],2); 

//     internal = grid->z[i].mass * grid->z[i].e_gas;

//     total_energy += grav_potential + kinetic + internal;

//   }

//   return total_energy;

// }





