#include <gsl/gsl_rng.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <limits>

#include "hydro_1D_lagrangian.h"
#include "physical_constants.h"
namespace pc = physical_constants;


void hydro_1D_lagrangian::init(ParameterReader *params, grid_general *g)
{
  grid = g;
  time_ = grid->t_now;

  z_start_ = 0;
  nz_ = grid->n_zones;
  r_out_.resize(nz_);
  v_out_.resize(nz_);
  visq_.resize(nz_);
  eden_.resize(nz_);
  mass_.resize(nz_);

  gamfac_      = params->getScalar<double>("hydro_gamma_index");
  cfl_         = params->getScalar<double>("hydro_cfl");
  mean_particle_mass_ = params->getScalar<double>("hydro_mean_particle_mass");
  C_q_         = params->getScalar<double>("hydro_viscosity_parameter");
  M_center_    = params->getScalar<double>("hydro_central_point_mass");
  use_gravity_ = params->getScalar<int>("hydro_use_gravity");
  use_transport_ = params->getScalar<int>("hydro_use_transport");
  boundary_outflow_ = params->getScalar<int>("hydro_boundary_outflow");
  boundary_rigid_outer_wall_ = params->getScalar<int>("hydro_boundary_rigid_outer_wall");
  r_accrete_   = params->getScalar<double>("hydro_accrete_radius");

  // handle conflicting input
  if (boundary_outflow_ == 1 && boundary_rigid_outer_wall_ == 1)
  {
    printf("WARNING!!! Hydro outer boundary cannot be simultaneously outflow and outer wall! Setting to outflow\n");
    boundary_rigid_outer_wall_ = 1;
  }

  grid->get_radial_edges(r_out_,r_min_,v_out_,v_min_);
  v_min_ = params->getScalar<double>("hydro_v_piston");

  // thermal bomb parameters
  double E_bomb = params->getScalar<double>("hydro_bomb_energy");
  double r_bomb = params->getScalar<double>("hydro_bomb_radius");
  // default if r_bomb undefined, put  in inner 5 zones
  if (r_bomb == 0)
  {
    if (nz_ <= 10)
      r_bomb = r_out_[0]; // if you use index 1 note that this will not work if there is only 1 zone
    else
      r_bomb = r_out_[10];
  }
  double vol_bomb = 4.0*pc::pi/3.0*r_bomb*r_bomb*r_bomb;
  double p_bomb = (E_bomb/vol_bomb/3.0);

  // make sure things are calculated
  for (int i = 0;i < nz_;i++)
  {
    double vol = grid->zone_volume(i);
    // radiation pressure EOS
    //	grid->z[i].p_gas = pc::a*pow(grid->z[i].T_gas,4)/3.0;
    // gas pressure EOS
    grid->z[i].p_gas = pc::k*grid->z[i].rho/(mean_particle_mass_ * pc::m_p)*grid->z[i].T_gas;

    // add in bomb thermal energy
    if (r_out_[i] <= 3*r_bomb && E_bomb > 0.)
    {
      grid->z[i].p_gas += p_bomb*exp(-r_out_[i]*r_out_[i]/r_bomb/r_bomb);
      //assuming radiation pressure dominates and T_gas = T_r
      grid->z[i].T_gas = pow(3.0*grid->z[i].p_gas/pc::a,0.25);
    }
    grid->z[i].cs = sqrt(gamfac_ * grid->z[i].p_gas/grid->z[i].rho);
    eden_[i] = grid->z[i].p_gas/(gamfac_- 1.0)/grid->z[i].rho; // Units energy per mass

    grid->z[i].e_gas = eden_[i]; // this is what will be used for constructing eps_imc. Units energy per mass
    visq_[i] = compute_artificial_viscosity(i);
    mass_[i] = grid->z[i].rho*vol;
  }

  // time to write out
  time_write_ = 0;
  // clear output file
  std::ofstream output;
  output.open("hydro_data.dat");
  output.close();

}


double hydro_1D_lagrangian::get_time_step()
{
  double tstep = std::numeric_limits<double>::infinity();

  // courant condition
  int z0 = 0;
  for (int i = z_start_;i < nz_;i++)
  {
  	double dr = get_dr(i);
  	double v;
  	if (i == 0) v = fabs(v_min_) + grid->z[i].cs;
  	else   v = fabs(v_out_[i-1]) + grid->z[i].cs;
  	double dt = cfl_*dr/v;
    if (dt < tstep) z0 = i;
  	if (dt < tstep) tstep = dt;
  }
  //  std::cout << z0 << " " << tstep << "\n";
  return tstep;
}



void hydro_1D_lagrangian::step(double dt)
{

  // add mass to inner zone
  //grid->z[0].mass  += Mdot*dt;
  //grid->z[0].e_gas += Mdot*e_add*dt;
  // recalculate density/pressure

  // split inner zone if it's too big
 // if (r_out_[0] > r_core_*1.10)


  // -----------------------------------------------------
  // update the velocities
  // -----------------------------------------------------
  double msum = M_center_;
  for (int i=z_start_;i<nz_;i++)
  {
    // get zones to do derivatives
    int z1, z2;
    if (nz_ == 1) {z1 = i; z2 = i;}
    else
    {
      if (i < nz_-1) { z1 = i; z2 = i+1;} // note that this will not work if there is only one zone
      else {z1 = i-1; z2 = i;}
    }

    // gas pressure and viscosity gradiant (times -1)
    double dp = grid->z[z1].p_gas - grid->z[z2].p_gas;
    double dq  = visq_[z1] - visq_[z2];
    // outer outflow boundary
    if (i == nz_-1 && boundary_outflow_) {
      dq = visq_[z2];
      dp = grid->z[z2].p_gas; }

    // radial element
    double dr1 = get_dr(z1);
    double dr2 = get_dr(z2);
    double rhomean = 0.5*(grid->z[z2].rho*dr2 + grid->z[z1].rho*dr1);

    // acceleration from pressure/visc gradients
    double accel = (dp + dq)/rhomean;

    // gravitational acceleration
    if (use_gravity_)
    {
      msum  += mass_[i];
      accel += -1*msum*pc::G/r_out_[i]/r_out_[i];
    }

   // acceleration from radiation flux
   if (use_transport_)
     accel += grid->z[i].fr_rad/grid->z[i].rho;

   // override accleration of outer zone boundary if using rigid wall
   if (i == nz_-1 && boundary_rigid_outer_wall_)
       accel = 0.;

    // update velocities and boudnaries
    v_out_[i] += accel*dt;
    r_out_[i] += v_out_[i]*dt;

    // accrete (i.e. destroy) zones below accretion radius
    if (r_out_[i] <= r_accrete_)
    {
      std::cout << "accreted " << i << "\n";
      M_center_ += mass_[i];
      r_out_[i] = r_min_;
      z_start_ = i+1;
    }

    if (i == z_start_)
    {
      if ((v_out_[i] < 0)&&(time_ > 20))
      {
        M_center_ += mass_[i];
        z_start_ = i+1;
        std::cout << "accreted " << i << "\n";
      }
    }

  }
  r_min_ += v_min_*dt;

  // -----------------------------------------------------
  // get new densities, energies, pressures
  // -----------------------------------------------------
  for (int i=z_start_;i<nz_;i++)
  {
    // new volume, density, artificial viscosity
    double r1 = r_out_[i];
    double r0 = r_min_;
    if (i > 0) r0 = r_out_[i-1];

    double new_vol = 4.0*pc::pi/3.0*(r1*r1*r1 - r0*r0*r0);
    double new_rho = mass_[i]/new_vol;
    double new_q   = compute_artificial_viscosity(i);

    // Calculate adiabatic and viscosity total energy input/output
    double dtau = (1.0/new_rho - 1.0/grid->z[i].rho);
    double qgam = -0.5*(new_q + visq_[i])*dtau;

    // implicitly calculate energy update
    eden_[i] = (eden_[i] - 0.5*(grid->z[i].p_gas + new_q + visq_[i])*dtau)/
      (1.0 + 0.5*(gamfac_ - 1.0)*new_rho*dtau);

    // add in energy from radiation
    // not implicit yet
    if (use_transport_)
    {
	     eden_[i] += (grid->z[i].e_abs  - grid->z[i].L_thermal * grid->z[i].eps_imc)*dt/new_rho; // e_abs already has eps_imc factor included
    }

// zvec[i].enrg = (zvec_old[i].enrg - 0.5 *
  //      (zvec_old[i].press + zvec[i].visc + zvec_old[i].visc) * rhofac)
    //  / (1.0 + 0.5*(param.gamma - 1.0)*zvec[i].rho * rhofac);

    //if (qgam < 0) qgam = 0;
    //double Sgam = qgam - grid->z[i].p_gas*dtau;
    //grid->z[i].Sgam = (1 - grid->z[i].eps_imc)*Sgam*grid->z[i].mass;
    // determine total energy sources/sink
    // eps_imc should be by default forced to 1 if thermal_equilibrium is enforced
    //double dE = grid->z[i].eps_imc*Sgam;
    //if (use_transport)

    //double dE = qgam - grid->z[i].p_gas*dtau;
    //dE += (grid->z[i].e_abs - grid->z[i].L_thermal)*new_vol*dt/new_rho; // - grid->z[i].e_emit)/new_rho;
    //eden_[i] += dE;
    //std::cout << grid->z[i].L_thermal*new_vol*dt << " " << eden_[i]*new_rho << " " << dE << "\n";

	// if (eden_[i] < 0)
	// 	printf("ERROR: e_gas < 0; %d %e %e %e %e\n",i,qgam,grid->z[i].p_gas,eden_[i],dE);
  // if (eden_[i] < 0)
//      exit(1);
   //   eden_[i] = pc::a*pow(100,4.0);


    // calculate other thermo properties, given the energy and density
    // Remember, e_gas has units of ergs per gram, not ergs per cc
    double p_gas = (gamfac_ - 1)*eden_[i]*new_rho;
    grid->z[i].p_gas  = p_gas;
    grid->z[i].e_gas  = eden_[i];
    grid->z[i].cs     = sqrt(gamfac_*p_gas/new_rho);
    grid->z[i].rho    = new_rho;
    // radiation pressure EOS
    //grid->z[i].T_gas  = pow(3.0*p_gas/pc::a,0.25);
    // gas pressure EOS
    grid->z[i].T_gas = p_gas/(pc::k*grid->z[i].rho/mean_particle_mass_/pc::m_p);

    visq_[i] = new_q;

    if (isnan(grid->z[i].T_gas)) {
      std::cout << i << " T_gas is nan " << new_rho << " " << p_gas << "\n";
      exit(0);
    }
  }
  grid->set_radial_edges(r_out_,r_min_,v_out_,v_min_);
  time_ += dt;


  //**************
  // integrated sums
  if (time_ > time_write_*1.01)
  {
    std::ofstream output;
    output.open("hydro_data.dat",std::ios_base::app | std::ios_base::out);

    double E_grav = 0, E_rad = 0, E_ke = 0, M_enc = M_center_;
    for (int i=z_start_;i<nz_;i++)
    {
      double vol = grid->zone_volume(i);
      M_enc  += grid->z[i].rho*vol;
      //double vol = grid->zone_volume(i);
      E_grav += -1*pc::G*M_enc*mass_[i]/r_out_[i]*use_gravity_;
      E_rad  += eden_[i]*mass_[i];
      double vm;
      if (i == 0) vm = v_out_[0]/2.0;
      else vm = 0.5*(v_out_[i] + v_out_[i-1]);
      E_ke   += 0.5*mass_[i]*vm*vm;
    }
    output << time_ << " " <<  M_enc  - M_center_ << " " <<  M_center_ << " " << z_start_;
    output << " " << E_grav << " " << E_rad << " " << E_ke << " " << E_grav + E_rad + E_ke << "\n";
    time_write_ = time_;
    output.close();
  }
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
