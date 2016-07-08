#include <math.h>
#include <gsl/gsl_rng.h>
#include "transport.h"
#include "radioactive.h"
#include "physical_constants.h"

namespace pc = physical_constants;
using std::cout;

//------------------------------------------------------------
// emit new particles
//------------------------------------------------------------
void transport::emit_particles(double dt)
{
  emit_inner_source(dt);
  emit_radioactive(dt);
}

//------------------------------------------------------------
// sample photon frequency from local emissivity
//------------------------------------------------------------
void transport::sample_photon_frequency(particle *p)
{
  if (p->type == photon)
  {
    int ilam  = emissivity_[p->ind].sample(gsl_rng_uniform(rangen));
    p->nu = nu_grid.sample(ilam,gsl_rng_uniform(rangen));
    if (p->nu > 1e20) std::cout << "pnu " << p->nu << "\n";
  }
  else if (p->type == gammaray)
  {
    p->nu = 1;
  }
  else
  {
    p->nu = 1;
  }



}

//------------------------------------------------------------
// General function to create a particle in zone i
// emitted isotropically in the comoving frame. 
// Useful for thermal radiation emitted all througout
// the grid
//------------------------------------------------------------
void transport::create_isotropic_particle
(int i, PType type, double Ep, double t)
{
  particle p;

  // particle index
  p.ind = i;

  // particle type
  p.type = type;
  
  // random sample position in zone
  std::vector<double> rand;
  rand.push_back(gsl_rng_uniform(rangen));
  rand.push_back(gsl_rng_uniform(rangen));
  rand.push_back(gsl_rng_uniform(rangen));
  double r[3];
  grid->sample_in_zone(i,rand,r);
  p.x[0] = r[0];
  p.x[1] = r[1];
  p.x[2] = r[2];

  // emit isotropically in comoving frame
  double mu  = 1 - 2.0*gsl_rng_uniform(rangen);
  double phi = 2.0*pc::pi*gsl_rng_uniform(rangen);
  double smu = sqrt(1 - mu*mu);
  p.D[0] = smu*cos(phi);
  p.D[1] = smu*sin(phi);
  p.D[2] = mu;

  // sample frequency from local emissivity
  sample_photon_frequency(&p);
  
  // set packet energy
  p.e  = Ep;

  // lorentz transform from the comoving to lab frame
  transform_comoving_to_lab(&p);

  // set time to current
  p.t  = t;

  // add to particle vector
  particles.push_back(p);

}


//------------------------------------------------------------
// Initialize a constant number of particles
// in each zone
//------------------------------------------------------------
void transport::initialize_particles(int init_particles)
{
  int n_add = 0;
  for (int i=0;i<grid->n_zones;i++)
  {
    // lab frame energy
    double E_zone = grid->z[i].e_rad*grid->zone_volume(i);
    // particle energy
    double Ep = E_zone/(init_particles);

    // create init_particles particles
    if (E_zone > 0) {
      for (int q=0;q<init_particles;q++) 
	create_isotropic_particle(i,photon,Ep,t_now_);
      n_add += init_particles; }
  }
  if (verbose) cout << "# initializing with " << n_add << " particles\n";  
}


//------------------------------------------------------------
// Emit gamma-rays from radioactive decay
//------------------------------------------------------------
void transport::emit_radioactive(double dt)
{
  // number of radioctive particles to emit
  int total_n_emit = params_->getScalar<int>("particles_n_emit_radioactive");
  if (total_n_emit == 0) return;
  int n_emit = total_n_emit/(1.0*MPI_nprocs);
  
  radioactive radio;
  double gfrac;

  // calculate the total decay energy on the grid
  double L_tot = 0;
  double *gamma_frac = new double[grid->n_zones];
  for (int i=0;i<grid->n_zones;i++)
  {
    double vol  = grid->zone_volume(i);
    double L_decay = 
      radio.decay(grid->elems_Z,grid->elems_A,grid->z[i].X_gas,t_now_,&gfrac);
    L_decay = grid->z[i].rho*L_decay*vol;
    grid->z[i].L_radio_emit = L_decay;
    gamma_frac[i] = gfrac;
    L_tot += L_decay;
  }
  
  if (verbose) cout << "# total radioactive L = " << L_tot << " ergs/sec\n";
  if (L_tot == 0) return;
  double E_p = L_tot*dt/(1.0*n_emit);

  int n_add_tot = 0;
  // loop over zones for emission
  for (int i=0;i<grid->n_zones;i++)
  {
    double E_emit =  grid->z[i].L_radio_emit*dt;
    int n_add = floor(E_emit/E_p);
    
    // randomize last one
    double extra = (E_emit - E_p*n_add)/E_p;
    if (gsl_rng_uniform(rangen) < extra) n_add++;

    // check that we have enough space to add these particles
    if (particles.size()+n_add > this->max_total_particles) {
      if (verbose) cout << "# Out of particle space; not adding in\n";
      return; }
    
    // setup particles
    for (int q=0;q<n_add;q++) 
    {
      double t  = t_now_ + dt*gsl_rng_uniform(rangen);

      // determine if make gamma-ray or positron
      if (gsl_rng_uniform(rangen) < gamma_frac[i])
	create_isotropic_particle(i,gammaray,E_p,t);
      else
      {
	// positrons are just immediately made into photons
	grid->z[i].L_radio_dep += E_p;
	create_isotropic_particle(i,photon,E_p,t);
      }
    }

    n_add_tot += n_add;
  }
  
  if (verbose) cout << "# added " << n_add_tot << " radiaoctive particles per proc\n";
  delete[] gamma_frac;
}



//------------------------------------------------------------
// inject particles from a central luminous source
//------------------------------------------------------------
void transport::emit_inner_source(double dt)
{
  // get the emisison propoerties from lua file
  // this could be set to be a function if we want
  int total_n_emit    = params_->getScalar<int>("core_n_emit");
  if (total_n_emit == 0) return;
  int n_emit = total_n_emit/(1.0*MPI_nprocs);
  double Ep  = L_core_*dt/n_emit;
  
  if (particles.size() + n_emit > this->max_total_particles)
    {cout << "# Not enough particle space\n"; return; }

  // inject particles from the source
  for (int i=0;i<n_emit;i++)
  {
    particle p;
   
    if (r_core_ == 0)
    {
      // central emission
      p.x[0] = 0;
      p.x[1] = 0;
      p.x[2] = 0;
      // emit isotropically in comoving frame
      double mu  = 1 - 2.0*gsl_rng_uniform(rangen);
      double phi = 2.0*pc::pi*gsl_rng_uniform(rangen);
      double smu = sqrt(1 - mu*mu);
      p.D[0] = smu*cos(phi);
      p.D[1] = smu*sin(phi);
      p.D[2] = mu;
    } 
    else 
    {
      // pick initial position on photosphere
      double phi_core   = 2*pc::pi*gsl_rng_uniform(rangen);
      double cosp_core  = cos(phi_core);
      double sinp_core  = sin(phi_core);
      double cost_core  = 1 - 2.0*gsl_rng_uniform(rangen);
      double sint_core  = sqrt(1-cost_core*cost_core);
      // real spatial coordinates    
      double a_phot = r_core_ + r_core_*1e-10;
      p.x[0] = a_phot*sint_core*cosp_core;
      p.x[1] = a_phot*sint_core*sinp_core;
      p.x[2] = a_phot*cost_core;

      // pick photon propagation direction wtr to local normal                   
      double phi_loc = 2*pc::pi*gsl_rng_uniform(rangen);
      // choose sqrt(R) to get outward, cos(theta) emission         
      double cost_loc  = sqrt(gsl_rng_uniform(rangen));
      double sint_loc  = sqrt(1 - cost_loc*cost_loc);
      // local direction vector                     
      double D_xl = sint_loc*cos(phi_loc);
      double D_yl = sint_loc*sin(phi_loc);
      double D_zl = cost_loc;
      // apply rotation matrix to convert D vector into overall frame        
      p.D[0] = cost_core*cosp_core*D_xl-sinp_core*D_yl+sint_core*cosp_core*D_zl;
      p.D[1] = cost_core*sinp_core*D_xl+cosp_core*D_yl+sint_core*sinp_core*D_zl;
      p.D[2] = -sint_core*D_xl+cost_core*D_zl;
    }

    // set energy of packet
    p.e = Ep;

    // get emission frequency
    if (core_frequency_ > 0)
    {
      // constant single frequency emission
      p.nu = core_frequency_;
    }
    else
    {
      // sample frequency from blackbody 
      int ilam = core_emission_spectrum_.sample(gsl_rng_uniform(rangen));
      p.nu = nu_grid.sample(ilam,gsl_rng_uniform(rangen));

      // straight bin emission
      //int ilam = gsl_rng_uniform(rangen)*nu_grid.size(); 
      //p.e *= core_emis.get_value(ilam)*nu_grid.size(); 

    } 

    // get index of current zone
    p.ind = grid->get_zone(p.x);

    // lorentz transform from the comoving to lab frame
    transform_comoving_to_lab(&p);

    // set time to current
    p.t  = t_now_;
    
    // set type to photon
    p.type = photon;
  
    // add to particle vector
    particles.push_back(p);
  }

  if (verbose) 
    printf("# Core emitted %d particles (%d per proc)\n",total_n_emit,n_emit);
}


