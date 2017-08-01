#ifndef _TRANSPORT_H
#define _TRANSPORT_H

#include <vector>
#include <list>
#include <mpi.h>

#include "particle.h"
#include "grid_general.h"
#include "cdf_array.h"
#include "locate_array.h"
#include "thread_RNG.h"
#include "spectrum_array.h"
#include "nlte_gas.h"
#include "ParameterReader.h"
#include "VoigtProfile.h"
#include "sedona.h"

using std::vector;

/* class RelativityData */
/* { */
/*   double gamma; */
/*   double beta; */
/*   double V[3]; */
/*   double D[3]; */
/*   double v_dot_D; */
/* }; */
  

class transport
{

 private:

  // arrays of particles
  std::list<particle> particles;
  int max_total_particles;

  // gas class for opacities
  nlte_gas gas;

  // pointer to parameter reader
  ParameterReader* params_;
  
  // MPI stuff
  int MPI_nprocs;
  int MPI_myID;
  MPI_Datatype MPI_real;
  int my_zone_start_, my_zone_stop_;
  double *src_MPI_block, *dst_MPI_block;
  double *src_MPI_zones, *dst_MPI_zones;

  // simulation parameters
  double step_size_;
  int    steady_state;
  int    radiative_eq;
  int    first_step_;
  int    verbose;
  int    use_detailed_lines_;
  int    omit_scattering_;
  int    store_Jnu_; 
  double maximum_opacity_;

  // current time in simulation
  double t_now_;

  // inner boundary
  double L_core_, r_core_, T_core_, time_core_, core_frequency_;
  // emission spectrum from inner boundary
  cdf_array<double> core_emission_spectrum_;
  // emission distribution across zones
  cdf_array<double> zone_emission_cdf_;

  // minimum and maximum temperatures
  double temp_max_value_, temp_min_value_;

  // class to hold output spectrum
  spectrum_array optical_spectrum;
  spectrum_array gamma_spectrum;

  // random number generator
  //mutable thread_RNG rangen;

   // random number generator
  gsl_rng *rangen;

  // Voigt profile class
  VoigtProfile voigt_profile_;

  // pointer to grid
  grid_general *grid;
  
  // the frequency grid for emissivity/opacity (Hz)
  locate_array nu_grid;

  // boundary conditions
  int boundary_in_reflect_, boundary_out_reflect_;

  // array to weight the emissivity(size of nu_grid)
  vector<real>  emissivity_weight_;

  // the zone opacity/emissivity variables
  vector< cdf_array<OpacityType> >    emissivity_;
  vector< vector<OpacityType> >       abs_opacity_;
  vector< vector<OpacityType> >       scat_opacity_;
  vector< vector<real> > J_nu_;
  vector<real> compton_opac;
  vector<real> photoion_opac;

  // the radiation quantities in the zone
  vector <real> e_rad;

  // -----------------------------------------
  // number of lines
  int n_lines_;
  // order list of the frequencies of each line
  vector<real> line_nu_;
  // the sqrt of ion mass for each line (needed for doppler width)
  vector<real> line_sqrt_Mion_;
  // line mean intensity
  vector< vector<real> > line_J_;
  // line extinction; multiply by 
  vector< vector<real> > line_opacity_;
  double line_velocity_width_;


  // setup functions
  void setup_core_emission();

  // opacity functions
  int   get_opacity(particle&, double, double&, double&);
  void   set_opacity();
  double klein_nishina(double);
  double blackbody_nu(double T, double nu);
  void   reduce_opacities();

  // creation of particles functions
  void   emit_particles(double dt);
  void   emit_inner_source(double dt);
  void   emit_radioactive(double dt);
  void   emit_thermal(double dt);
  void   emit_heating_source(double dt);

  void   create_isotropic_particle(int,PType,double,double);
  void   initialize_particles(int);
  void sample_photon_frequency(particle*);
  
  // special relativistic functions
  void   transform_comoving_to_lab(particle*);
  void   transform_lab_to_comoving(particle*);
  void   lorentz_transform(particle*, int);
  double dshift_comoving_to_lab(particle*);
  double dshift_lab_to_comoving(particle*);
  double do_dshift(particle*, int);

  //propagation of particles functions
  ParticleFate propagate(particle &p, double tstop);

  // scattering functions
  ParticleFate do_scatter(particle*, double);
  void compton_scatter(particle*);
  void isotropic_scatter(particle*, int);

  // radiation quantities functions
  void wipe_radiation();
  void reduce_radiation(double);
  void reduce_Tgas();

  // solve equilibrium temperature
  void solve_eq_temperature();
  double rad_eq_function(int,double);
  double temp_brent_method(int);


 public:
  

  //--------------------------------
  // constructor and defaults
  //------------------------------

  transport()
  {
    time_core_ = 0;
  }

  //----- functions ---------------
  
  // set things up
  void init(ParameterReader*, grid_general*);
  
  // run a transport step
  void step(double dt);

  // return 
  int n_particles() { return particles.size(); }

  // finalize and output spectra
  void output_spectrum();
  void output_spectrum(int);

  // print out functions
  void write_radiation_file(int, int);  

};


#endif
