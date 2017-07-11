#include <mpi.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>
#include <vector>
#include <cassert>
#include <list>
#include <algorithm>

#include "transport.h"
#include "ParameterReader.h"
#include "physical_constants.h"

using std::cout;
namespace pc = physical_constants;


//----------------------------------------------------------------------------
// Initialize the transport module
// Includes setting up the grid, particles,
// and MPI work distribution
//----------------------------------------------------------------------------
void transport::init(ParameterReader* par, grid_general *g)
{ 
  params_  = par;
  grid = g;

  // get mpi rank
  MPI_Comm_size( MPI_COMM_WORLD, &MPI_nprocs );
  MPI_Comm_rank( MPI_COMM_WORLD, &MPI_myID  );
  MPI_real = ( sizeof(real)==4 ? MPI_FLOAT : MPI_DOUBLE );
  verbose = (MPI_myID==0);

  // determine my zones to work on
  int nz = grid->n_zones;
  int blocks = floor(nz/MPI_nprocs);
  int remainder = nz - blocks*MPI_nprocs;

  int rcount = 0;
  for (int i=0;i<MPI_nprocs;i++)
  {
    int start = i*blocks + rcount;
    int stop  = start + blocks;
    if (rcount < remainder) { stop += 1; rcount += 1;}
    if (i == MPI_myID) 
    {
      my_zone_start_ = start;
      my_zone_stop_  = stop;
    }
  }
  // arrays for communication 
  src_MPI_block = new double[Max_MPI_Blocksize];
  dst_MPI_block = new double[Max_MPI_Blocksize];
  src_MPI_zones = new double[nz];
  dst_MPI_zones = new double[nz];

//  std::cout << MPI_myID <<  " " << my_zone_start_ << " " << my_zone_stop_ << 
//   " " << my_zone_stop_ - my_zone_start_ << "\n";

  // setup and seed random number generator
  const gsl_rng_type * TypeR;
  gsl_rng_env_setup();
  gsl_rng_default_seed = (unsigned int)time(NULL) + MPI_myID;
  TypeR = gsl_rng_default;
  rangen = gsl_rng_alloc (TypeR);
  
  // read relevant parameters
  max_total_particles = params_->getScalar<int>("particles_max_total");
  radiative_eq    = params_->getScalar<int>("transport_radiative_equilibrium");
  steady_state    = (params_->getScalar<int>("transport_steady_iterate") > 0);
  temp_max_value_ = params_->getScalar<double>("limits_temp_max");
  temp_min_value_ = params_->getScalar<double>("limits_temp_min");

  // initialize the frequency grid  
  std::vector<double> nu_dims = params_->getVector<double>("transport_nu_grid");
  if ((nu_dims.size() != 4)&&(nu_dims.size() != 3)) {
    cout << "# improperly defined nu_grid; need {nu_1, nu_2, dnu, (log?)}; exiting\n";
    exit(1); }
  if (nu_dims.size() == 3)
    nu_grid.init(nu_dims[0],nu_dims[1],nu_dims[2]);
  if (nu_dims.size() == 4)
  {
    if (nu_dims[3] == 1) nu_grid.log_init(nu_dims[0],nu_dims[1],nu_dims[2]);
    else nu_grid.init(nu_dims[0],nu_dims[1],nu_dims[2]);
  }
  if (verbose)
  {
     std::cout << "# frequency grid: n = " << nu_grid.size() << "\n";
  }

  // intialize output spectrum
  std::vector<double>stg = params_->getVector<double>("spectrum_time_grid");
  std::vector<double>sng = params_->getVector<double>("spectrum_nu_grid");
  int nmu  = params_->getScalar<int>("spectrum_n_mu");
  int nphi = params_->getScalar<int>("spectrum_n_phi");
  optical_spectrum.init(stg,sng,nmu,nphi);
  std::vector<double>gng = params_->getVector<double>("gamma_nu_grid");  
  gamma_spectrum.init(stg,sng,nmu,nphi);
  
  // initialize nlte_gas class
  std::string atomdata = params_->getScalar<string>("data_atomic_file");  
  int use_nlte = params_->getScalar<int>("opacity_use_nlte");
  gas.initialize(atomdata,grid->elems_Z,grid->elems_A,nu_grid,use_nlte);

  std::string fuzzfile = params_->getScalar<string>("data_fuzzline_file");  
  int nl = gas.read_fuzzfile(fuzzfile);
  if (verbose) std::cout << "# From fuzzfile \"" << fuzzfile << "\" " << 
		 nl << " lines used\n";
  
  // set gas opacity flags and parameters
  gas.epsilon_      = params_->getScalar<double>("opacity_epsilon");
  gas.grey_opacity_ = params_->getScalar<double>("opacity_grey_opacity");
  gas.use_electron_scattering_opacity 
    = params_->getScalar<int>("opacity_electron_scattering");
  gas.use_line_expansion_opacity  
    = params_->getScalar<int>("opacity_line_expansion");
  gas.use_fuzz_expansion_opacity  
    = params_->getScalar<int>("opacity_fuzz_expansion");
  gas.use_bound_free_opacity  
    = params_->getScalar<int>("opacity_bound_free");
  gas.use_bound_bound_opacity  
    = params_->getScalar<int>("opacity_bound_bound");
  gas.use_free_free_opacity  
    = params_->getScalar<int>("opacity_free_free");
  double min_ext = params_->getScalar<double>("opacity_minimum_extinction");
  //maximum_opacity_ = params_->getScalar<double>("opacity_maximum_opacity");

  gas.set_minimum_extinction(min_ext);
  first_step_ = 1;

  if (verbose) gas.print_properties();

  // set up outer inner boundary condition
  boundary_in_reflect_ = params_->getScalar<int>("transport_boundary_in_reflect");
  boundary_out_reflect_ = params_->getScalar<int>("transport_boundary_out_reflect");

  // parameters for treatment of detailed lines
  use_detailed_lines_  = params_->getScalar<double>("opacity_lines");
  line_velocity_width_ = params_->getScalar<double>("line_velocity_width");
  gas.line_velocity_width_ = line_velocity_width_;

  // get line frequencies and ion masses
  line_nu_ = gas.get_line_frequency_list();
  line_sqrt_Mion_ = gas.get_line_ion_mass_list();
  for (size_t i=0;i<line_sqrt_Mion_.size();i++)
    line_sqrt_Mion_[i] = sqrt(line_sqrt_Mion_[i]);

  // parameters for storing opacities
  n_lines_ = gas.get_number_of_lines();
  omit_scattering_ = params_->getScalar<int>("opacity_no_scattering");
  store_Jnu_ = params_->getScalar<int>("transport_store_Jnu");
  // sanity check
  if ((!store_Jnu_)&&(use_nlte)) 
    std::cout << "WARNING: not storing Jnu while using NLTE; Bad idea!\n";

  line_opacity_.resize(grid->n_zones);
  abs_opacity_.resize(grid->n_zones);
  if (!omit_scattering_) scat_opacity_.resize(grid->n_zones);
  emissivity_.resize(grid->n_zones);
  J_nu_.resize(grid->n_zones);

  for (int i=0; i<grid->n_zones;  i++)
  {
    if (use_detailed_lines_)
      line_opacity_[i].resize(n_lines_);
   
    // allocate absorptive opacity
    try {
      abs_opacity_[i].resize(nu_grid.size()); }
    catch (std::bad_alloc const&) {
      cout << "Memory allocation fail!" << std::endl; }
   
    // allocate scattering opacity
    if (!omit_scattering_)
    {
      try {
        scat_opacity_[i].resize(nu_grid.size()); }
     catch (std::bad_alloc const&) {
        cout << "Memory allocation fail!" << std::endl; }
    }

    // allocate emissivity
    emissivity_[i].resize(nu_grid.size());

    // allocate Jnu (radiation field)
    if (store_Jnu_)
      J_nu_[i].resize(nu_grid.size());
    else 
      J_nu_[i].resize(1);
  }
  compton_opac.resize(grid->n_zones);
  photoion_opac.resize(grid->n_zones);

  // setup emissivity weight  -- debug
  emissivity_weight_.resize(nu_grid.size());
  double norm = 0;
  for (int j=0;j<nu_grid.size();j++)
  { 
    //double nu  = nu_grid.center(j);
    double w = 1.0; // - 1.0/(1.0 + pow(nu/3e15,2.0));
    //if (nu > 2e15) w = 100;
    emissivity_weight_[j] = w;
    norm += w;
  }
  for (int j=0;j<nu_grid.size();j++) emissivity_weight_[j] *= nu_grid.size()/norm;
 // for (int j=0;j<nu_grid.size();j++) std::cout << nu_grid.center(j) << " " << emissivity_weight_[j] << "\n";

  // allocate space for emission distribution function across zones
  zone_emission_cdf_.resize(grid->n_zones);

  // read pamaeters for core emission and setup
  setup_core_emission();

  // initialize time
  t_now_ = g->t_now;

  // initialize particles
  int n_parts = params_->getScalar<int>("particles_n_initialize");
  initialize_particles(n_parts);
}


void transport::setup_core_emission()
{

  // -----------------------------------------------------------
  // set up inner boundary emission properties
  r_core_         = params_->getScalar<double>("core_radius");
  T_core_         = params_->getScalar<double>("core_temperature");
  core_frequency_ = params_->getScalar<double>("core_photon_frequency");
  L_core_         = params_->getFunction("core_luminosity", 0);
  time_core_      = params_->getScalar<double>("core_timescale");

  // set blackbody from L and R if appropriate
  if ((L_core_ !=0)&&(r_core_ != 0)&&(T_core_ == 0))
    T_core_ = pow(L_core_/(4.0*pc::pi*r_core_*r_core_*pc::sb),0.25);


  int total_n_emit    = params_->getScalar<int>("core_n_emit");
  if (total_n_emit > 0)
  {
    // allocate and set up core emission spectrum
    core_emission_spectrum_.resize(nu_grid.size());

    std::string core_spectrum_filename = params_->getScalar<string>("core_spectrum_file");  
    vector<double> cspec_nu, cspec_Lnu;
    if (core_spectrum_filename != "") 
    {
     double x1,x2;
      std::ifstream specfile;
      specfile.open(core_spectrum_filename.c_str());
      if (!specfile.is_open()) 
      {
        if (verbose) std::cout << "Can't open core_spectrum_file " << core_spectrum_filename << "\n";
        core_spectrum_filename = "";
      }
      else while (!specfile.eof( ))   
      {
        specfile >> x1;
        specfile >> x2;
        cspec_nu.push_back(x1);
        cspec_Lnu.push_back(x2);
     }
   }

    // set up emission spectrum 
   double L_sum = 0;
   for (int j=0;j<nu_grid.size();j++)
   { 
      double nu  = nu_grid.center(j);
      double dnu = nu_grid.delta(j);

      // read in spectrum
      if (core_spectrum_filename != "") 
      {
        double Lnu;
        int ind = lower_bound(cspec_nu.begin(),cspec_nu.end(),nu)- cspec_nu.begin() - 1;
        if (ind < 0) Lnu = 0;
        else if ((size_t)ind >= cspec_nu.size()-1) Lnu = 0;
       else 
       {
          //double slope = 0; //(cspec_Lnu[ind+1] - cspec_Lnu[ind-1])/(cspec_nu[ind+1] - cspec_nu[ind]);
          Lnu = cspec_Lnu[ind]; // + slope*(nu - cspec_nu[ind]);
        }
        core_emission_spectrum_.set_value(j,Lnu*dnu*emissivity_weight_[j]); 
       L_sum += Lnu*dnu;
      }
      else
      // blackbody spectrum 
     {
        double bb;
        if (T_core_ <= 0) bb = 1;
        else bb = blackbody_nu(T_core_,nu);
        core_emission_spectrum_.set_value(j,bb*dnu*emissivity_weight_[j]);
        // blackbody flux is pi*B(T)
       L_sum += 4.0*pc::pi*r_core_*r_core_*pc::pi*bb*dnu;
      }
    }
   core_emission_spectrum_.normalize(); 
   if (L_core_ == 0) L_core_ = L_sum;

    if (verbose) 
    {
      if (core_spectrum_filename != "") 
        cout << "# Inner source luminosity (at t = 0) = " << L_core_ << 
       " erg/s, read from file " << core_spectrum_filename << "\n";
      else
        cout << "# Inner source luminosity = " << L_core_ << 
        "  erg/s, from a blackbody T = " << T_core_ << "\n";
    }
  } 
  // -----------------------------------------------------------

}

