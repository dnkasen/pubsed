#include <mpi.h>
#include <fstream>
#include <iostream>
#include <iomanip>   
#include <math.h>
#include <cassert>

#include "grid_1D_sphere.h"
#include "physical_constants.h"

namespace pc = physical_constants;

using std::string;
using std::cout;
using std::endl;

//------------------------------------------------------------
// initialize the zone geometry from model file
//------------------------------------------------------------
void grid_1D_sphere::read_model_file(ParameterReader* params)
{
  // verbocity
  int my_rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
  const int verbose = (my_rank == 0);

  // open up the model file, complaining if it fails to open
  string model_file = params->getScalar<string>("model_file");
  std::ifstream infile;
  infile.open(model_file.c_str());
  if(infile.fail())
  {
    if (verbose) cout << "Err: can't read model file: " << model_file << endl;
    exit(4);
  }

  // geometry of model
  infile >> grid_type;
  if(grid_type != "1D_sphere") 
  {
    if (verbose) cout << "Err: grid_type param disagrees with the model file\n";
    exit(4);
  }
  if (verbose) cout << "# Model is a 1D_sphere\n";
  
  // type of system
  string system;
  infile >> system;

  // number of zones
  infile >> n_zones;
  z.resize(n_zones);
  r_out.resize(n_zones);
  vol.resize(n_zones);

  // read zone properties for a supernova remnant
  if(system == "SNR") read_SNR_file(infile,verbose);

  infile.close();
}
    

void grid_1D_sphere::read_SNR_file(std::ifstream &infile, int verbose)
{
  // general properties
  double texp;
  infile >> r_out.min;
  infile >> texp;

  this->t_now = texp;

  // read element isotopes
  infile >> this->n_elems;
  for (int k=0;k<n_elems;k++) 
  {
    double species;
    infile >> species;
    int el_Z  = floor(species);
    int el_A = round(100*(species - el_Z));
    elems_Z.push_back(el_Z);
    elems_A.push_back(el_A);
  }

  // loop over zones and read
  for (int i=0; i<n_zones; i++)
  {
    // read state variables
    infile >> z[i].v[0];
    infile >> z[i].rho;
    infile >> z[i].T_gas;

    // read composition
    for (int k=0;k<n_elems;k++)
    {
      double x;
      infile >> x;
      z[i].X_gas.push_back(x);
    }

    // assume homology for radius
    r_out[i] = z[i].v[0]*texp;
      
    // assume LTE radiation field to start
    z[i].e_rad = pc::a*pow(z[i].T_gas,4);
      
    // calculate shell volume
    double r0;
    if(i==0) r0 = r_out.min;
    else     r0 = r_out[i-1];
    vol[i] = 4.0*pc::pi/3.0*(r_out[i]*r_out[i]*r_out[i] - r0*r0*r0);
  }
  
  // print out properties of the model
  if (verbose) 
  {
    cout << "#\n####### SNR MODEL ##########\n";
    cout << "# n_x = " << n_zones << endl;
    cout << "# elems (n=" << n_elems << ") ";
    for (int k=0;k<n_elems;k++) cout << elems_Z[k] << "." << elems_A[k] << " ";
    cout << "\n#\n";

    // summed properties
    double tmass = 0;
    double ke    = 0;
    double re    = 0;
    double *elem_mass = new double[n_elems];
    for (int k=0;k<n_elems;k++) elem_mass[k] = 0;

    // calculate some useful summed properties
    for (int i=0;i<n_zones;i++)
    {
      tmass += z[i].rho*vol[i];
      for (int k=0;k<n_elems;k++) elem_mass[k] += vol[i]*z[i].rho*z[i].X_gas[k];
      ke += 0.5*z[i].rho*vol[i]*z[i].v[0]*z[i].v[0];
      re += z[i].e_rad*vol[i];
    }

    printf("# mass = %.4e (%.4e Msun)\n",tmass,tmass/pc::m_sun);
    for (int k=0;k<n_elems;k++) {
      cout << "# " << elems_Z[k] << "." << elems_A[k] <<  "\t";
      cout << elem_mass[k] << " (" << elem_mass[k]/pc::m_sun << " Msun)\n"; }
    printf("# kinetic energy   = %.4e\n",ke);
    printf("# radiation energy = %.4e\n",re);
    cout << "##############################\n#\n";
    
    delete elem_mass;
   
  }
}

  

//************************************************************
// expand the grid
//************************************************************
void grid_1D_sphere::expand(double e) 
{
    for (int i=0;i<n_zones;i++) r_out[i] *= e; 
    r_inner   *=e;
    r_out.min *=e;

    // recalculate shell volume
    for (int i=0;i<n_zones;i++) 
    {
      double r0;
      if(i==0) r0 = r_out.min;
      else     r0 = r_out[i-1];
      vol[i] = 4.0*pc::pi/3.0*(r_out[i]*r_out[i]*r_out[i] - r0*r0*r0);
    }
  
}

//************************************************************
// Overly simple search to find zone
//************************************************************
int grid_1D_sphere::get_zone(const double *x) const
{
  double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);

  // check if off the boundaries
  if(r < r_out.min             ) return -1;
  if(r > r_out[r_out.size()-1] ) return -2;

  // find in zone array using stl algorithm up_bound and subtracting iterators
  return r_out.locate(r);
}




//************************************************************
// Write out the file
//************************************************************
void grid_1D_sphere::write_out(int iw)
{
  char zonefile[1000];
  char base[1000];

  if (iw < 10) sprintf(base,"_0000%d",iw);
  else if (iw < 100) sprintf(base,"_000%d",iw);
  else if (iw < 1000) sprintf(base,"_00%d",iw);
  else if (iw < 10000) sprintf(base,"_0%d",iw);
  else sprintf(base,"_%d",iw);
  sprintf(zonefile,"ray%s",base);

  std::ofstream outf;
  outf.open(zonefile);
  outf << std::setprecision(4);
  outf << std::scientific;

  for (int i=0;i<n_zones;i++)
  {
    double rin = r_inner;
    if (i > 0) rin = r_out[i-1];
    double T_rad = pow(z[i].e_rad/pc::a,0.25);
    double rc = 0.5*(r_out[i] + rin);

    outf << rc << "\t";
    outf << z[i].rho << "\t";
    outf << z[i].v[0] << "\t";
    outf << z[i].T_gas << "\t";
    outf << T_rad << "\t";
    outf << endl;
  }


}

//************************************************************
// return volume of zone (precomputed)
//************************************************************
double  grid_1D_sphere::zone_volume(const int i) const
{
  assert(i >= 0);
  return vol[i];
}


//************************************************************
// return length of zone
//************************************************************
double  grid_1D_sphere::zone_min_length(const int i) const
{
  assert(i >= 0);
  if (i == 0) return (r_out[i] - r_out.min);
  else return (r_out[i] - r_out[i-1]);
}





//************************************************************
// sample a random position within the spherical shell
//************************************************************
void grid_1D_sphere::sample_in_zone
(int i, std::vector<double> ran, double r[3])
{
  // inner radius of shell
  double r_0;
  if (i == 0) r_0 = r_inner; 
  else r_0 = r_out[i-1];

  // thickness of shell
  double dr = r_out[i] - r_0;
  // sample radial position in shell
  r_0 = r_0 + dr*ran[0];

  // random spatial angles
  double mu  = 1 - 2.0*ran[1];
  double phi = 2.0*pc::pi*ran[2];
  double sin_theta = sqrt(1 - mu*mu);

  // set the real 3-d coordinates
  r[0] = r_0*sin_theta*cos(phi);
  r[1] = r_0*sin_theta*sin(phi);
  r[2] = r_0*mu;
}



//************************************************************
// get the velocity vector 
//************************************************************
void grid_1D_sphere::velocity_vector(int i, double x[3], double v[3])
{
  // radius in zone
  double rr = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);

  // linearly interpolate velocity here
  double v_0, r_0;
  if (i == 0) {v_0 = 0; r_0 = r_inner; }
  else {v_0 = z[i-1].v[0]; r_0 = r_out[i-1]; }
  double dr = rr - r_0;
  double dv_dr = (z[i].v[0] - v_0)/(r_out[i] - r_0);

  double vv = v_0 + dv_dr*dr;

  // assuming radial velocity (may want to interpolate here)
  v[0] = x[0]/rr*vv;
  v[1] = x[1]/rr*vv;
  v[2] = x[2]/rr*vv;

  // check for pathological case
  if (rr == 0)
  {
    v[0] = 0;
    v[1] = 0;
    v[2] = 0;
  }

}

