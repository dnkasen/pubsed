#include <cstdlib>
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
  if (verbose) {
    cout << "# model file = " << model_file << "\n";
    cout << "# Model is a 1D_sphere\n"; }
  
  // type of system
  string system;
  infile >> system;

  // number of zones
  infile >> n_zones;
  z.resize(n_zones);
  r_out.resize(n_zones);
  vol.resize(n_zones);
  
  // read zone properties for a supernova remnant
  if (system == "SNR") 
    read_SNR_file(infile,verbose,1);
  else if (system == "standard") 
    read_SNR_file(infile,verbose,0);
  else {
    if (verbose) std::cout << " Don't recognize model type " << system << "; Exiting\n";
    exit(1); }

  infile.close();
}
    


void grid_1D_sphere::read_SNR_file(std::ifstream &infile, int verbose, int snr)
{
  // read header, general properties
  double texp;
  infile >> r_out.min;
  infile >> texp;
  this->t_now = texp;

  // set v at inner boundary = 0
  v_inner_ = 0;

  // read element isotopes, format is Z.A
  infile >> this->n_elems;
  for (int k=0;k<n_elems;k++) 
  {
    std::string species;
    infile >> species;
    int ind = species.find(".");
    std::string el_Z = species.substr(0,ind);
    std::string el_A = species.substr(ind+1,species.size() - ind);
    elems_Z.push_back(std::stoi(el_Z));
    elems_A.push_back(std::stoi(el_A));
  }

  // loop over zones and read
  for (int i=0; i<n_zones; i++)
  {
    // read state variables
    if (snr)
    {
      infile >> z[i].v[0];
      infile >> z[i].rho;
      infile >> z[i].T_gas;
      // assume homology for radius
      r_out[i] = z[i].v[0]*t_now;
    }
    else
    {
      infile >> r_out[i];
      infile >> z[i].v[0];
      infile >> z[i].rho;
      infile >> z[i].T_gas;
      std::cout << i << " " << r_out[i] << "\n";
    }
    // read composition
    z[i].mu = 0;
    for (int k=0;k<n_elems;k++)
    {
      double x;
      infile >> x;
      z[i].X_gas.push_back(x);
      z[i].mu += x*elems_A[k];
    }

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
    if (snr) cout << "#\n####### 1D SNR MODEL ##########\n";
    else cout << "#\n####### 1D STANDARD MODEL ##########\n";
    cout << "# n_x = " << n_zones << endl;
    cout << "# elems (n=" << n_elems << ") ";
    for (int k=0;k<n_elems;k++) cout << elems_Z[k] << "." << elems_A[k] << " ";
    cout << "\n#\n";

    // summed properties
    double tmass = 0;
    double ke    = 0;
    double re    = 0;
    std::vector<double>elem_mass(n_elems);
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
       
  }
}

  

//************************************************************
// expand the grid
//************************************************************
void grid_1D_sphere::expand(double e) 
{
  for (int i=0;i<n_zones;i++) r_out[i] *= e; 
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
  if(r >= r_out[r_out.size()-1] ) return -2;

  // find in zone array using stl algorithm up_bound and subtracting iterators
  int ind = r_out.locate(r);
  return ind;
}


//************************************************************
// Overly simple search to find zone
//************************************************************
int grid_1D_sphere::get_next_zone(const double *x, const double *D, int i, double r_core, double *l) const
{
  double rsq   = (x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
  double xdotD = (D[0]*x[0] + D[1]*x[1] + D[2]*x[2]);

  // distance to outer shell edge
  double r_o = r_out[i];
  double l_out = -1*xdotD + sqrt(xdotD*xdotD + r_o*r_o - rsq);


  double r_i = 0;
  int ind_in = i-1;
  if (i != 0)  r_i = r_out[i-1];
  if (r_core >= r_i) {
      r_i = r_core;
      ind_in = -1; }


  double l_in;
  if ((i == 0)&&(r_core == 0)) l_in = -1;
  else
  {
    double rad = xdotD*xdotD + r_i*r_i - rsq;
    if   (rad < 0)  l_in = -1;
    else l_in = -1*xdotD - sqrt(rad);
  } 

  // find shortest positive distance
  int ind;
  //double tiny = 1 + 1e-10;
  if ((l_out < l_in)||(l_in < 0))
  {
    ind = i + 1;
    if (ind == n_zones) ind = -2;
    *l = l_out; //*tiny;
  }
  else
  {
    ind = ind_in;
    *l = l_in; //*tiny;
  }

  return ind;
}



//************************************************************
// Write out the file
//************************************************************
void grid_1D_sphere::write_plotfile(int iw, double tt)
{
  // write ascii version, for convenience
  char zonefile[1000];
  sprintf(zonefile,"plt_%05d.dat",iw);

  std::ofstream outf;
  outf.open(zonefile);
  outf << std::setprecision(4);
  outf << std::scientific;

  outf << "# t = " << tt << "  ; rmin = " << r_out.min << "\n";
  outf << "#   r              rho              v              T_gas          T_rad         L_dep_nuc       L_emit_nuc\n";

  for (int i=0;i<n_zones;i++)
  {
    double rin = r_out.min;
    if (i > 0) rin = r_out[i-1];
    double T_rad = pow(z[i].e_rad/pc::a,0.25);

    outf << std::setprecision(8);
    outf << r_out[i] << "\t";
    outf << z[i].rho << "\t";
    outf << z[i].v[0] << "\t";
    outf << z[i].T_gas << "\t";
    outf << T_rad << "\t";
    outf << z[i].L_radio_dep << "\t";
    outf << z[i].L_radio_emit << "\t";
    outf << endl;
  }

  // write hdf5 file
  sprintf(zonefile,"plt_%05d.h5",iw);

  // open hdf5 file
  hid_t file_id = H5Fcreate( zonefile, H5F_ACC_TRUNC, H5P_DEFAULT,  H5P_DEFAULT);
 
  // print out r array
  hsize_t  dims_x[1]={(hsize_t)n_zones};
  float *xarr = new float[n_zones];
  for (int i=0;i<n_zones;i++) xarr[i] = r_out[i];
  H5LTmake_dataset(file_id,"r",1,dims_x,H5T_NATIVE_FLOAT,xarr);
  delete [] xarr;

  // print out r min
  // print out time
  hsize_t  dims_r[1]={1};
  float r0 = r_out.min;
  H5LTmake_dataset(file_id,"r_inner",1,dims_r,H5T_NATIVE_FLOAT,&r0); 

  hsize_t  dims_g[1]={(hsize_t) n_zones};
  write_hdf5_plotfile_zones(file_id, dims_g, 1, tt);

  write_integrated_quantities(iw,tt);

  H5Fclose (file_id);
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
// sample a random position within the spherical shell
//************************************************************
void grid_1D_sphere::sample_in_zone
(int i, std::vector<double> ran, double r[3])
{
  // inner radius of shell
  double r_0;
  if (i == 0) r_0 = r_out.min; 
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
void grid_1D_sphere::get_velocity(int i, double x[3], double D[3], double v[3], double *dvds)
{
  // radius in zone
  double rr = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);

  // linearly interpolate velocity here
  double v_0, r_0;
  if (i == 0) {v_0 = v_inner_; r_0 = r_out.min; }
  else {v_0 = z[i-1].v[0]; r_0 = r_out[i-1]; }
  double dr = rr - r_0;
  double dv_dr = (z[i].v[0] - v_0)/(r_out[i] - r_0);

  double vv = v_0 + dv_dr*dr;

  // assuming radial velocity 
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

  *dvds = dv_dr;  // not quite right, but upper limit

}

void grid_1D_sphere::get_radial_edges
(std::vector<double> &r, double &r0, std::vector<double> &v, double &v0) const
{
  for (int i=0;i<n_zones;i++)
  {
    r[i] = r_out[i];
    v[i] = z[i].v[0];
  }
  r0 = r_out.min;
  v0 = v_inner_;
}
void grid_1D_sphere::set_radial_edges
(const std::vector<double> r, const double r0, 
const std::vector<double> v, const double v0) 
{
  r_out.min = r0;
  v_inner_ = v0;
  for (int i=0;i<n_zones;i++)
  {
    r_out[i] = r[i];
    z[i].v[0] = v[i];

    // calculate shell volume
    double r0;
    if(i==0) r0 = r_out.min;
    else     r0 = r_out[i-1];
    vol[i] = 4.0*pc::pi/3.0*(r_out[i]*r_out[i]*r_out[i] - r0*r0*r0);
  }


}


