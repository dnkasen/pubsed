#include <cstdlib>
#include <mpi.h>
#include <fstream>
#include <iostream>
#include <iomanip>   
#include <math.h>
#include <cassert>
#include "hdf5.h"
#include "hdf5_hl.h"

#include "grid_2D_cyln.h"
#include "physical_constants.h"

namespace pc = physical_constants;

using std::string;
using std::cout;
using std::endl;

//------------------------------------------------------------
// initialize the zone geometry from model file
//------------------------------------------------------------
void grid_2D_cyln::read_model_file(ParameterReader* params)
{
  // verbocity
  int my_rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
  const int verbose = (my_rank == 0);
  hsize_t     dims[3];
  double dr[2];
  
  // open up the model file, complaining if it fails to open
  string model_file = params->getScalar<string>("model_file");

  // open hdf5 file 
  hid_t file_id = H5Fopen (model_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  herr_t status;

  // get dimensions
  status = H5LTget_dataset_info(file_id,"/comp",dims, NULL, NULL);
  status = H5LTread_dataset_double(file_id,"/dr",dr);

  nx_     = dims[0];
  nz_     = dims[1];
  n_elems = dims[2];
  dx_ = dr[0];
  dz_ = dr[1];
  n_zones = nz_*nx_;
  z.resize(n_zones);
  vol_.resize(n_zones);

  int *etmp = new int[n_elems];
  status = H5LTread_dataset_int(file_id,"/Z",etmp);
  for (int k=0;k<n_elems;k++) elems_Z.push_back(etmp[k]);
  status = H5LTread_dataset_int(file_id,"/A",etmp);
  for (int k=0;k<n_elems;k++) elems_A.push_back(etmp[k]);
  delete [] etmp;

  double *tmp = new double[n_zones];
  // read density
  status = H5LTread_dataset_double(file_id,"/rho",tmp);
  for (int i=0; i < n_zones; i++) z[i].rho = tmp[i];
  // read temperature
  status = H5LTread_dataset_double(file_id,"/temp",tmp);
  for (int i=0; i < n_zones; i++) z[i].T_gas = tmp[i];
  // read vx
  status = H5LTread_dataset_double(file_id,"/vx",tmp);
  for (int i=0; i < n_zones; i++) z[i].v[0] = tmp[i];
  // read vz
  status = H5LTread_dataset_double(file_id,"/vz",tmp);
  for (int i=0; i < n_zones; i++) z[i].v[2] = tmp[i];
  // read erad
  status = H5LTread_dataset_double(file_id,"/erad",tmp);
  for (int i=0; i < n_zones; i++) z[i].e_rad = tmp[i];
  delete [] tmp;

  // get mass fractions
  double *ctmp = new double[n_zones*n_elems];
  status = H5LTread_dataset_double(file_id,"/comp",ctmp);
  int cnt = 0;
  for (int i=0; i < n_zones; i++)
    for (int k=0; k < n_elems;  k++)
    {
      z[i].X_gas.push_back(ctmp[cnt]);
      cnt++;
    }
  delete [] ctmp;

  // Print out properties
  if (verbose)
  {
    std::cout << "# grid: n_zones = " << n_zones << " ; nx = " << nx_ << "; nz = " << nz_ << "\n";
    double totmass = 0, totke = 0, totrad = 0;
    double *elem_mass = new double[n_elems];
    for (int k = 0;k < n_elems; k++) elem_mass[k] = 0;
    int cnt = 0;
    for (int i=0;i<nx_;i++)
      for (int j=0;j<nz_;j++)
      {
        double rx   = (i+1)*dx_;
        double vrsq = z[cnt].v[0]*z[cnt].v[0] + z[cnt].v[2]*z[cnt].v[2];
        vol_[cnt]   = 2*pc::pi*rx*dx_*dz_;
        index_x_.push_back(i);
        index_z_.push_back(j);

        // compute integral quantities
        totmass    += vol_[cnt]*z[cnt].rho;
        totke      += 0.5*vol_[cnt]*z[cnt].rho*vrsq;
        totrad     += vol_[cnt]*z[cnt].e_rad;
        for (int k = 0;k < n_elems; k++) elem_mass[k] += vol_[cnt]*z[cnt].rho*z[cnt].X_gas[k];
        cnt++;
      }
    printf("# mass = %.4e (%.4e Msun)\n",totmass,totmass/pc::m_sun);
    for (int k=0;k<n_elems;k++) {
      cout << "# " << elems_Z[k] << "." << elems_A[k] <<  "\t";
      cout << elem_mass[k] << " (" << elem_mass[k]/pc::m_sun << " Msun)\n"; 
    }
    printf("# kinetic energy   = %.4e\n",totke);
    printf("# radiation energy = %.4e\n",totrad);
    cout << "##############################\n#\n";
  } 
}
    


//************************************************************
// expand the grid
//************************************************************
void grid_2D_cyln::expand(double e) 
{
  /*

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
  */

}

//************************************************************
// Overly simple search to find zone
//************************************************************
int grid_2D_cyln::get_zone(const double *x) const
{
    double p = sqrt(x[0]*x[0] + x[1]*x[1]);
    int ix = floor(p/dx_);
    int iz = floor(x[2]/dz_) + nz_*dz_/2;

    // check if off boundaries
    if (ix >= nx_) return -2;
    if (iz <= nz_) return -2;
    if (iz >= nz_) return -2;

    return ix*nx_ + iz;
}


//************************************************************
// Overly simple search to find zone
//************************************************************
int grid_2D_cyln::get_next_zone(const double *x, const double *D, int i, double r_core, double *l) const
{
 /* 
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
    *l = l_out; //tiny;
  }
  else
  {
    ind = ind_in;
    *l = l_in; //tiny;
  }

  return ind; 
  */
}



//************************************************************
// Write out the file
//************************************************************
void grid_2D_cyln::write_out(int iw, double tt)
{
  
}



//************************************************************
// return volume of zone (precomputed)
//************************************************************
double  grid_2D_cyln::zone_volume(const int i) const
{
  return vol_[i];
}


//************************************************************
// sample a random position within the spherical shell
//************************************************************
void grid_2D_cyln::sample_in_zone
(int i, std::vector<double> ran, double r[3])
{
    int ix = floor(i/nx_);
    int iz = i - ix*nx_;
    double p   = ix*dx_ + dx_*ran[0];
    double phi = 2*pc::pi*ran[1];
    r[0] = p*cos(phi);
    r[1] = p*sin(phi);
    r[2] = iz*dz_ - dz_*iz/2.0 + dz_*ran[2];
}



//************************************************************
// get the velocity vector 
//************************************************************
void grid_2D_cyln::get_velocity(int i, double x[3], double D[3], double v[3], double *dvds)
{
 /*
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
*/

}

// void grid_1D_sphere::get_radial_edges
// (std::vector<double> &r, double &r0, std::vector<double> &v, double &v0) const
// {
//   for (int i=0;i<n_zones;i++)
//   {
//     r[i] = r_out[i];
//     v[i] = z[i].v[0];
//   }
//   r0 = r_out.min;
//   v0 = v_inner_;
// }
// void grid_1D_sphere::set_radial_edges
// (const std::vector<double> r, const double r0, 
// const std::vector<double> v, const double v0) 
// {
//   r_out.min = r0;
//   v_inner_ = v0;
//   for (int i=0;i<n_zones;i++)
//   {
//     r_out[i] = r[i];
//     z[i].v[0] = v[i];

//     // calculate shell volume
//     double r0;
//     if(i==0) r0 = r_out.min;
//     else     r0 = r_out[i-1];
//     vol[i] = 4.0*pc::pi/3.0*(r_out[i]*r_out[i]*r_out[i] - r0*r0*r0);
//   }


//}


