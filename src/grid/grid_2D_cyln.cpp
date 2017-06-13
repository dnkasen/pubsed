#include <cstdlib>
#include <mpi.h>
#include <fstream>
#include <iostream>
#include <iomanip>   
#include <math.h>
#include <cassert>
#include <limits>

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
  
  // open up the model file, complaining if it fails to open
  string model_file = params->getScalar<string>("model_file");

  // open hdf5 file 
  hid_t file_id = H5Fopen (model_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  herr_t status;

  // get time
  double tt[1];
  status = H5LTread_dataset_double(file_id,"/time",tt);
  if (status < 0) if (verbose) std::cout << "# Grid Err; can't find time\n";
  t_now = tt[0];

  // get grid size and dimensions
  hsize_t     dims[3];
  double dr[2];
  status = H5LTget_dataset_info(file_id,"/comp",dims, NULL, NULL);
  if (status < 0) if (verbose) std::cout << "# Grid Err; can't find comp\n";
  status = H5LTread_dataset_double(file_id,"/dr",dr);
  if (status < 0) if (verbose) std::cout << "# Grid Err; can't find dr\n";

  nx_     = dims[0];
  nz_     = dims[1];
  n_elems = dims[2];
  dx_ = dr[0];
  dz_ = dr[1];
  n_zones = nz_*nx_;
  zcen_ = dz_*nz_/2.0;
  z.resize(n_zones);
  vol_.resize(n_zones);

  int *etmp = new int[n_elems];
  status = H5LTread_dataset_int(file_id,"/Z",etmp);
  if (status < 0) if (verbose) std::cout << "# Grid Err; can't find Z\n";
  for (int k=0;k<n_elems;k++) elems_Z.push_back(etmp[k]);
  status = H5LTread_dataset_int(file_id,"/A",etmp);
  if (status < 0) if (verbose) std::cout << "# Grid Err; can't find A\n";
  for (int k=0;k<n_elems;k++) elems_A.push_back(etmp[k]);
  delete [] etmp;

  double *tmp = new double[n_zones];
  // read density
  status = H5LTread_dataset_double(file_id,"/rho",tmp);
  if (status < 0) if (verbose) std::cout << "# Grid Err; can't find rho\n";
  for (int i=0; i < n_zones; i++) z[i].rho = tmp[i];
  // read temperature
  status = H5LTread_dataset_double(file_id,"/temp",tmp);
  if (status < 0) if (verbose) std::cout << "# Grid Err; can't find temp\n";
  for (int i=0; i < n_zones; i++) z[i].T_gas = tmp[i];
  // read vx
  status = H5LTread_dataset_double(file_id,"/vx",tmp);
  if (status < 0) if (verbose) std::cout << "# Grid Err; can't find vx\n";
  for (int i=0; i < n_zones; i++) z[i].v[0] = tmp[i];
  // read vz
  status = H5LTread_dataset_double(file_id,"/vz",tmp);
  if (status < 0) if (verbose) std::cout << "# Grid Err; can't find vz\n";
  for (int i=0; i < n_zones; i++) z[i].v[2] = tmp[i];
  // read erad
  status = H5LTread_dataset_double(file_id,"/erad",tmp);
  if (status < 0) if (verbose) std::cout << "# Grid Err; can't find erad\n";
  for (int i=0; i < n_zones; i++) z[i].e_rad = tmp[i];
  delete [] tmp;

  // get mass fractions
  double *ctmp = new double[n_zones*n_elems];
  status = H5LTread_dataset_double(file_id,"/comp",ctmp);
  int cnt = 0;
  for (int i=0; i < n_zones; i++)
  {
    z[i].X_gas.resize(n_elems);
    for (int k=0; k < n_elems;  k++)
    {
      z[i].X_gas[k] = ctmp[cnt];
      cnt++;
    }
  }
  delete [] ctmp;

  //---------------------------------------------------
  // Calculate volume, indices, model properties
  //---------------------------------------------------
  double totmass = 0, totke = 0, totrad = 0;
  double *elem_mass = new double[n_elems];
  for (int k = 0;k < n_elems; k++) elem_mass[k] = 0;
  cnt = 0;
  for (int i=0;i<nx_;++i)
    for (int j=0;j<nz_;++j)
    {
      double r0   = i*dx_;
      double r1   = (i+1)*dx_;
      double vrsq = z[cnt].v[0]*z[cnt].v[0] + z[cnt].v[2]*z[cnt].v[2];
      vol_[cnt]   = pc::pi*(r1*r1 - r0*r0)*dz_;
      index_x_.push_back(i);
      index_z_.push_back(j);

      // compute integral quantities
      totmass    += vol_[cnt]*z[cnt].rho;
      totke      += 0.5*vol_[cnt]*z[cnt].rho*vrsq;
      totrad     += vol_[cnt]*z[cnt].e_rad;
      for (int k = 0;k < n_elems; k++) elem_mass[k] += vol_[cnt]*z[cnt].rho*z[cnt].X_gas[k];
      cnt++;
    }

  //--------------------------------------------------- 
  // Printout model properties
  //---------------------------------------------------
  if (verbose)
  {
    std::cout << "# grid: n_zones = " << n_zones << " ; nx = " << nx_ << "; nz = " << nz_ << "\n";
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
// Write out the file
//************************************************************
void grid_2D_cyln::write_out(int iw, double tt)
{
  // get file name
  char zonefile[1000];
  sprintf(zonefile,"zone_%05d.h5",iw);

  // open hdf5 file
  hid_t file_id = H5Fcreate( zonefile, H5F_ACC_TRUNC, H5P_DEFAULT,  H5P_DEFAULT);

  // print out time
  hsize_t  dims_t[1]={1};
  float time_a[1];
  time_a[0] = tt;
  H5LTmake_dataset(file_id,"time",1,dims_t,H5T_NATIVE_FLOAT,time_a);

  // print out radial size
  hsize_t  dims_dr[1]={2};
  float dr[2];
  dr[0] = dx_;
  dr[1] = dz_;
  H5LTmake_dataset(file_id,"dr",1,dims_dr,H5T_NATIVE_FLOAT,dr);

  // print out z array
  hsize_t  dims_z[1]={(hsize_t)nz_};
  float *zarr = new float[nz_];
  for (int i=0;i<nz_;i++) zarr[i] = i*dz_ - zcen_;
  H5LTmake_dataset(file_id,"z",1,dims_z,H5T_NATIVE_FLOAT,zarr);
  delete [] zarr;

  // print out x array
  hsize_t  dims_x[1]={(hsize_t)nx_};
  float *xarr = new float[nx_];
  for (int i=0;i<nx_;i++) xarr[i] = i*dx_;
  H5LTmake_dataset(file_id,"r",1,dims_x,H5T_NATIVE_FLOAT,xarr);
  delete [] xarr;


  // print out zone arrays
  hsize_t  dims_g[2]={(hsize_t) nx_,(hsize_t) nz_};
  float *arr = new float[n_zones];

  // print out rho
  for (int i=0;i<n_zones;++i) arr[i] = z[i].rho; 
  H5LTmake_dataset(file_id,"rho",2,dims_g,H5T_NATIVE_FLOAT,arr);

  // print out T_rad
  for (int i=0;i<n_zones;++i) arr[i] = pow(z[i].e_rad/pc::a,0.25);
  H5LTmake_dataset(file_id,"T_rad",2,dims_g,H5T_NATIVE_FLOAT,arr);

  delete [] arr;

}



//************************************************************
// expand the grid
//************************************************************
void grid_2D_cyln::expand(double e) 
{
  dx_ *= e;
  dz_ *= e;
  zcen_ = dz_*nz_/2.0;

  // recalculate shell volume
  for (int i=0;i<n_zones;i++) 
  {
    vol_[i] = vol_[i]*e*e*e;
  }

}

//************************************************************
// Overly simple search to find zone
//************************************************************
int grid_2D_cyln::get_zone(const double *x) const
{
  double p = sqrt(x[0]*x[0] + x[1]*x[1]);
  int ix = floor(p/dx_);
  int iz = floor((x[2] + zcen_)/dz_);

  // check if off boundaries
  if (ix >= nx_) return -2;
  if (iz <    0) return -2;
  if (iz >= nz_) return -2;

  return ix*nz_ + iz;
}


//************************************************************
// Overly simple search to find zone
//************************************************************
int grid_2D_cyln::get_next_zone
(const double *x, const double *D, int i, double r_core, double *l) const
{
  // impact parameter of particle position
  double  p = sqrt(x[0]*x[0] + x[1]*x[1]);
  // z position and direction vector
  double  z = x[2];
  double Dz = D[2];

  int ix = index_x_[i]; //floor(p/dx_);
  int iz = index_z_[i]; //floor(z/dz_); // note this index can be negative

  double lp,lz;
  int d_ip = 0;
  int d_iz = 0;

  //std::cout << "p: " << ix*dx_ << " " << p << " " << (ix+1)*dx_ << "\n";
  //std::cout << "z: " << iz*dz_ - zcen_ << " " << z << " " << (iz+1)*dz_ - zcen_ << "\n";
  //std::cout << "i: " << i << " " << " ix: " << index_x_[i] << " iz:" << index_z_[i] << " ";
  //std::cout << "ii: " << get_zone(x) << "\n";

  // tiny offset so we don't land exactly on boundaries
  double tiny = 1e-10;
  // distance to z interface
  if (Dz > 0)
  {
    // up interface
    double zt = dz_*(iz + 1) - zcen_ + dz_*tiny;
    lz   = (zt - z)/Dz;
    //std::cout << "iz_up = "<< iz << "; Dz = " << Dz << "; zt = " << zt << "; z = " << z << "; lz = " << lz << "\n";
    d_iz = 1;
  }
  else
  {
    // down interface
    double zt = dz_*iz - zcen_ - dz_*tiny;
    lz   = (zt - z)/Dz;
    if (Dz == 0) lz = std::numeric_limits<double>::infinity();
    //std::cout << "Dz_dn = " << Dz << ";zt = " << zt << "; z = " << z << "; lz = " << lz << "\n";
    d_iz = -1;
  }

  // distance to  p interfaces (annulus)
  double pt,c,lp_out,lp_in,det;
  double a = D[0]*D[0] + D[1]*D[1];
  double b = 2*(x[0]*D[0] + x[1]*D[1]);
  
  if (a == 0) lp = std::numeric_limits<double>::infinity();
  else
  {
    // outer interface
    pt = dx_*(ix + 1) + dx_*tiny;
    c  = p*p - pt*pt;
    det = b*b - 4*a*c;
    if (det < 0) lp_out =  std::numeric_limits<double>::infinity();
    else lp_out = (-1.0*b + sqrt(det))/(2*a);
    if (lp_out < 0) lp_out =  std::numeric_limits<double>::infinity();
 
    // inner interface
    pt = dx_*(ix) - dx_*tiny;
    c = p*p - pt*pt; 
    det = b*b - 4*a*c;
    if (det < 0) lp_in =  std::numeric_limits<double>::infinity();
    else lp_in = (-1.0*b - sqrt(det))/(2*a);
    if (lp_in < 0) lp_in =  std::numeric_limits<double>::infinity();
    if (ix == 0) lp_in   =  std::numeric_limits<double>::infinity();

    if (lp_in < lp_out) 
    {
     lp  = lp_in;
     d_ip = -1;
    }
    else
    {
     lp = lp_out;
     d_ip = 1;
    }
  }

  int new_iz = index_z_[i];
  int new_ip = index_x_[i];
  // find smallest interface
  if (lz < lp)
  {
    *l = lz;
    new_iz += d_iz;
    //std::cout << "step z " << new_iz << " " << d_iz << " " << lz << "\n";
  }
  else
  {
    *l = lp;
    new_ip += d_ip;
    //std::cout << "step p " << new_ip << " " << d_ip << " " << lp << "\n";
  }

  if (isnan(*l)) 
  {
  //   std::cout << "step p " << new_ip << " " << d_ip << " " << lp << "\n";
   // std::cout << "step z " << new_iz << " " << d_iz << " " << lz << "\n";

  }

  // escaped 
  if ((new_iz < 0)||(new_iz >= nz_)||(new_ip >= nx_)) return -2;

  return new_ip*nz_ + new_iz;

; /* 
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
  double p   = index_x_[i]*dx_ + dx_*ran[0];
  double phi = 2*pc::pi*ran[1];
  r[0] = p*cos(phi);
  r[1] = p*sin(phi);
  r[2] = index_z_[i]*dz_ - zcen_ + dz_*ran[2];
}



//************************************************************
// get the velocity vector 
//************************************************************
void grid_2D_cyln::get_velocity(int i, double x[3], double D[3], double v[3], double *dvds)
{
  v[0] = 0;
  v[1] = 0;
  v[2] = 0;
  *dvds = 0;

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


