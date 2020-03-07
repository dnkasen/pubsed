#include <cstdlib>
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

#ifdef MPI_PARALLEL
#include "mpi.h"
#endif



namespace pc = physical_constants;

using std::string;
using std::cout;
using std::cerr;
using std::endl;

//------------------------------------------------------------
// initialize the zone geometry from model file
//------------------------------------------------------------
void grid_2D_cyln::read_model_file(ParameterReader* params)
{
  // verbocity
#ifdef MPI_PARALLEL
  int my_rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
  const int verbose = (my_rank == 0);
#else
  const int verbose = 1;
#endif


  // open up the model file, complaining if it fails to open
  string model_file = params->getScalar<string>("model_file");

  // open hdf5 file
  hid_t file_id = H5Fopen (model_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  herr_t status;

  // get time
  double tt[1];
  status = H5LTread_dataset_double(file_id,"/time",tt);
  if (status < 0) if (verbose) std::cerr << "# Grid Err; can't find time" << endl;
  t_now = tt[0];

  // get grid size and dimensions
  hsize_t     dims[3];
  status = H5LTget_dataset_info(file_id,"/comp",dims, NULL, NULL);
  if (status < 0) if (verbose) std::cerr << "# Grid Err; can't find comp" << endl;

  nx_     = dims[0];
  nz_     = dims[1];
  n_elems = dims[2];
  n_zones = nx_*nz_;
  x_out_.resize(nx_);
  z_out_.resize(nz_);
  z.resize(n_zones);
  dx_.resize(nx_);
  dz_.resize(nz_);
  vol_.resize(n_zones);

  int *etmp = new int[n_elems];
  status = H5LTread_dataset_int(file_id,"/Z",etmp);
  if (status < 0) if (verbose) std::cerr << "# Grid Err; can't find Z" << endl;
  for (int k=0;k<n_elems;k++) elems_Z.push_back(etmp[k]);
  status = H5LTread_dataset_int(file_id,"/A",etmp);
  if (status < 0) if (verbose) std::cerr << "# Grid Err; can't find A" << endl;
  for (int k=0;k<n_elems;k++) elems_A.push_back(etmp[k]);
  delete [] etmp;

  // read minima in each direction
  herr_t status_rmin = H5LTfind_dataset(file_id,"rmin");
  if (status_rmin == 1){
    double rmin[2];
    status = H5LTread_dataset_double(file_id,"/rmin",rmin);
    x_out_.min = rmin[0];
    z_out_.min = rmin[1];
  }
  else {if (verbose) std::cerr << "# Grid Err; can't find rmin" << endl;}

  herr_t status_dr = H5LTfind_dataset(file_id,"dr");
  if (status_dr == 1) {
    double dr[2];
    status = H5LTread_dataset_double(file_id,"/dr",dr);
    if (status_rmin == 0){
      if (verbose) std::cerr << "# Setting minimum values to x_min = 0, z_min = -dz*nz/2" << endl;
      x_out_.min = 0;
      z_out_.min = -dr[1]*nz_/2.0;
    }
    for (int i=0; i < nx_; i++) x_out_[i] = x_out_.min + (i+1.)*dr[0];
    for (int k=0; k < nz_; k++) z_out_[k] = z_out_.min + (k+1.)*dr[1];
  }

  // read x bins
  herr_t status_x = H5LTfind_dataset(file_id,"x_out");
  if (status_x == 1) {
    double *btmp = new double[nx_];
    status = H5LTread_dataset_double(file_id,"/x_out",btmp);
    for (int i=0; i < nx_; i++) x_out_[i] = btmp[i];
    delete [] btmp;
  }
  // read z bins
  herr_t status_z = H5LTfind_dataset(file_id,"z_out");
  if (status_z == 1) {
    double *btmp = new double[nz_];
    status = H5LTread_dataset_double(file_id,"/z_out",btmp);
    for (int i=0; i < nz_; i++) z_out_[i] = btmp[i];
    delete [] btmp;
  }

  if ( (status_dr == 0) && ( (status_x == 0) || (status_z == 0) ) ) {
    if (verbose) std::cerr << "# Grid Err; can't find one of the following inputs to define the grid: 1) dr or 2) x_out, z_out" << endl;
  }

  for (int i=0; i < nx_; i++) dx_[i] = x_out_.delta(i);
  for (int i=0; i < nz_; i++) dz_[i] = z_out_.delta(i);

  // read zone properties
  double *tmp = new double[n_zones];
  // read density
  status = H5LTread_dataset_double(file_id,"/rho",tmp);
  if (status < 0) if (verbose) std::cerr << "# Grid Err; can't find rho" << endl;
  for (int i=0; i < n_zones; i++) z[i].rho = tmp[i];
  // read temperature
  status = H5LTread_dataset_double(file_id,"/temp",tmp);
  if (status < 0) if (verbose) std::cerr << "# Grid Err; can't find temp" << endl;
  for (int i=0; i < n_zones; i++) z[i].T_gas = tmp[i];
  // read vx
  status = H5LTread_dataset_double(file_id,"/vx",tmp);
  if (status < 0) if (verbose) std::cerr << "# Grid Err; can't find vx" << endl;
  for (int i=0; i < n_zones; i++) z[i].v[0] = tmp[i];
  // read vz
  status = H5LTread_dataset_double(file_id,"/vz",tmp);
  if (status < 0) if (verbose) std::cerr << "# Grid Err; can't find vz" << endl;
  for (int i=0; i < n_zones; i++) z[i].v[2] = tmp[i];
  // read erad
  status = H5LTread_dataset_double(file_id,"/erad",tmp);
  if (status < 0) if (verbose) std::cerr << "# Grid Err; can't find erad" << endl;
  for (int i=0; i < n_zones; i++) z[i].e_rad = tmp[i];
  // read grey opacity if the user defines a zone-specific grey opacity
  int use_zone_specific_grey_opacity = params->getScalar<int>("opacity_zone_specific_grey_opacity");
  if(use_zone_specific_grey_opacity != 0){
    status = H5LTread_dataset_double(file_id,"/grey_opacity",tmp);
    if (status < 0) if (verbose) std::cerr << "# Grid Err; can't find grey_opacity" << endl;
    for (int i=0; i < n_zones; i++) z[i].zone_specific_grey_opacity = tmp[i];
  }
  // set bulk grey opacity (note: this parameter is set in the param file, not in the hdf5 file)
  // and set total grey opacity
  double bulk_grey_opacity = params->getScalar<double>("opacity_grey_opacity");
  for (int i=0; i < n_zones; i++){
    z[i].bulk_grey_opacity = bulk_grey_opacity;
    z[i].total_grey_opacity = z[i].bulk_grey_opacity + z[i].zone_specific_grey_opacity;
  }
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

  // close HDF5 input file
  H5Fclose (file_id);

  // allocate indexs
  index_x_.resize(n_zones);
  index_z_.resize(n_zones);

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
      double r0 = x_out_.left(i);
      double r1 = x_out_.right(i);
      vol_[cnt] = pc::pi*(r1*r1 - r0*r0)*dz_[j];

      index_x_[cnt] = i;
      index_z_[cnt] = j;

      double vrsq = z[cnt].v[0]*z[cnt].v[0] + z[cnt].v[2]*z[cnt].v[2];

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
    std::cout << "# gridtype: 2D cylindrical\n";
    std::cout << "# n_zones = " << n_zones << "\n";
    std::cout << "# (nx,nz) = (" << nx_ << ", " << nz_ << ")\n";
    // std::cout << "# (dx,dz) = (" << dx_ << ", " << dz_ << ")\n";

    printf("# mass = %.4e (%.4e Msun)\n",totmass,totmass/pc::m_sun);
    for (int k=0;k<n_elems;k++) {
      cout << "# " << elems_Z[k] << "." << elems_A[k] <<  "\t";
      cout << elem_mass[k] << " (" << elem_mass[k]/pc::m_sun << " Msun)\n";
    }
    printf("# kinetic energy   = %.4e\n",totke);
    printf("# radiation energy = %.4e\n",totrad);
    cout << "##############################\n#" << endl;
  }
  delete [] elem_mass;
}



//************************************************************
// Write out the file
//************************************************************
void grid_2D_cyln::write_plotfile(int iw, double tt, int write_mass_fractions)
{
  // get file name
  char zonefile[1000];
  sprintf(zonefile,"plt_%05d.h5",iw);

  // open hdf5 file
  hid_t file_id = H5Fcreate( zonefile, H5F_ACC_TRUNC, H5P_DEFAULT,  H5P_DEFAULT);

  // // print out radial size
  // hsize_t  dims_dr[1]={2};
  // float dr[2];
  // dr[0] = dx_;
  // dr[1] = dz_;
  // H5LTmake_dataset(file_id,"dr",1,dims_dr,H5T_NATIVE_FLOAT,dr);

  // print out x array
  hsize_t  dims_x[1]={(hsize_t)nx_};
  float *xarr = new float[nx_];
  for (int i=0;i<nx_;i++) xarr[i] = x_out_.left(i);
  H5LTmake_dataset(file_id,"r",1,dims_x,H5T_NATIVE_FLOAT,xarr);
  delete [] xarr;

  // print out z array
  hsize_t  dims_z[1]={(hsize_t)nz_};
  float *zarr = new float[nz_];
  for (int i=0;i<nz_;i++) zarr[i] = z_out_.left(i);
  H5LTmake_dataset(file_id,"z",1,dims_z,H5T_NATIVE_FLOAT,zarr);
  delete [] zarr;

  hsize_t  dims_min[1]={1};
  float x0 = x_out_.min;
  H5LTmake_dataset(file_id,"x_min",1,dims_min,H5T_NATIVE_FLOAT,&x0);
  float z0 = z_out_.min;
  H5LTmake_dataset(file_id,"z_min",1,dims_min,H5T_NATIVE_FLOAT,&z0);

  hsize_t  dims_g[2]={(hsize_t) nx_,(hsize_t) nz_};
  write_hdf5_plotfile_zones(file_id, dims_g, 2, tt);
  write_integrated_quantities(iw,tt);

  // close HDF5 input file
  H5Fclose (file_id);
}



//************************************************************
// expand the grid
//************************************************************
void grid_2D_cyln::expand(double e)
{
  for (int i=0; i < nx_; i++){
    x_out_[i] *= e;
    dx_[i] = dx_[i]*e;
  }
  for (int k=0; k < nz_; k++){
    z_out_[k] *= e;
    dz_[k] = dz_[k]*e;
  }
  x_out_.min *= e;
  z_out_.min *= e;

  for (int i=0; i < n_zones; i++) vol_[i] = vol_[i]*e*e*e;
}

//************************************************************
// Overly simple search to find zone
//************************************************************
int grid_2D_cyln::get_zone(const double *x) const
{
  // double p = sqrt(x[0]*x[0] + x[1]*x[1]);
  // int ix = floor(p/dx_);
  // int iz = floor((x[2] + zcen_)/dz_);

  // // check if off boundaries
  // if (ix >= nx_) return -2;
  // if (iz <    0) return -2;
  // if (iz >= nz_) return -2;

  // return ix*nz_ + iz;

  double p = sqrt(x[0]*x[0] + x[1]*x[1]);

  if (p < x_out_.min) return -2;
  if (p > x_out_[nx_-1]) return -2;
  if (x[2] < z_out_.min) return -2;
  if (x[2] > z_out_[nz_-1]) return -2;

  int i = x_out_.locate_within_bounds(p);
  int k = z_out_.locate_within_bounds(x[2]);

  int ind =  i*nz_ + k;
  return ind;
}


//************************************************************
// Find distance to next zone along path
//************************************************************
int grid_2D_cyln::get_next_zone
(const double *x, const double *D, int i, double r_core, double *l) const
{
  // impact parameter of particle position
  double  p = sqrt(x[0]*x[0] + x[1]*x[1]);
  // z position and direction vector
  double  z = x[2];
  double Dz = D[2];

  int ix = index_x_[i];
  int iz = index_z_[i];

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
    double zt = z_out_.right(iz) + dz_[iz]*tiny;
    lz   = (zt - z)/Dz;
    //std::cout << "iz_up = "<< iz << "; Dz = " << Dz << "; zt = " << zt << "; z = " << z << "; lz = " << lz << "\n";
    d_iz = 1;
  }
  else
  {
    // down interface
    double zt = z_out_.left(iz) - dz_[iz]*tiny;
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
    pt = x_out_.right(ix) + dx_[ix]*tiny;
    c  = p*p - pt*pt;
    det = b*b - 4*a*c;
    if (det < 0) lp_out =  std::numeric_limits<double>::infinity();
    else lp_out = (-1.0*b + sqrt(det))/(2*a);
    if (lp_out < 0) lp_out =  std::numeric_limits<double>::infinity();

    // inner interface
    pt = x_out_.left(ix) - dx_[ix]*tiny;
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
// sample a random position within the annulus weighted by volume
//************************************************************
void grid_2D_cyln::sample_in_zone(int i, std::vector<double> ran, double r[3])
{
  int ix = index_x_[i];
  int iz = index_z_[i];

  double phi = 2*pc::pi*ran[1];
  double p_in = x_out_.left(ix);
  double p_out = x_out_.right(ix);
  double p_samp = sqrt( p_in*p_in + ran[0]*( p_out*p_out-p_in*p_in ) );

  r[0] = p_samp*cos(phi);
  r[1] = p_samp*sin(phi);
  r[2] = z_out_.left(iz) + ran[2]*dz_[iz];
}

//************************************************************
// get the velocity vector
//************************************************************
void grid_2D_cyln::get_velocity(int i, double x[3], double D[3], double v[3], double *dvds)
{

  if (use_homologous_velocities_ == 1) {
    v[0] = x[0]/t_now;
    v[1] = x[1]/t_now;
    v[2] = x[2]/t_now;
    *dvds = 1.0/t_now;
  }
  else {
    double p = sqrt(x[0]*x[0] + x[1]*x[1]);
    v[0] = z[i].v[0]*x[0]/p;
    v[1] = z[i].v[0]*x[1]/p;
    v[2] = z[i].v[2];

    if (p == 0) {
      v[0] = 0;
      v[1] = 0;
    }

    *dvds = 0;
  }

}

void grid_2D_cyln::writeCheckpointGrid(std::string fname) {
  if (my_rank == 0) {
    /* Write out geometry-independent quantities */
    writeCheckpointGeneralGrid(fname);
    hsize_t single_val = 1;
    
    createDataset(fname, "grid", "nx", 1, &single_val, H5T_NATIVE_INT);
    createDataset(fname, "grid", "nz", 1, &single_val, H5T_NATIVE_INT);
    createDataset(fname, "grid", "dx", 1, &single_val, H5T_NATIVE_DOUBLE);
    createDataset(fname, "grid", "dz", 1, &single_val, H5T_NATIVE_DOUBLE);
    createDataset(fname, "grid", "zcen", 1, &single_val, H5T_NATIVE_DOUBLE);

    writeSimple(fname, "grid", "nx", &nx_, H5T_NATIVE_INT);
    writeSimple(fname, "grid", "nz", &nz_, H5T_NATIVE_INT);

    x_out_.writeCheckpoint(fname, "grid", "x_out");
    z_out_.writeCheckpoint(fname, "grid", "z_out");

    writeVector(fname, "grid", "dx", dx_, H5T_NATIVE_DOUBLE);
    writeVector(fname, "grid", "dz", dz_, H5T_NATIVE_DOUBLE);
    writeVector(fname, "grid", "index_x", index_x_, H5T_NATIVE_INT);
    writeVector(fname, "grid", "index_z", index_z_, H5T_NATIVE_INT);
    writeVector(fname, "grid", "vol", vol_, H5T_NATIVE_DOUBLE);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

void grid_2D_cyln::readCheckpointGrid(std::string fname, bool test) {
  for (int rank = 0; rank < nproc; rank++) {
    if (my_rank == rank) {
      readCheckpointGeneralGrid(fname, test);
      readSimple(fname, "grid", "nx", &nx_new_, H5T_NATIVE_INT);
      readSimple(fname, "grid", "nz", &nz_new_, H5T_NATIVE_INT);

      x_out_new_.readCheckpoint(fname, "grid", "x_out");
      z_out_new_.readCheckpoint(fname, "grid", "z_out");

      readVector(fname, "grid", "dx", dx_new_, H5T_NATIVE_DOUBLE);
      readVector(fname, "grid", "dz", dz_new_, H5T_NATIVE_DOUBLE);
      readVector(fname, "grid", "index_x", index_x_new_, H5T_NATIVE_INT);
      readVector(fname, "grid", "index_z", index_z_new_, H5T_NATIVE_INT);
      readVector(fname, "grid", "vol", vol_new_, H5T_NATIVE_DOUBLE);

      if (not test) {
        nx_ = nx_new_;
        nz_ = nz_new_;

        x_out_ = x_out_new_;
        z_out_ = z_out_new_;

        index_x_ = index_x_new_;
        index_z_ = index_z_new_;
        vol_ = vol_new_;
        dx_ = dx_new_;
        dz_ = dz_new_;
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
}


void grid_2D_cyln::restartGrid(ParameterReader* params) {
#ifdef MPI_PARALLEL
  int my_rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
  const int verbose = (my_rank == 0);
#else
  const int verbose = 1;
#endif
  string restart_file = params->getScalar<string>("run_restart_file");

  // geometry of model
  if(params->getScalar<string>("grid_type") != "grid_2D_cyln")
  {
    if (verbose) cerr << "Err: grid_type param disagrees with the model file" << endl;
    exit(4);
  }
  if (verbose) {
    cout << "# model file = " << restart_file << "\n";
    cout << "# Model is a 2D_sphere\n"; }

  readCheckpointGrid(restart_file);
  readCheckpointZones(restart_file);
}
