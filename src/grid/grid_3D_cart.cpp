#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cassert>
#include "grid_3D_cart.h"
#include "physical_constants.h"

#include "hdf5.h"
#include "hdf5_hl.h"

#ifdef MPI_PARALLEL
#include "mpi.h"
#endif

namespace pc = physical_constants;

using std::string;
using std::cout;
using std::cerr;
using std::endl;

//------------------------------------------------------------
// Read in a cartesian model file
//------------------------------------------------------------
void grid_3D_cart::read_model_file(ParameterReader* params)
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
  hsize_t     dims[4];
  status = H5LTget_dataset_info(file_id,"/comp",dims, NULL, NULL);
  if (status < 0) if (verbose) std::cerr << "# Grid Err; can't find comp" << endl;

  nx_     = dims[0];
  ny_     = dims[1];
  nz_     = dims[2];
  n_elems = dims[3];
  n_zones = nx_*ny_*nz_;
  x_out_.resize(nx_);
  y_out_.resize(ny_);
  z_out_.resize(nz_);
  z.resize(n_zones);
  dx_.resize(nx_);
  dy_.resize(ny_);
  dz_.resize(nz_);
  vol_.resize(n_zones);

  // read elements Z and A
  int *etmp = new int[n_elems];
  status = H5LTread_dataset_int(file_id,"/Z",etmp);
  if (status < 0) if (verbose) std::cerr << "# Grid Err; can't find Z" << endl;
  for (int k=0;k<n_elems;k++) elems_Z.push_back(etmp[k]);
  status = H5LTread_dataset_int(file_id,"/A",etmp);
  if (status < 0) if (verbose) std::cerr << "# Grid Err; can't find A" << endl;
  for (int k=0;k<n_elems;k++) elems_A.push_back(etmp[k]);
  delete [] etmp;

  double* xtmp;
  double* ytmp;
  double* ztmp;
  double xmin;
  double ymin;
  double zmin;
  double dr[3];
  double rmin[3];
  // Check which grid inputs exist
  herr_t status_x = H5LTfind_dataset(file_id,"x_out");
  herr_t status_y = H5LTfind_dataset(file_id,"y_out");
  herr_t status_z = H5LTfind_dataset(file_id,"z_out");
  status_dr_ = H5LTfind_dataset(file_id,"dr");
  status_rmin_ = H5LTfind_dataset(file_id,"rmin");
  status_xyz_ = (status_x && status_y && status_z);
  // Read in data if exists
  if (status_dr_) {
    status = H5LTread_dataset_double(file_id,"/dr",dr);
  }
  if (status_rmin_) {
    status = H5LTread_dataset_double(file_id,"/rmin",rmin);
    xmin = rmin[0];
    ymin = rmin[1];
    zmin = rmin[2];
  }
  else {
    if (status_dr_) {
      xmin = -dr[0]*nx_/2.0;
      ymin = -dr[1]*ny_/2.0;
      zmin = -dr[2]*nz_/2.0;
    }
  }
  if (status_xyz_) {
    xtmp = new double[nx_];
    status = H5LTread_dataset_double(file_id,"/x_out",xtmp);
    ytmp = new double[ny_];
    status = H5LTread_dataset_double(file_id,"/y_out",ytmp);
    ztmp = new double[nz_];
    status = H5LTread_dataset_double(file_id,"/z_out",ztmp);
  }
  // Initialize x_out_, y_out_ and z_out_ if possible
  if (status_xyz_ && status_rmin_) {
    x_out_.init(xtmp, nx_, xmin);
    y_out_.init(ytmp, ny_, ymin);
    z_out_.init(ztmp, nz_, zmin);
    delete [] xtmp;
    delete [] ytmp;
    delete [] ztmp;
  }
  else if (status_xyz_ && (!status_rmin_)) {
    std::cerr << "Error: Missing rmin in grid. Exiting." << std::endl;
    delete [] xtmp;
    delete [] ytmp;
    delete [] ztmp;
    exit(99);
  }
  else if (status_dr_ && status_rmin_) {
    double xmax = xmin + dr[0] * nx_;
    double ymax = ymin + dr[1] * ny_;
    double zmax = zmin + dr[2] * nz_;
    x_out_.init(xmin, xmax, dr[0]);
    y_out_.init(ymin, ymax, dr[1]);
    z_out_.init(zmin, zmax, dr[2]);
  }
  else if (status_dr_ && (!status_rmin_)) {
    xmin = -dr[0]*nx_/2.0;
    ymin = -dr[1]*ny_/2.0;
    zmin = -dr[2]*nz_/2.0;
    std::cerr << xmin << " " << ymin << " " << zmin << std::endl;
    double xmax = xmin + dr[0] * nx_;
    double ymax = ymin + dr[1] * ny_;
    double zmax = zmin + dr[2] * nz_;
    x_out_.init(xmin, xmax, dr[0]);
    y_out_.init(ymin, ymax, dr[1]);
    z_out_.init(zmin, zmax, dr[2]);
  }
  else {
    if (verbose) std::cerr << "# Grid Err; can't find one of the following inputs to define the grid: 1) dr or 2) x_out, y_out, z_out" << endl;
    exit(10);
  }

  for (int i=0; i < nx_; i++) dx_[i] = x_out_.delta(i);
  for (int i=0; i < ny_; i++) dy_[i] = y_out_.delta(i);
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
  // read vy
  status = H5LTread_dataset_double(file_id,"/vy",tmp);
  if (status < 0) if (verbose) std::cerr << "# Grid Err; can't find vx" << endl;
  for (int i=0; i < n_zones; i++) z[i].v[1] = tmp[i];
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
  index_y_.resize(n_zones);
  index_z_.resize(n_zones);

  //---------------------------------------------------
  // Calculate volume, indices, model properties
  //---------------------------------------------------
  double totmass = 0, totke = 0, totrad = 0;
  double *elem_mass = new double[n_elems];
  for (int l = 0;l < n_elems; l++) elem_mass[l] = 0;
  cnt = 0;
  for (int i=0;i<nx_;++i)
  {
    for (int j=0;j<ny_;++j)
    {
      for (int k=0;k<nz_;++k)
      {
        vol_[cnt] = dx_[i]*dy_[j]*dz_[k];

        index_x_[cnt] = i;
        index_y_[cnt] = j;
        index_z_[cnt] = k;

        double vrsq = z[cnt].v[0]*z[cnt].v[0] + z[cnt].v[1]*z[cnt].v[1] + z[cnt].v[2]*z[cnt].v[2];

        // compute integral quantities
        totmass    += vol_[cnt]*z[cnt].rho;
        totke      += 0.5*vol_[cnt]*z[cnt].rho*vrsq;
        totrad     += vol_[cnt]*z[cnt].e_rad;
        for (int l = 0;l < n_elems; l++) elem_mass[l] += vol_[cnt]*z[cnt].rho*z[cnt].X_gas[l];
        cnt++;
      }
    }
  }

  //---------------------------------------------------
  // Printout model properties
  //---------------------------------------------------
  if (verbose)
  {
    std::cout << "# gridtype: 3D cartesian\n";
    std::cout << "# n_zones = " << n_zones << "\n";
    std::cout << "# (nx,ny,nz) = (";
    std::cout << nx_ << ", " << ny_ << ", " << nz_ << ")\n";
    // std::cout << "# (dx,dy,dz) = (";
    // std::cout << dx_ << ", " << dy_ << ", " << dz_ << ")\n";
    std::cout << "# (x_min,y_min,z_min) = (";
    std::cout << x_out_.minval() << ", " << y_out_.minval() << ", " << z_out_.minval() << ")\n";
    printf("# mass = %.4e (%.4e Msun)\n",totmass,totmass/pc::m_sun);
    for (int k=0;k<n_elems;k++) {
      cout << "# " << elems_Z[k] << "." << elems_A[k] <<  "\t";
      cout << elem_mass[k] << " (" << elem_mass[k]/pc::m_sun << " Msun)\n";
    }
    printf("# kinetic energy   = %.4e\n",totke);
    printf("# radiation energy = %.4e\n",totrad);
    cout << "##############################\n#" << endl;
  }

}


//************************************************************
// Write out the file
//************************************************************
void grid_3D_cart::write_plotfile(int iw, double tt, int write_mass_fractions)
{
	// get file name
	char plotfilename[1000];
	sprintf(plotfilename,"plt_%05d.h5",iw);

	// open hdf5 file
	hid_t file_id = H5Fcreate(plotfilename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  // If fixed dr's exist, output them. Some tests want it
  if (status_dr_) {
  	hsize_t  dims_dr[1]={3};
  	float dr[3];
  	dr[0] = x_out_.delval();
    dr[1] = y_out_.delval();
    dr[2] = z_out_.delval();
    H5LTmake_dataset(file_id,"dr",1,dims_dr,H5T_NATIVE_FLOAT,dr);
  }

	// print out x array
	hsize_t  dims_x[1]={(hsize_t)nx_};
	float *xarr = new float[nx_];
	for (int i=0;i<nx_;i++) xarr[i] = x_out_.left(i);
	H5LTmake_dataset(file_id,"x",1,dims_x,H5T_NATIVE_FLOAT,xarr);
  delete [] xarr;

  // print out y array
	hsize_t  dims_y[1]={(hsize_t)ny_};
	float *yarr = new float[ny_];
	for (int i=0;i<ny_;i++) yarr[i] = y_out_.left(i);
	H5LTmake_dataset(file_id,"y",1,dims_y,H5T_NATIVE_FLOAT,yarr);
  delete [] yarr;

  // print out z array
	hsize_t  dims_z[1]={(hsize_t)nz_};
	float *zarr = new float[nz_];
	for (int i=0;i<nz_;i++) zarr[i] = z_out_.left(i);
	H5LTmake_dataset(file_id,"z",1,dims_z,H5T_NATIVE_FLOAT,zarr);
  delete [] zarr;

  hsize_t  dims_g[3]={(hsize_t) nx_,(hsize_t) ny_,(hsize_t) nz_};
  write_hdf5_plotfile_zones(file_id, dims_g, 3, tt);
  write_integrated_quantities(iw,tt);

  // Close the output file
  H5Fclose (file_id);

}



//************************************************************
// expand the grid
//************************************************************
void grid_3D_cart::expand(double e)
{
  x_out_.scale(e);
  y_out_.scale(e);
  z_out_.scale(e);

  for (int i=0; i < nx_; i++) dx_[i] = dx_[i]*e;
  for (int j=0; j < ny_; j++) dy_[j] = dy_[j]*e;
  for (int k=0; k < nz_; k++) dz_[k] = dz_[k]*e;

  for (int i=0; i < n_zones; i++) vol_[i] = vol_[i]*e*e*e;
}



//------------------------------------------------------------
// Overly simple search to find zone
//------------------------------------------------------------
int grid_3D_cart::get_zone(const double *x) const
{
  // int i = floor((x[0]-x0_)/dx_);
  // int j = floor((x[1]-y0_)/dy_);
  // int k = floor((x[2]-z0_)/dz_);

  // // check for off grid
  // if ((i < 0)||(i > nx_-1)) return -2;
  // if ((j < 0)||(j > ny_-1)) return -2;
  // if ((k < 0)||(k > nz_-1)) return -2;

  // int ind =  i*ny_*nz_ + j*nz_ + k;
  // return ind;

  if (x[0] < x_out_.minval()) return -2;
  if (x[0] > x_out_[nx_-1]) return -2;
  if (x[1] < y_out_.minval()) return -2;
  if (x[1] > y_out_[ny_-1]) return -2;
  if (x[2] < z_out_.minval()) return -2;
  if (x[2] > z_out_[nz_-1]) return -2;

  int i = x_out_.locate_within_bounds(x[0]);
  int j = y_out_.locate_within_bounds(x[1]);
  int k = z_out_.locate_within_bounds(x[2]);

  int ind =  i*ny_*nz_ + j*nz_ + k;
  return ind;
}


//************************************************************
// Find distance to next zone along path
//************************************************************
int grid_3D_cart::get_next_zone
(const double *x, const double *D, int i, double r_core, double *l) const
{
  int ix = index_x_[i];
  int iy = index_y_[i];
  int iz = index_z_[i];

  // tiny offset so we don't land exactly on boundaries
  double tiny = 1e-10;

  double bn, len[3];

  //---------------------------------
  // distance to x interfaces
  //---------------------------------
  if (D[0] > 0)
    bn = x_out_.right(ix) + dx_[ix]*tiny;
  else
    bn = x_out_.left(ix) - dx_[ix]*tiny;
  //std::cout << bn << " ";
  len[0] = (bn - x[0])/D[0];

  //---------------------------------
  // distance to y interfaces
  //---------------------------------
  if (D[1] > 0)
    bn = y_out_.right(iy) + dy_[iy]*tiny;
  else
    bn = y_out_.left(iy) - dy_[iy]*tiny;
  len[1] = (bn - x[1])/D[1];
//  std::cout << bn << " ";

   //---------------------------------
  // distance to z interfaces
  //---------------------------------
  if (D[2] > 0)
    bn = z_out_.right(iz) + dz_[iz]*tiny;
  else
    bn = z_out_.left(iz) - dz_[iz]*tiny;
  len[2] = (bn - x[2])/D[2];
  //std::cout << bn << "\n";

  // find shortest distance
  int idx=0,idy=0,idz=0;
  if ((len[0] < len[1])&&(len[0] < len[2]))
  {
    if (D[0] < 0) idx = -1;
    else idx = 1;
    *l = len[0];
  }
  else if (len[1] < len[2])
  {
    if (D[1] < 0) idy = -1;
    else idy = 1;
    *l = len[1];
  }
  else
  {
    if (D[2] < 0) idz = -1;
    else idz = 1;
    *l = len[2];
  }

  int new_ix = index_x_[i] + idx;
  int new_iy = index_y_[i] + idy;
  int new_iz = index_z_[i] + idz;

  // check for off grid
  if ((new_ix < 0)||(new_ix > nx_-1)) return -2;
  if ((new_iy < 0)||(new_iy > ny_-1)) return -2;
  if ((new_iz < 0)||(new_iz > nz_-1)) return -2;

  return new_ix*(ny_*nz_) + new_iy*(nz_) + new_iz;
}



//------------------------------------------------------------
// return volume of zone (precomputed)
//------------------------------------------------------------
double grid_3D_cart::zone_volume(const int i) const
{
  return vol_[i];
}

//------------------------------------------------------------
// sample a random position within the cubical cell
//------------------------------------------------------------
void grid_3D_cart::sample_in_zone
(const int i, const std::vector<double> ran,double r[3])
{
  int ix = index_x_[i];
  int iy = index_y_[i];
  int iz = index_z_[i];

  r[0] = x_out_.left(ix) + ran[2]*dx_[ix];
  r[1] = y_out_.left(iy) + ran[2]*dy_[iy];
  r[2] = z_out_.left(iz) + ran[2]*dz_[iz];
}

//************************************************************
// get coordinates of zone center
//************************************************************
void grid_3D_cart::coordinates(int i,double r[3])
{
  int ix = index_x_[i];
  int iy = index_y_[i];
  int iz = index_z_[i];

  r[0] = x_out_.left(ix) + 0.5*dx_[i];
  r[1] = y_out_.left(iy) + 0.5*dy_[i];
  r[2] = z_out_.left(iz) + 0.5*dz_[i];
}

//------------------------------------------------------------
// get the velocity vector
//------------------------------------------------------------
void grid_3D_cart::get_velocity(int i, double x[3], double D[3], double v[3], double *dvds)
{

  if (use_homologous_velocities_ == 1) {
    v[0] = x[0]/t_now;
    v[1] = x[1]/t_now;
    v[2] = x[2]/t_now;
    *dvds = 1.0/t_now;
  }
  else {

    // trilinear interpolation
    // e.g., https://en.wikipedia.org/wiki/Trilinear_interpolation

    // indices along each dimensions for this zone
    int ic[3];
    ic[0] = index_x_[i];
    ic[1] = index_y_[i];
    ic[2] = index_z_[i];

    // put grid data into array form
    // this should probably be stored this way to start...
    int npts[3];
    npts[0] = nx_;
    npts[1] = ny_;
    npts[2] = nz_;
    double rmin[3];
    rmin[0] = x_out_.minval();
    rmin[1] = y_out_.minval();
    rmin[2] = z_out_.minval();
    double del[3];
    del[0] = dx_[ic[0]];
    del[1] = dy_[ic[1]];
    del[2] = dz_[ic[2]];

    // find lattice points to interpolate from
    int    i_lo[3], i_hi[3];
    double r_lo[3];
    for (int j=0;j<3;j++)
    {
      double cen = rmin[j] + (ic[j] + 0.5)*del[j];
      if (x[j] > cen)
      {
        r_lo[j] = cen;
        i_lo[j] = ic[j];
        i_hi[j] = ic[j] + 1;
        if (i_hi[j] >= npts[j]) i_hi[j] = i_lo[j];
      }
      else
      {
        r_lo[j] = rmin[j] + (ic[j] - 0.5)*del[j];
        i_lo[j] = ic[j] - 1;
        i_hi[j] = ic[j];
        if (i_lo[j] < 0) i_lo[j] = i_hi[j];
      }
    }

    // differences in each dimension
    double xd = (x[0] - r_lo[0])/dx_[ic[0]];
    double yd = (x[1] - r_lo[1])/dy_[ic[1]];
    double zd = (x[2] - r_lo[2])/dz_[ic[2]];

    // interpolate each of the 3 velocity components
    for (int j = 0;j < 3; j++)
    {
      // values at 8 lattice points
      double c000 = z[get_index(i_lo[0],i_lo[1],i_lo[2])].v[j];
      double c001 = z[get_index(i_lo[0],i_lo[1],i_hi[2])].v[j];
      double c010 = z[get_index(i_lo[0],i_hi[1],i_lo[2])].v[j];
      double c011 = z[get_index(i_lo[0],i_hi[1],i_hi[2])].v[j];
      double c100 = z[get_index(i_hi[0],i_lo[1],i_lo[2])].v[j];
      double c101 = z[get_index(i_hi[0],i_lo[1],i_hi[2])].v[j];
      double c110 = z[get_index(i_hi[0],i_hi[1],i_lo[2])].v[j];
      double c111 = z[get_index(i_hi[0],i_hi[1],i_hi[2])].v[j];

      // interpolation along x
      double c00 = c000*(1-xd) + c100*xd;
      double c01 = c001*(1-xd) + c101*xd;
      double c10 = c010*(1-xd) + c110*xd;
      double c11 = c011*(1-xd) + c111*xd;

      // interpolation along y
      double c0 = c00*(1-yd) + c10*yd;
      double c1 = c01*(1-yd) + c11*yd;

      // interpolation along z
      double c = c0*(1-zd) + c1*zd;
      v[j] = c;
    }

    // debug -- need to calculate this
    *dvds = 0;
  }
  
}
