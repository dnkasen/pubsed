#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cassert>
#include "Lua.h"
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
  if (status < 0) if (verbose) std::cout << "# Grid Err; can't find time\n";
  t_now = tt[0];

  // get grid size and dimensions
  hsize_t     dims[4];
  double dr[3], rmin[3];
  status = H5LTget_dataset_info(file_id,"/comp",dims, NULL, NULL);
  if (status < 0) if (verbose) std::cout << "# Grid Err; can't find comp\n";
  status = H5LTread_dataset_double(file_id,"/dr",dr);
  if (status < 0) if (verbose) std::cout << "# Grid Err; can't find dr\n";
  status = H5LTread_dataset_double(file_id,"/rmin",rmin);
  if (status < 0) if (verbose) std::cout << "# Grid Err; can't find rmin\n";

  nx_     = dims[0];
  ny_     = dims[1];
  nz_     = dims[2];
  n_elems = dims[3];
  dx_ = dr[0];
  dy_ = dr[1];
  dz_ = dr[2];
  x0_ = rmin[0];
  y0_ = rmin[1];
  z0_ = rmin[2];
  n_zones = nz_*ny_*nx_;
  z.resize(n_zones);

  // volume
  vol_ = dx_*dy_*dz_;

  // read elements Z and A
  int *etmp = new int[n_elems];
  status = H5LTread_dataset_int(file_id,"/Z",etmp);
  if (status < 0) if (verbose) std::cout << "# Grid Err; can't find Z\n";
  for (int k=0;k<n_elems;k++) elems_Z.push_back(etmp[k]);
  status = H5LTread_dataset_int(file_id,"/A",etmp);
  if (status < 0) if (verbose) std::cout << "# Grid Err; can't find A\n";
  for (int k=0;k<n_elems;k++) elems_A.push_back(etmp[k]);
  delete [] etmp;

  // read zone properties
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
  // read vy
  status = H5LTread_dataset_double(file_id,"/vy",tmp);
  if (status < 0) if (verbose) std::cout << "# Grid Err; can't find vx\n";
  for (int i=0; i < n_zones; i++) z[i].v[1] = tmp[i];
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

  // close HDF5 input file
  H5Fclose (file_id);

  // allocate indexs
  index_x_ = new int[n_zones];
  index_y_ = new int[n_zones];
  index_z_ = new int[n_zones];

  //---------------------------------------------------
  // Calculate volume, indices, model properties
  //---------------------------------------------------
  double totmass = 0, totke = 0, totrad = 0;
  double *elem_mass = new double[n_elems];
  for (int l = 0;l < n_elems; l++) elem_mass[l] = 0;
  cnt = 0;
  for (int i=0;i<nx_;++i)
    for (int j=0;j<ny_;++j)
      for (int k=0;k<nz_;++k)
      {
        index_x_[cnt] = i;
        index_y_[cnt] = j;
        index_z_[cnt] = k;

        double vrsq = z[cnt].v[0]*z[cnt].v[0] + z[cnt].v[1]*z[cnt].v[1] + z[cnt].v[2]*z[cnt].v[2];

        // compute integral quantities
        totmass    += vol_*z[cnt].rho;
        totke      += 0.5*vol_*z[cnt].rho*vrsq;
        totrad     += vol_*z[cnt].e_rad;
        for (int l = 0;l < n_elems; l++) elem_mass[l] += vol_*z[cnt].rho*z[cnt].X_gas[l];
        cnt++;
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
    std::cout << "# (dx,dy,dz) = (";
    std::cout << dx_ << ", " << dy_ << ", " << dz_ << ")\n";
    std::cout << "# (x0,y0,z0) = (";
    std::cout << x0_ << ", " << y0_ << ", " << z0_ << ")\n";
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
void grid_3D_cart::write_plotfile(int iw, double tt, int write_mass_fractions)
{
	// get file name
	char plotfilename[1000];
	sprintf(plotfilename,"plt_%05d.h5",iw);

	// open hdf5 file
	hid_t file_id = H5Fcreate(plotfilename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	// print out zone sizes
	hsize_t  dims_dr[1]={3};
	float dr[3];
	dr[0] = dx_;
  dr[1] = dy_;
  dr[2] = dz_;
  H5LTmake_dataset(file_id,"dr",1,dims_dr,H5T_NATIVE_FLOAT,dr);

	// print out x array
	hsize_t  dims_x[1]={(hsize_t)nx_};
	float *xarr = new float[nx_];
	for (int i=0;i<nx_;i++) xarr[i] = i*dx_ + x0_;
	H5LTmake_dataset(file_id,"x",1,dims_x,H5T_NATIVE_FLOAT,xarr);
  delete [] xarr;

  // print out y array
	hsize_t  dims_y[1]={(hsize_t)ny_};
	float *yarr = new float[ny_];
	for (int i=0;i<ny_;i++) yarr[i] = i*dy_ + y0_;
	H5LTmake_dataset(file_id,"y",1,dims_y,H5T_NATIVE_FLOAT,yarr);
  delete [] yarr;

  // print out z array
	hsize_t  dims_z[1]={(hsize_t)nz_};
	float *zarr = new float[nz_];
	for (int i=0;i<nz_;i++) zarr[i] = i*dz_ + z0_;
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
  dx_ *= e;
  dy_ *= e;
  dz_ *= e;
  x0_ *= e;
  y0_ *= e;
  z0_ *= e;
  vol_ = vol_*e*e*e;
}



//------------------------------------------------------------
// Overly simple search to find zone
//------------------------------------------------------------
int grid_3D_cart::get_zone(const double *x) const
{
  int i = floor((x[0]-x0_)/dx_);
  int j = floor((x[1]-y0_)/dy_);
  int k = floor((x[2]-z0_)/dz_);

  // check for off grid
  if ((i < 0)||(i > nx_-1)) return -2;
  if ((j < 0)||(j > ny_-1)) return -2;
  if ((k < 0)||(k > nz_-1)) return -2;

  int ind =  i*ny_*nz_ + j*nz_ + k;
  return ind;
}


//************************************************************
// Find distance to next zone along path
//************************************************************
int grid_3D_cart::get_next_zone
(const double *x, const double *D, int i, double r_core, double *l) const
{
  // tiny offset so we don't land exactly on boundaries
  double tiny = 1e-10;

  double bn, len[3];

  //---------------------------------
  // distance to x interfaces
  //---------------------------------
  if (D[0] > 0)
    bn = dx_*(index_x_[i] + 1 + tiny) + x0_;
  else
    bn = dx_*(index_x_[i] + tiny) + x0_;
  //std::cout << bn << " ";
  len[0] = (bn - x[0])/D[0];

  //---------------------------------
  // distance to y interfaces
  //---------------------------------
  if (D[1] > 0)
    bn = dy_*(index_y_[i] + 1 + tiny) + y0_;
  else
    bn = dy_*(index_y_[i] + tiny) + y0_;
  len[1] = (bn - x[1])/D[1];
//  std::cout << bn << " ";

   //---------------------------------
  // distance to z interfaces
  //---------------------------------
  if (D[2] > 0)
    bn = dz_*(index_z_[i] + 1 + tiny) + z0_;
  else
    bn = dz_*(index_z_[i] + tiny) + z0_;
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
  return vol_;
}

//------------------------------------------------------------
// sample a random position within the cubical cell
//------------------------------------------------------------
void grid_3D_cart::sample_in_zone
(const int i, const std::vector<double> ran,double r[3])
{
  r[0] = x0_ + (index_x_[i] + ran[0])*dx_;
  r[1] = y0_ + (index_y_[i] + ran[1])*dy_;
  r[2] = z0_ + (index_z_[i] + ran[2])*dz_;
}

//************************************************************
// expand the grid
//************************************************************
void grid_3D_cart::coordinates(int i,double r[3])
{
  r[0] = x0_ + (index_x_[i] + 0.5)*dx_;
  r[1] = y0_ + (index_y_[i] + 0.5)*dy_;
  r[2] = z0_ + (index_z_[i] + 0.5)*dz_;
}


//------------------------------------------------------------
// get the velocity vector
//------------------------------------------------------------
void grid_3D_cart::get_velocity(int i, double x[3], double D[3], double v[3], double *dvds)
{
  // may want to interpolate here
  v[0] = z[i].v[0];
  v[1] = z[i].v[1];
  v[2] = z[i].v[2];

  // debug -- need to calculate this
  *dvds = 0;
}
