#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cassert>
#include "grid_3D_sphere.h"
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
// Read in a spherical model file
//------------------------------------------------------------
void grid_3D_sphere::read_model_file(ParameterReader* params)
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

  nr_       = dims[0];
  ntheta_   = dims[1];
  nphi_     = dims[2];
  n_elems = dims[3];
  n_zones = nr_*ntheta_*nphi_;
  r_out_.resize(nr_);
  theta_out_.resize(ntheta_);
  phi_out_.resize(nphi_);
  z.resize(n_zones);
  dr_.resize(nr_);
  dtheta_.resize(ntheta_);
  dphi_.resize(nphi_);
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

  double* rtmp;
  double* thetatmp;
  double* phitmp;
  double r_minimum;
  double theta_minimum;
  double phi_minimum;
  double dr[3];
  double rmin[1];
  // Check which grid inputs exist
  herr_t status_r = H5LTfind_dataset(file_id,"r_out");
  herr_t status_theta = H5LTfind_dataset(file_id,"theta_out");
  herr_t status_phi = H5LTfind_dataset(file_id,"phi_out");
  status_dr_ = H5LTfind_dataset(file_id,"dr");
  status_rmin_ = H5LTfind_dataset(file_id,"rmin");
  status_rthetaphi_ = (status_r && status_theta && status_phi);
  // Read in data if exists
  if (status_dr_) {
    status = H5LTread_dataset_double(file_id,"/dr",dr);
  }
  if (status_rmin_) {
    status = H5LTread_dataset_double(file_id,"/rmin",rmin);
    r_minimum = rmin[0];
    theta_minimum = 0;
    phi_minimum = 0;
  }
  else {
    if (status_dr_) {
      r_minimum = 0;
      theta_minimum = 0;
      phi_minimum = 0;
    }
  }
  if (status_rthetaphi_) {
    rtmp = new double[nr_];
    status = H5LTread_dataset_double(file_id,"/r_out",rtmp);
    thetatmp = new double[ntheta_];
    status = H5LTread_dataset_double(file_id,"/theta_out",thetatmp);
    phitmp = new double[nphi_];
    status = H5LTread_dataset_double(file_id,"/phi_out",phitmp);
  }
  // Initialize r_out_, theta_out_ and phi_out_ if possible
  if (status_rthetaphi_ && status_rmin_) {
    if (fabs(thetatmp[ntheta_-1] - pc::pi)/pc::pi > 1e-6){
      if (verbose) std::cerr << "# Grid Err; rightmost value of theta grid not equal to pi" << endl;
      exit(99);
    }
    else{
      thetatmp[ntheta_-1] = pc::pi;
    }
    if (fabs(phitmp[nphi_-1] - 2.*pc::pi)/(2.*pc::pi) > 1e-6){
      if (verbose) std::cerr << "# Grid Err; rightmost value of phi grid not equal to 2*pi" << endl;
      exit(99);
    }
    else{
      phitmp[nphi_-1] = 2.*pc::pi;
    }
    r_out_.init(rtmp, nr_, r_minimum);
    theta_out_.init(thetatmp, ntheta_, theta_minimum);
    phi_out_.init(phitmp, nphi_, phi_minimum);
    delete [] rtmp;
    delete [] thetatmp;
    delete [] phitmp;
  }
  else if (status_rthetaphi_ && (!status_rmin_)) {
    std::cerr << "Error: Missing rmin in grid. Exiting." << std::endl;
    delete [] rtmp;
    delete [] thetatmp;
    delete [] phitmp;
    exit(99);
  }
  else if (status_dr_ && status_rmin_) {
    double r_maximum = r_minimum + dr[0] * nr_;
    double theta_maximum = theta_minimum + dr[1] * ntheta_;
    double phi_maximum = phi_minimum + dr[2] * nphi_;
    r_out_.init(r_minimum, r_maximum, dr[0]);
    theta_out_.init(theta_minimum, theta_maximum, dr[1]);
    phi_out_.init(phi_minimum, phi_maximum, dr[2]);
  }
  else if (status_dr_ && (!status_rmin_)) {
    r_minimum = 0;
    theta_minimum = 0;
    phi_minimum = 0;
    std::cerr << r_minimum << " " << theta_minimum << " " << phi_minimum << std::endl;
    double r_maximum = r_minimum + dr[0] * nr_;
    double theta_maximum = theta_minimum + dr[1] * ntheta_;
    double phi_maximum = phi_minimum + dr[2] * nphi_;
    r_out_.init(r_minimum, r_maximum, dr[0]);
    theta_out_.init(theta_minimum, theta_maximum, dr[1]);
    phi_out_.init(phi_minimum, phi_maximum, dr[2]);
  }
  else {
    if (verbose) std::cerr << "# Grid Err; can't find one of the following inputs to define the grid: 1) dr or 2) r_out, theta_out, phi_out" << endl;
    exit(10);
  }

  for (int i=0; i < nr_; i++) dr_[i] = r_out_.delta(i);
  for (int i=0; i < ntheta_; i++) dtheta_[i] = theta_out_.delta(i);
  for (int i=0; i < nphi_; i++) dphi_[i] = phi_out_.delta(i);

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
  // read vr
  status = H5LTread_dataset_double(file_id,"/vr",tmp);
  if (status < 0) if (verbose) std::cerr << "# Grid Err; can't find vr" << endl;
  for (int i=0; i < n_zones; i++) z[i].v[0] = tmp[i];
  // read vtheta
  status = H5LTread_dataset_double(file_id,"/vtheta",tmp);
  if (status < 0) if (verbose) std::cerr << "# Grid Err; can't find vtheta" << endl;
  for (int i=0; i < n_zones; i++) z[i].v[1] = tmp[i];
  // read vphi
  status = H5LTread_dataset_double(file_id,"/vphi",tmp);
  if (status < 0) if (verbose) std::cerr << "# Grid Err; can't find vphi" << endl;
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
  index_r_.resize(n_zones);
  index_theta_.resize(n_zones);
  index_phi_.resize(n_zones);

  //---------------------------------------------------
  // Calculate volume, indices, model properties
  //---------------------------------------------------
  double totmass = 0, totke = 0, totrad = 0;
  double *elem_mass = new double[n_elems];
  for (int l = 0;l < n_elems; l++) elem_mass[l] = 0;
  cnt = 0;
  for (int i=0;i<nr_;++i)
  {
    for (int j=0;j<ntheta_;++j)
    {
      for (int k=0;k<nphi_;++k)
      {
        double r0 = r_out_.left(i);
        double r1 = r_out_.right(i);
        double theta0 = theta_out_.left(j);
        double theta1 = theta_out_.right(j);
        double phi0 = phi_out_.left(k);
        double phi1 = phi_out_.right(k);
        vol_[cnt] = (4.*pc::pi/3.)*(r1*r1*r1 - r0*r0*r0) * ((theta1 - theta0)/pc::pi) * ((phi1 - phi0)/(2.*pc::pi));

        index_r_[cnt] = i;
        index_theta_[cnt] = j;
        index_phi_[cnt] = k;

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
    std::cout << "# gridtype: 3D spherical\n";
    std::cout << "# n_zones = " << n_zones << "\n";
    std::cout << "# (nr,ntheta,nphi) = (";
    std::cout << nr_ << ", " << ntheta_ << ", " << nphi_ << ")\n";
    // std::cout << "# (dx,dy,dz) = (";
    // std::cout << dx_ << ", " << dy_ << ", " << dz_ << ")\n";
    std::cout << "# (r_min,theta_min,phi_min) = (";
    std::cout << r_out_.minval() << ", " << theta_out_.minval() << ", " << phi_out_.minval() << ")\n";
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
void grid_3D_sphere::write_plotfile(int iw, double tt, int write_mass_fractions)
{
	// get file name
	char plotfilename[1000];
	sprintf(plotfilename,"plt_%05d.h5",iw);

	// open hdf5 file
	hid_t file_id = H5Fcreate(plotfilename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	// print out r array
	hsize_t  dims_r[1]={(hsize_t)nr_};
	float *rarr = new float[nr_];
	for (int i=0;i<nr_;i++) rarr[i] = r_out_.left(i);
	H5LTmake_dataset(file_id,"r",1,dims_r,H5T_NATIVE_FLOAT,rarr);
  delete [] rarr;

  // print out theta array
	hsize_t  dims_theta[1]={(hsize_t)ntheta_};
	float *thetaarr = new float[ntheta_];
	for (int i=0;i<ntheta_;i++) thetaarr[i] = theta_out_.left(i);
	H5LTmake_dataset(file_id,"theta",1,dims_theta,H5T_NATIVE_FLOAT,thetaarr);
  delete [] thetaarr;

  // print out phi array
	hsize_t  dims_phi[1]={(hsize_t)nphi_};
	float *phiarr = new float[nphi_];
	for (int i=0;i<nphi_;i++) phiarr[i] = phi_out_.left(i);
	H5LTmake_dataset(file_id,"phi",1,dims_phi,H5T_NATIVE_FLOAT,phiarr);
  delete [] phiarr;

  hsize_t  dims_g[3]={(hsize_t) nr_,(hsize_t) ntheta_,(hsize_t) nphi_};
  write_hdf5_plotfile_zones(file_id, dims_g, 3, tt);
  write_integrated_quantities(iw,tt);

  // Close the output file
  H5Fclose (file_id);

}



//************************************************************
// expand the grid
//************************************************************
void grid_3D_sphere::expand(double e)
{
  r_out_.scale(e);

  for (int i=0; i < nr_; i++) dr_[i] = dr_[i]*e;

  for (int i=0; i < n_zones; i++) vol_[i] = vol_[i]*e*e*e;
}



//------------------------------------------------------------
// Overly simple search to find zone
//------------------------------------------------------------
int grid_3D_sphere::get_zone(const double *x) const
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

  double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
  double theta = atan2(sqrt(x[0]*x[0] + x[1]*x[1]), x[2]);
  double phi = atan2(x[1], x[0]);
  if (phi < 0){
    phi += 2.*pc::pi;
  }

  if (r < r_out_.minval()) return -1;
  if (r > r_out_[nr_-1]) return -2;

  int i = r_out_.locate_within_bounds(r);
  int j = theta_out_.locate_within_bounds(theta);
  int k = phi_out_.locate_within_bounds(phi);

  int ind =  i*ntheta_*nphi_ + j*nphi_ + k;
  return ind;
}


//************************************************************
// Find distance to next zone along path
//************************************************************
int grid_3D_sphere::get_next_zone
(const double *x, const double *D, int i, double r_core, double *l) const
{
  double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
  // double theta = atan2(sqrt(x[0]*x[0] + x[1]*x[1]), x[2]);
  // double phi = atan2(x[1], x[0]);
  // if (phi < 0){
  //   phi += 2.*pc::pi;
  // }

  int ir = index_r_[i];
  int itheta = index_theta_[i];
  int iphi = index_phi_[i];

  double lr, ltheta, lphi;
  int d_ir = 0;
  int d_itheta = 0;
  int d_iphi = 0;

  // tiny offset so we don't land exactly on boundaries
  double tiny = 1e-10;

  //---------------------------------
  // distance to r interfaces
  // one must calculate the intersection of a ray and a sphere, since the inner and outer boundaries of constant radius are spheres
  //---------------------------------

  double r_bnd_out, r_bnd_in, lr_out, lr_in;
  double a = D[0]*D[0] + D[1]*D[1] + D[2]*D[2];
  double b = 2*(x[0]*D[0] + x[1]*D[1] + x[2]*D[2]);
  double c, det;

  // outer interface
  r_bnd_out = r_out_.right(ir);
  c = r*r - r_bnd_out*r_bnd_out;
  det = b*b - 4*a*c;
  if (det < 0) lr_out = std::numeric_limits<double>::infinity();
  else lr_out = (-1.*b + sqrt(det))/(2.*a);
  if (lr_out < 0) lr_out = std::numeric_limits<double>::infinity();

  // inner interface
  r_bnd_in = r_out_.left(ir);
  // check if inner boundary is from a core
  if (r_bnd_in <= r_core) r_bnd_in = r_core;
  c = r*r - r_bnd_in*r_bnd_in;
  det = b*b - 4*a*c;
  if (det < 0) lr_in = std::numeric_limits<double>::infinity();
  else lr_in = (-1.*b - sqrt(det))/(2.*a);
  if (lr_in < 0) lr_in = std::numeric_limits<double>::infinity();
  // if in innermost zone and there is no inner boundary (i.e. both r_min = 0 and r_core = 0), then we never hit the inner shell
  if ((ir == 0) && (r_bnd_in == 0)) lr_in = std::numeric_limits<double>::infinity();

  // if moving inward
  if(lr_in < lr_out){
    d_ir = -1;
    // if inner boundary is from a core
    if ((r_core > 0) && (r_bnd_in <= r_core)) lr = lr_in - tiny*dr_[ir];
    // if in innermost zone and there is an inner boundary
    else if ((ir == 0) && (r_bnd_in > 0)) lr = lr_in - tiny*dr_[ir];
    // if normal boundary crossing
    else lr = lr_in + tiny*dr_[ir+d_ir];
  }
  // if moving outward
  else{
    d_ir = 1;
    // if in outermost boundary
    if (ir == nr_-1) lr = lr_out - tiny*dr_[ir];
    // if normal boundary crossing
    else lr = lr_out + tiny*dr_[ir+d_ir];
  }

  //---------------------------------
  // distance to theta interfaces
  // one must calculate the intersection of a ray and a cone, since the inner and outer boundaries of constant theta are cones
  //---------------------------------

  double theta_bnd, theta_op;
  double t, t1, t2, tint1, tint2;
  bool is_correct_cone;
  double ltheta_out, ltheta_in;
  a=0, b=0, c=0, det=0;

  // outer interface
  theta_bnd = theta_out_.right(itheta);

  if (theta_bnd <= pc::pi/2.) {theta_op = theta_bnd;}
  else {theta_op = pc::pi - theta_bnd;}

  a = D[2]*D[2] - cos(theta_op)*cos(theta_op);
  b = 2.*( D[2]*x[2] - (D[0]*x[0] + D[1]*x[1] + D[2]*x[2])*cos(theta_op)*cos(theta_op) );
  c = x[2]*x[2] - (x[0]*x[0] + x[1]*x[1] + x[2]*x[2])*cos(theta_op)*cos(theta_op);
  det = b*b - 4.*a*c;

  if (theta_op == 0){
    ltheta_out = std::numeric_limits<double>::infinity();
  }
  else if (theta_op == pc::pi/2.){
    if (D[2] == 0) ltheta_out = std::numeric_limits<double>::infinity();
    else ltheta_out = -x[2]/D[2];
    if (ltheta_out < 0) ltheta_out = std::numeric_limits<double>::infinity();    
  }
  else if (a == 0){
    if (b == 0){
      ltheta_out = std::numeric_limits<double>::infinity();
    }
    else{
      t = -c/b;
      if (theta_bnd < pc::pi/2.) {is_correct_cone = (x[2] + t*D[2] > 0);}
      else {is_correct_cone = (x[2] + t*D[2] < 0);}
      if ((t>0) && is_correct_cone) ltheta_out = t;
      else ltheta_out = std::numeric_limits<double>::infinity();
    }
  }
  else{
    if (det <= 0) ltheta_out = std::numeric_limits<double>::infinity();
    else{
      t1 = (-1.*b - sqrt(det))/(2.*a);
      if (theta_bnd < pc::pi/2.) {is_correct_cone = (x[2] + t1*D[2] > 0);}
      else {is_correct_cone = (x[2] + t1*D[2] < 0);}
      if ((t1>0) && is_correct_cone) tint1 = t1;
      else tint1 = std::numeric_limits<double>::infinity();

      t2 = (-1.*b + sqrt(det))/(2.*a);
      if (theta_bnd < pc::pi/2.) {is_correct_cone = (x[2] + t2*D[2] > 0);}
      else {is_correct_cone = (x[2] + t2*D[2] < 0);}
      if ((t2>0) && is_correct_cone) tint2 = t2;
      else tint2 = std::numeric_limits<double>::infinity();

      ltheta_out = fmin(tint1,tint2);
    }
  }

  // inner interface
  theta_bnd = theta_out_.left(itheta);

  if (theta_bnd <= pc::pi/2.) {theta_op = theta_bnd;}
  else {theta_op = pc::pi - theta_bnd;}

  a = D[2]*D[2] - cos(theta_op)*cos(theta_op);
  b = 2.*( D[2]*x[2] - (D[0]*x[0] + D[1]*x[1] + D[2]*x[2])*cos(theta_op)*cos(theta_op) );
  c = x[2]*x[2] - (x[0]*x[0] + x[1]*x[1] + x[2]*x[2])*cos(theta_op)*cos(theta_op);
  det = b*b - 4.*a*c;

  if (theta_op == 0){
    ltheta_in = std::numeric_limits<double>::infinity();
  }
  else if (theta_op == pc::pi/2.){
    if (D[2] == 0) ltheta_in = std::numeric_limits<double>::infinity();
    else ltheta_in = -x[2]/D[2];
    if (ltheta_in < 0) ltheta_in = std::numeric_limits<double>::infinity();    
  }
  else if (a == 0){
    if (b == 0){
      ltheta_in = std::numeric_limits<double>::infinity();
    }
    else{
      t = -c/b;
      if (theta_bnd < pc::pi/2.) {is_correct_cone = (x[2] + t*D[2] > 0);}
      else {is_correct_cone = (x[2] + t*D[2] < 0);}
      if ((t>0) && is_correct_cone) ltheta_in = t;
      else ltheta_in = std::numeric_limits<double>::infinity();
    }
  }
  else{
    if (det <= 0) ltheta_in = std::numeric_limits<double>::infinity();
    else{
      t1 = (-1.*b - sqrt(det))/(2.*a);
      if (theta_bnd < pc::pi/2.) {is_correct_cone = (x[2] + t1*D[2] > 0);}
      else {is_correct_cone = (x[2] + t1*D[2] < 0);}
      if ((t1>0) && is_correct_cone) tint1 = t1;
      else tint1 = std::numeric_limits<double>::infinity();

      t2 = (-1.*b + sqrt(det))/(2.*a);
      if (theta_bnd < pc::pi/2.) {is_correct_cone = (x[2] + t2*D[2] > 0);}
      else {is_correct_cone = (x[2] + t2*D[2] < 0);}
      if ((t2>0) && is_correct_cone) tint2 = t2;
      else tint2 = std::numeric_limits<double>::infinity();

      ltheta_in = fmin(tint1,tint2);
    }
  }

  // if moving inward
  if(ltheta_in < ltheta_out){
    d_itheta = -1;
    double x_new[3] = {x[0] + ltheta_in*D[0], x[1] + ltheta_in*D[1], x[2] + ltheta_in*D[2]};
    double r_new = sqrt(x_new[0]*x_new[0] + x_new[1]*x_new[1] + x_new[2]*x_new[2]);
    int new_itheta = itheta + d_itheta;
    if (new_itheta == -1) new_itheta = 1;
    ltheta = ltheta_in + tiny*r_new*dtheta_[new_itheta];
  }
  // if moving outward
  else{
    d_itheta = 1;
    double x_new[3] = {x[0] + ltheta_out*D[0], x[1] + ltheta_out*D[1], x[2] + ltheta_out*D[2]};
    double r_new = sqrt(x_new[0]*x_new[0] + x_new[1]*x_new[1] + x_new[2]*x_new[2]);
    int new_itheta = itheta + d_itheta;
    if (new_itheta == ntheta_) new_itheta = ntheta_-2;
    ltheta = ltheta_out + tiny*r_new*dtheta_[new_itheta];
  }

  //---------------------------------
  // distance to phi interfaces
  // one must calculate the intersection of a ray and a plane, since the inner and outer boundaries of constant phi are planes
  //---------------------------------

  double phi_bnd;
  double n[3];
  a=0, b=0;
  double lphi_out, lphi_in;

  // outer interface
  phi_bnd = phi_out_.right(iphi);
  n[0] = -sin(phi_bnd);
  n[1] = cos(phi_bnd);
  n[2] = 0;
  a = -(x[0]*n[0] + x[1]*n[1] + x[2]*n[2]);
  b = D[0]*n[0] + D[1]*n[1] + D[2]*n[2];
  if (b == 0) lphi_out = std::numeric_limits<double>::infinity();
  else lphi_out = a/b;
  if (lphi_out < 0) lphi_out = std::numeric_limits<double>::infinity();

  if (lphi_out == 0){
    cout << x[0] << " " << x[1] << " " << x[2] << " " << D[0] << " " << D[1] << " " << D[2] << endl;
  }

  // inner interface
  phi_bnd = phi_out_.left(iphi);
  n[0] = -sin(phi_bnd);
  n[1] = cos(phi_bnd);
  n[2] = 0;
  a = -(x[0]*n[0] + x[1]*n[1] + x[2]*n[2]);
  b = D[0]*n[0] + D[1]*n[1] + D[2]*n[2];
  if (b == 0) lphi_in = std::numeric_limits<double>::infinity();
  else lphi_in = a/b;
  if (lphi_in < 0) lphi_in = std::numeric_limits<double>::infinity();

  if (lphi_in == 0){
    cout << x[0] << " " << x[1] << " " << x[2] << " " << D[0] << " " << D[1] << " " << D[2] << endl;
  }

  // if moving inward
  if(lphi_in < lphi_out){
    d_iphi = -1;
    double x_new[3] = {x[0] + lphi_in*D[0], x[1] + lphi_in*D[1], x[2] + lphi_in*D[2]};
    double r_new = sqrt(x_new[0]*x_new[0] + x_new[1]*x_new[1] + x_new[2]*x_new[2]);
    double theta_new = atan2(sqrt(x[0]*x[0] + x[1]*x[1]), x[2]);
    int new_iphi = iphi + d_iphi;
    if (new_iphi == -1) new_iphi = nphi_-1;
    lphi = lphi_in + tiny*r_new*sin(theta_new)*dphi_[new_iphi];
  }
  // if moving outward
  else{
    d_iphi = 1;
    double x_new[3] = {x[0] + lphi_out*D[0], x[1] + lphi_out*D[1], x[2] + lphi_out*D[2]};
    double r_new = sqrt(x_new[0]*x_new[0] + x_new[1]*x_new[1] + x_new[2]*x_new[2]);
    double theta_new = atan2(sqrt(x[0]*x[0] + x[1]*x[1]), x[2]);
    int new_iphi = iphi + d_iphi;
    if (new_iphi == nphi_) new_iphi = 0;
    lphi = lphi_out + tiny*r_new*sin(theta_new)*dphi_[new_iphi];
  }

  //---------------------------------
  // find shortest distance
  //---------------------------------

  int new_ir = ir;
  int new_itheta = itheta;
  int new_iphi = iphi;
  // if particle hits a r interface first
  if ((lr < ltheta) && (lr < lphi)){
    *l = lr;
    new_ir += d_ir;
    // if moving inward and inner boundary is from a core
    if ((d_ir == -1) && (r_core > 0) && (r_bnd_in <= r_core)) return -1;
    // if moving inward and in innermost zone and there is an inner boundary
    else if ((new_ir == -1) && (r_bnd_in > 0)) return -1;
    // if moving outward and in outermost zone
    else if (new_ir == nr_) return -2;
  }
  // if particles hits a theta interface first
  else if (ltheta < lphi){
    *l = ltheta;
    new_itheta += d_itheta;
    if (new_itheta == -1) new_itheta = 1;
    else if (new_itheta == ntheta_) new_itheta = ntheta_-2;
  }
  // if particles hits a phi interface first
  else{
    *l = lphi;
    new_iphi += d_iphi;
    if (new_iphi == -1) new_iphi = nphi_-1;
    else if (new_iphi == nphi_) new_iphi = 0;
  }

  int i_new = new_ir*(ntheta_*nphi_) + new_itheta*(nphi_) + new_iphi;
  return i_new;
}



//------------------------------------------------------------
// return volume of zone (precomputed)
//------------------------------------------------------------
double grid_3D_sphere::zone_volume(const int i) const
{
  return vol_[i];
}

//------------------------------------------------------------
// sample a random position within the spherical cell
//------------------------------------------------------------
void grid_3D_sphere::sample_in_zone
(const int i, const std::vector<double> ran,double r[3])
{
  int ir = index_r_[i];
  int itheta = index_theta_[i];
  int iphi = index_phi_[i];

  double r_inner = r_out_.left(ir);
  double r_outer = r_out_.right(ir);
  double r_samp = cbrt( r_inner*r_inner*r_inner + ran[0]*( r_outer*r_outer*r_outer - r_inner*r_inner*r_inner ) );
  double theta_samp = acos( cos(theta_out_.left(itheta)) - ran[1]*( cos(theta_out_.left(itheta)) - cos(theta_out_.right(itheta)) ) );
  double phi_samp = phi_out_.left(iphi) + ran[2]*dphi_[iphi];

  // cout << "theta_left, theta_right = " << theta_out_.left(itheta) << ", " << theta_out_.right(itheta) << endl;
  // cout << "cos(theta_left), cos(theta_right) = " << cos(theta_out_.left(itheta)) << ", " << cos(theta_out_.right(itheta)) << endl;
  // cout << "val = " << cos(theta_out_.left(itheta)) - ran[1]*( cos(theta_out_.left(itheta)) - cos(theta_out_.right(itheta)) ) << endl;
  // cout << "theta_samp = " << theta_samp << endl;

  r[0] = r_samp*sin(theta_samp)*cos(phi_samp);
  r[1] = r_samp*sin(theta_samp)*sin(phi_samp);
  r[2] = r_samp*cos(theta_samp);

  // cout << "x[0], x[1], x[2] = " << r[0] << ", " << r[1] << ", " << r[2] << endl;
}

//------------------------------------------------------------
// get the velocity vector
//------------------------------------------------------------
void grid_3D_sphere::get_velocity(int i, double x[3], double D[3], double v[3], double *dvds)
{

  // verbocity
#ifdef MPI_PARALLEL
  int my_rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
  const int verbose = (my_rank == 0);
#else
  const int verbose = 1;
#endif

  if (use_homologous_velocities_ == 1) {
    v[0] = x[0]/t_now;
    v[1] = x[1]/t_now;
    v[2] = x[2]/t_now;
    *dvds = 1.0/t_now;
  }
  // I have to add this functionality for non-homologous expansion
  // Note: v[3] is expecting output of the form (v_x, v_y, v_z), whereas z[i].v[3] is in the form (v_r, v_theta, v_phi)
  else{
    if (verbose) cerr << "# Error: grid_3D_sphere currently only works with homologous expansion. Exiting." << endl;
    exit(1);
  }
  // else {

  //   // trilinear interpolation
  //   // e.g., https://en.wikipedia.org/wiki/Trilinear_interpolation

  //   // indices along each dimensions for this zone
  //   int ic[3];
  //   ic[0] = index_x_[i];
  //   ic[1] = index_y_[i];
  //   ic[2] = index_z_[i];

  //   // put grid data into array form
  //   // this should probably be stored this way to start...
  //   int npts[3];
  //   npts[0] = nx_;
  //   npts[1] = ny_;
  //   npts[2] = nz_;
  //   double rmin[3];
  //   rmin[0] = x_out_.minval();
  //   rmin[1] = y_out_.minval();
  //   rmin[2] = z_out_.minval();
  //   double del[3];
  //   del[0] = dx_[ic[0]];
  //   del[1] = dy_[ic[1]];
  //   del[2] = dz_[ic[2]];

  //   // find lattice points to interpolate from
  //   int    i_lo[3], i_hi[3];
  //   double r_lo[3];
  //   for (int j=0;j<3;j++)
  //   {
  //     double cen = rmin[j] + (ic[j] + 0.5)*del[j];
  //     if (x[j] > cen)
  //     {
  //       r_lo[j] = cen;
  //       i_lo[j] = ic[j];
  //       i_hi[j] = ic[j] + 1;
  //       if (i_hi[j] >= npts[j]) i_hi[j] = i_lo[j];
  //     }
  //     else
  //     {
  //       r_lo[j] = rmin[j] + (ic[j] - 0.5)*del[j];
  //       i_lo[j] = ic[j] - 1;
  //       i_hi[j] = ic[j];
  //       if (i_lo[j] < 0) i_lo[j] = i_hi[j];
  //     }
  //   }

  //   // differences in each dimension
  //   double xd = (x[0] - r_lo[0])/dx_[ic[0]];
  //   double yd = (x[1] - r_lo[1])/dy_[ic[1]];
  //   double zd = (x[2] - r_lo[2])/dz_[ic[2]];

  //   // interpolate each of the 3 velocity components
  //   for (int j = 0;j < 3; j++)
  //   {
  //     // values at 8 lattice points
  //     double c000 = z[get_index(i_lo[0],i_lo[1],i_lo[2])].v[j];
  //     double c001 = z[get_index(i_lo[0],i_lo[1],i_hi[2])].v[j];
  //     double c010 = z[get_index(i_lo[0],i_hi[1],i_lo[2])].v[j];
  //     double c011 = z[get_index(i_lo[0],i_hi[1],i_hi[2])].v[j];
  //     double c100 = z[get_index(i_hi[0],i_lo[1],i_lo[2])].v[j];
  //     double c101 = z[get_index(i_hi[0],i_lo[1],i_hi[2])].v[j];
  //     double c110 = z[get_index(i_hi[0],i_hi[1],i_lo[2])].v[j];
  //     double c111 = z[get_index(i_hi[0],i_hi[1],i_hi[2])].v[j];

  //     // interpolation along x
  //     double c00 = c000*(1-xd) + c100*xd;
  //     double c01 = c001*(1-xd) + c101*xd;
  //     double c10 = c010*(1-xd) + c110*xd;
  //     double c11 = c011*(1-xd) + c111*xd;

  //     // interpolation along y
  //     double c0 = c00*(1-yd) + c10*yd;
  //     double c1 = c01*(1-yd) + c11*yd;

  //     // interpolation along z
  //     double c = c0*(1-zd) + c1*zd;
  //     v[j] = c;
  //   }

  //   // debug -- need to calculate this
  //   *dvds = 0;
  // }
  
}