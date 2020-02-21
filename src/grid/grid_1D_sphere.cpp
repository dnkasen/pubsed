#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <cassert>

#include "grid_1D_sphere.h"
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
void grid_1D_sphere::read_model_file(ParameterReader* params)
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

  // determine the model file format
  std::string mod_extension(".mod");
  std::string hdf5_extension(".h5");

  std::size_t found = model_file.find(hdf5_extension);
  if (found!=std::string::npos)
  {
    std::cout << "# model file is an hdf5 file (.h5)" << endl;
    read_hdf5_file(model_file,params,verbose);
  }
  else
  {
    found = model_file.find(mod_extension);
    if (verbose)
    {
      if (found!=std::string::npos)
        std::cout << "# model file is ASCII format (.mod)" << endl;
      else
        cerr << "# Don't recognize model file extension, assuming ascii" << endl;
    }
    read_ascii_file(model_file,verbose);
  }
}

//------------------------------------------------------------
//------------------------------------------------------------
// Read model data from an hdf5 input file
//------------------------------------------------------------
//------------------------------------------------------------
void grid_1D_sphere::read_hdf5_file(std::string model_file, ParameterReader* params, int verbose)
{
  // open hdf5 file
  hid_t file_id = H5Fopen(model_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  herr_t status;

  // get time
  double tt[1];
  status = H5LTread_dataset_double(file_id,"/time",tt);
  if (status < 0) if (verbose) std::cerr << "# Grid Err; can't find time" << endl;
  t_now = tt[0];

  //get inner radius
  double rm[1];
  status = H5LTread_dataset_double(file_id,"/r_min",rm);
  if (status < 0) if (verbose) std::cerr << "# Grid Err; can't find r_min" << endl;
  r_out.min = rm[0];

  v_inner_ = 0.; // just like in .mod case

  // get grid size and dimensions
  hsize_t dims[2];
  status = H5LTget_dataset_info(file_id,"/comp",dims, NULL, NULL);
  if (status < 0) if (verbose) std::cerr << "# Grid Err; can't find comp" << endl;

  n_zones = dims[0];
  n_elems = dims[1];
  z.resize(n_zones);
  r_out.resize(n_zones);
  vol.resize(n_zones);

  int *etmp = new int[n_elems];
  status = H5LTread_dataset_int(file_id,"/Z",etmp);
  if (status < 0) if (verbose) std::cerr << "# Grid Err; can't find Z" << endl;
  for (int k=0;k<n_elems;k++) elems_Z.push_back(etmp[k]);
  status = H5LTread_dataset_int(file_id,"/A",etmp);
  if (status < 0) if (verbose) std::cerr << "# Grid Err; can't find A" << endl;
  for (int k=0;k<n_elems;k++) elems_A.push_back(etmp[k]);
  delete [] etmp;

  double *tmp = new double[n_zones];
  // read radii
  status = H5LTread_dataset_double(file_id,"/r_out",tmp);
  if (status < 0) if (verbose) std::cerr << "# Grid Err; can't find r_out" << endl;
  for (int i=0; i < n_zones; i++) r_out[i] = tmp[i];
  // read density
  status = H5LTread_dataset_double(file_id,"/rho",tmp);
  if (status < 0) if (verbose) std::cerr << "# Grid Err; can't find rho" << endl;
  for (int i=0; i < n_zones; i++) z[i].rho = tmp[i];
  // read temperature
  status = H5LTread_dataset_double(file_id,"/temp",tmp);
  if (status < 0) if (verbose) std::cerr << "# Grid Err; can't find temp" << endl;
  for (int i=0; i < n_zones; i++) z[i].T_gas = tmp[i];
  // read v
  status = H5LTread_dataset_double(file_id,"/v",tmp);
  if (status < 0) if (verbose) std::cerr << "# Grid Err; can't find v" << endl;
  for (int i=0; i < n_zones; i++) z[i].v[0] = tmp[i];
  // read erad
  status = H5LTread_dataset_double(file_id,"/erad",tmp);
  if (status < 0)
  {
    if (verbose) std::cout << "# Grid warning: Can't find erad. Using gas temperature and assuming blackbody radiation field." << endl;
    for (int i=0; i < n_zones; i++) z[i].e_rad = pc::a * pow(z[i].T_gas,4.);
  }
  else
  {
    for (int i=0; i < n_zones; i++) z[i].e_rad = tmp[i];
  }
  // read grey opacity if the user defines a zone-specific grey opacity
  int use_zone_specific_grey_opacity = params->getScalar<int>("opacity_zone_specific_grey_opacity");
  if(use_zone_specific_grey_opacity != 0){
    status = H5LTread_dataset_double(file_id,"/grey_opacity",tmp);
    if (status < 0)
    {
      if (verbose) std::cerr << "# Grid warning: Can't find grey_opacity. Setting zone-specific component of grey opacity to zero." << endl;
      for (int i=0; i < n_zones; i++) z[i].zone_specific_grey_opacity = 0;
    }
    else
    {
      for (int i=0; i < n_zones; i++) z[i].zone_specific_grey_opacity = tmp[i];
    }
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
    double norm = 0;
    for (int k=0; k < n_elems;  k++)
    {
      z[i].X_gas[k] = ctmp[cnt];
      norm += z[i].X_gas[k];
      cnt++;
    }

    // Make sure initial compositions are normalized, and compute mu
    double inverse_mu_sum = 0.;
    for (int k = 0; k < n_elems; k++)
    {
      z[i].X_gas[k] /= norm;
      inverse_mu_sum += z[i].X_gas[k]/elems_A[k];
    }
    z[i].mu_I = 1./inverse_mu_sum;
  }
  delete [] ctmp;

  // close HDF5 input file
  H5Fclose (file_id);

  for (int i=0; i < n_zones; i++)
  {
    // calculate shell volume
    double r0;
    if(i==0) r0 = r_out.min;
    else     r0 = r_out[i-1];
    vol[i] = 4.0*pc::pi/3.0*(r_out[i]*r_out[i]*r_out[i] - r0*r0*r0);
  }


  // print out properties of the model
  if (verbose)
  {
    cout << "#\n####### 1D SPHERE STANDARD MODEL ##########\n";
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
    cout << "##############################\n#" << endl;

  }

}

//------------------------------------------------------------
//------------------------------------------------------------
// Read model data from an ascii input file
//------------------------------------------------------------
//------------------------------------------------------------
void grid_1D_sphere::read_ascii_file(std::string model_file, int verbose)
{
  std::ifstream infile;
  infile.open(model_file.c_str());
  if(infile.fail())
  {
    if (verbose) cerr << "Err: can't read model file: " << model_file << endl;
    exit(4);
  }

  // geometry of model
  infile >> grid_type;
  if(grid_type != "1D_sphere")
  {
    if (verbose) cerr << "Err: grid_type param disagrees with the model file" << endl;
    exit(4);
  }
  if (verbose)
  {
    cout << "# model file = " << model_file << "\n";
    cout << "# Model is a 1D_sphere\n";
  }

  // type of system
  string system;
  infile >> system;

  // number of zones
  infile >> n_zones;
  z.resize(n_zones);
  r_out.resize(n_zones);
  vol.resize(n_zones);

  // read style of this model file
  int snr = 0;
  if (system == "SNR")
    snr = 1;
  else if (system == "standard")
    snr = 0;
  else
  {
    if (verbose) cerr << " Don't recognize model type " << system << "; Exiting" << endl;
    exit(1);
  }

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

  // read bulk grey opacity (note: this parameter is set in the param file, not in the ascii file)
  double bulk_grey_opacity = params->getScalar<double>("opacity_grey_opacity");

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
    }
    // read composition

    double norm = 0.;
    for (int k=0;k<n_elems;k++)
    {
      double x;
      infile >> x;
      z[i].X_gas.push_back(x);
      norm += x;
    }

    // Make sure initial compositions are normalized, and compute mu
    double inverse_mu_sum = 0.;
    for (int k = 0; k < n_elems; k++)
    {
      z[i].X_gas[k] /= norm;
      inverse_mu_sum += z[i].X_gas[k]/elems_A[k];
    }
    z[i].mu_I = 1./inverse_mu_sum;

    // assume LTE radiation field to start
    z[i].e_rad = pc::a*pow(z[i].T_gas,4);
    // DEBUG - this was left over from something...
    //z[i].e_rad = pc::a* pow(3.4e6,4);

    // assume zone-specific grey opacity is zero
    z[i].zone_specific_grey_opacity = 0;

    // set bulk grey opacity and total grey opacity
    z[i].bulk_grey_opacity = bulk_grey_opacity;
    z[i].total_grey_opacity = z[i].bulk_grey_opacity + z[i].zone_specific_grey_opacity;

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
    cout << "##############################\n#" << endl;

  }
}

//************************************************************
// expand the grid
//************************************************************
void grid_1D_sphere::expand(double e)
{
  for (int i=0;i<n_zones;i++)
    r_out[i] *= e;
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



void grid_1D_sphere::restartGrid(ParameterReader* params) {
  // verbocity
#ifdef MPI_PARALLEL
  int my_rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
  const int verbose = (my_rank == 0);
#else
  const int verbose = 1;
#endif
  string restart_file = params->getScalar<string>("run_restart_file");

  // geometry of model
  if(params->getScalar<string>("grid_type") != "grid_1D_sphere")
  {
    if (verbose) cerr << "Err: grid_type param disagrees with the model file" << endl;
    exit(4);
  }
  if (verbose) {
    cout << "# model file = " << restart_file << "\n";
    cout << "# Model is a 1D_sphere\n"; }

  readCheckpointGrid(restart_file);
  readCheckpointZones(restart_file);


  // print out properties of the model
  if (verbose)
  {
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
    cout << "##############################\n#" << endl;

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
  int ind = r_out.locate_within_bounds(r);
  return ind;
}

//************************************************************
// Overly simple search to find zone
//************************************************************
int grid_1D_sphere::get_next_zone(const double *x, const double *D, int i, double r_core, double *l) const
{
  double rsq   = (x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
  double xdotD = (D[0]*x[0] + D[1]*x[1] + D[2]*x[2]);

  // Calculate distance to the outer shell edge
  // using quadratic formula
  double r_o = r_out[i];
  double l_out = -1*xdotD + sqrt(xdotD*xdotD + r_o*r_o - rsq);


  double r_in;    // radius of the inner shell edge
  int ind_in;    // index of interior shell
  // get radius of inner shell edge
  if (i != 0)
  {
    r_in = r_out[i-1];
    ind_in = i-1;
  }
  // for innermost shell, use minimum r
  else
  {
    r_in = r_out.min;
    ind_in = -1;
  }

  // check for a core boundary
  if (r_core >= r_in)
  {
    r_in = r_core;
    ind_in = -1;
  }


  // find distance to inner shell
  double l_in;
  // if in innermost zone and there is no inner boundary,
  // (i.e., r_in = 0) then we never hit the inner shell
  // so set l_in = -1 as an indicator
  if ((i == 0)&&(r_in == 0)) l_in = -1;
  // otherwise calculate the distance
  else
  {
    double rad = xdotD*xdotD + r_in*r_in - rsq;
    if   (rad < 0)  l_in = -1;
    else l_in = -1*xdotD - sqrt(rad);
  }

  int ind;
   // offset so we don't land *exactly on a boundary
  double tiny_offset = 1 + 1e-6;

  // if l_out is shortest positive distance, set this as distance
  if ((l_out < l_in)||(l_in < 0))
  {
    ind = i + 1;

    // if in outermost zone, move to just inside the outer edge
    if (ind == n_zones)
    {
        ind = -2;
        *l = l_out/tiny_offset;
    }
    // else move a tiny bit past outer zone edge
    else
        *l = l_out*tiny_offset;
  }
  // otherwise set inward as distance to move
  else
  {
    ind = ind_in;

    // if at inner boundary, muve to just outside of it
    if (ind_in == -1)
        *l = l_in/tiny_offset;
    // otherwise move a tiny past inner zone edge
    else
        *l = l_in*tiny_offset;
  }
  return ind;
}




//************************************************************
// Write out the file
//************************************************************
void grid_1D_sphere::write_plotfile(int iw, double tt, int write_mass_fracs)
{
  // write ascii version, for convenience
  char zonefile[1000];
  sprintf(zonefile,"plt_%05d.dat",iw);

  FILE *outfile;
  outfile = fopen(zonefile,"w");

  fprintf(outfile,"# t = %8.4e ; rmin = %8.4e\n",tt, r_out.min);
  fprintf(outfile, "#  %-12.12s %-15.15s %-15.15s %-15.15s %-15.15s %-15.15s %-15.15s %-15.15s","r", "rho","v", "T_gas", "T_rad", "n_elec", "L_dep_nuc","L_emit_nuc");
  if (write_mass_fracs) // output mass fractions
  {
    for (int j =0; j < n_elems; j++)
    {
      char elem_id[10];
      sprintf(elem_id,"%d.%d",elems_Z[j],elems_A[j]);
      fprintf(outfile," %-15.15s",elem_id);
    }

  }
  fprintf(outfile,"\n");

  for (int i=0;i<n_zones;i++)
  {
    double rin = r_out.min;
    if (i > 0) rin = r_out[i-1];
    double T_rad = pow(z[i].e_rad/pc::a,0.25);

    fprintf(outfile, "%12.8e  %12.8e  %12.8e  %12.8e  %12.8e  %12.8e  %12.8e  %12.8e", r_out[i], z[i].rho, z[i].v[0], z[i].T_gas, T_rad, z[i].n_elec, z[i].L_radio_dep, z[i].L_radio_emit);
    if (write_mass_fracs) // output mass fractions
    {
      for (int j =0; j < n_elems; j++)
        fprintf(outfile,"  %12.8e", z[i].X_gas[j]);
    }
    fprintf(outfile,"\n");
  }

  fclose(outfile);

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
// sample a random position within the spherical shell weighted by volume
//************************************************************
void grid_1D_sphere::sample_in_zone(int i, std::vector<double> ran, double r[3])
{
  // inner radius of shell
  double r_0;
  if (i == 0) r_0 = r_out.min;
  else r_0 = r_out[i-1];

  // sample radial position in shell weighted by volume
  double r_samp = pow( r_0*r_0*r_0 + ran[0]*( r_out[i]*r_out[i]*r_out[i]-r_0*r_0*r_0 ), 1.0/3.0);

  // random spatial angles
  double mu  = 1 - 2.0*ran[1];
  double phi = 2.0*pc::pi*ran[2];
  double sin_theta = sqrt(1 - mu*mu);

  // set the real 3-d coordinates
  r[0] = r_samp*sin_theta*cos(phi);
  r[1] = r_samp*sin_theta*sin(phi);
  r[2] = r_samp*mu;
}



//************************************************************
// get the velocity vector
//************************************************************
void grid_1D_sphere::get_velocity(int i, double x[3], double D[3], double v[3], double *dvds)
{

  if (use_homologous_velocities_ == 1) {
    v[0] = x[0]/t_now;
    v[1] = x[1]/t_now;
    v[2] = x[2]/t_now;
    *dvds = 1.0/t_now;
  }
  else {
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

void grid_1D_sphere::writeCheckpointGrid(std::string fname) {
  if (my_rank == 0) {
    writeCheckpointGeneralGrid(fname);
    /* Specific to 1D */
    hsize_t single_val = 1;
    hsize_t zone_size = n_zones;

    createDataset(fname, "grid", "v_inner", 1, &single_val, H5T_NATIVE_DOUBLE);
    writeSimple(fname, "grid", "v_inner", &v_inner_, H5T_NATIVE_DOUBLE);

    r_out.writeCheckpoint(fname, "grid", "r_out");

    writeVector(fname, "grid", "vol", vol, H5T_NATIVE_DOUBLE);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

void grid_1D_sphere::readCheckpointGrid(std::string fname, bool test) {
  for (int rank = 0; rank < nproc; rank++) {
    if (my_rank == rank) {
      readCheckpointGeneralGrid(fname, test);
      /* Specific to 1D */
      readSimple(fname, "grid", "v_inner", &v_inner_new, H5T_NATIVE_DOUBLE);
      r_out_new.readCheckpoint(fname, "grid", "r_out");
      readVector(fname, "grid", "vol", vol_new, H5T_NATIVE_DOUBLE);
      if (not test) {
        v_inner_ = v_inner_new;
        r_out = r_out_new;
        vol = vol_new;
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
}
