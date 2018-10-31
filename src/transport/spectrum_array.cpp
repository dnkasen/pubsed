#include <stdio.h>
#include <iostream>
#include <math.h>
#include <string.h>
#include <vector>
#include "hdf5.h"
#include "hdf5_hl.h"
#include "sedona.h"

#include "spectrum_array.h"
#include "physical_constants.h"

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

using std::vector;
namespace pc = physical_constants;

//--------------------------------------------------------------
// Constructors
//--------------------------------------------------------------
spectrum_array::spectrum_array()
{
  strcpy(name,DEFAULT_NAME);
}

void spectrum_array::set_name(std::string n)
{
  strcpy(name,n.c_str());
}

//--------------------------------------------------------------
// Initialization and Allocation
//--------------------------------------------------------------
void spectrum_array::init(std::vector<double> t, std::vector<double> w,
		    int n_mu, int n_phi)
{
  // assign time grid
  double t_start = t[0];
  double t_stop  = t[1];
  double t_del   = t[2];
  this->time_grid.init(t_start,t_stop,t_del);
  int n_times  = this->time_grid.size();

  // assign frequency grid
  if ((w.size() != 4)&&(w.size() != 3)) {
    std::cout << "# improperly defined spectrum_nu_grid; need {nu_1, nu_2, dnu, (log?)}; exiting\n";
    exit(1); }

  double w_start = w[0];
  double w_stop  = w[1];
  double w_del   = w[2];

  // initialize the frequency grid
  if (w.size() == 3)
   this->wave_grid.init(w_start,w_stop,w_del);
  if (w.size() == 4)
  {
    if (w[3] == 1) wave_grid.log_init(w_start,w_stop,w_del);
    else wave_grid.init(w_start,w_stop,w_del);
  }
  int n_wave   = this->wave_grid.size();

  // asign mu grid
  this->mu_grid.init(-1,1,n_mu);

  // asign phi grid
  this->phi_grid.init(0,2*pc::pi,n_phi);

  // index parameters
  this->n_elements  = n_times*n_wave*n_mu*n_phi;
  this->a3 = n_phi;
  this->a2 = n_mu*a3;
  this->a1 = n_wave*a2;

  // allocate memory
  this->click.resize(n_elements);
  this->flux.resize(n_elements);

  // clear
  wipe();
}


//--------------------------------------------------------------
// Functional procedure: Wipe
//--------------------------------------------------------------
void spectrum_array::wipe()
{
  for (size_t i=0;i<click.size();i++)
  {
    flux[i]   = 0;
    click[i]  = 0;
  }
}



//--------------------------------------------------------------
// handles the indexing: should be called in this order
//    time, wavelength, mu, phi
//--------------------------------------------------------------
int spectrum_array::index(int t, int l, int m, int p)
{
  return t*a1 + l*a2 + m*a3 + p;
}


//--------------------------------------------------------------
// count a particle
////--------------------------------------------------------------
void spectrum_array::count(double t, double w, double E, double *D)
{
  double mu  = D[2];
  double phi = atan2(D[1],D[0]);

  // locate bin number in all dimensions
  int t_bin = time_grid.locate(t);
  int l_bin = wave_grid.locate(w);
  int m_bin = mu_grid.locate(mu);
  int p_bin = phi_grid.locate(phi);

  // keep all photons, even if off wavelength grid
  if (l_bin < 0) l_bin = 0;
  if (l_bin >= wave_grid.size()) l_bin = 0;

  // if off the grids, just return without counting
  if ((t_bin < 0)||(l_bin < 0)||(m_bin < 0)||(p_bin < 0)) return;
  if (t_bin >= time_grid.size()) return;
  if (m_bin >= mu_grid.size())   return;
  if (p_bin >= phi_grid.size())  return;

  // add to counters
  int ind      = index(t_bin,l_bin,m_bin,p_bin);

  flux[ind]  += E;
  click[ind] += 1;
}



//--------------------------------------------------------------
// print out
//--------------------------------------------------------------
void spectrum_array::print()
{
   // get file name
  char specfile[1000];
  sprintf(specfile,"%s.dat",name);

  FILE *out = fopen(specfile,"w");

  int n_times  = this->time_grid.size();
  int n_wave   = this->wave_grid.size();
  int n_mu     = this->mu_grid.size();
  int n_phi    = this->phi_grid.size();

  fprintf(out,"# %d %d %d %d\n",n_times,n_wave,n_mu,n_phi);

  double *darray = new double[n_elements];

  // unitize and printout
  for (int k=0;k<n_mu;k++)
    for (int m=0;m<n_phi;m++)
      for (int i=0;i<n_times;i++)
	     for (int j=0;j<n_wave;j++)
        {
	       int id = index(i,j,k,m);
	       if (n_times > 1)  fprintf(out,"%12.4e ",time_grid.center(i));;
	       if (n_wave > 1)   fprintf(out,"%12.4e ",wave_grid.center(j));
	       if (n_mu > 1)     fprintf(out,"%12.4f ",mu_grid.center(k));
	       if (n_phi> 1)     fprintf(out,"%12.4f ",phi_grid.center(m));

    	   double norm = 1.0/(n_mu*n_phi);
	       if (n_wave > 1)  norm *= wave_grid.delta(j);
	       if (n_times > 1) norm *= time_grid.delta(i);

         // normalize it
         darray[id] = flux[id]/norm;

    	   fprintf(out,"%12.5e %10d\n", darray[id],click[id]);
    	 }
  fclose(out);

  // write hdf5 spectrum file
  sprintf(specfile,"%s.h5",name);
  hid_t file_id = H5Fcreate( specfile, H5F_ACC_TRUNC, H5P_DEFAULT,  H5P_DEFAULT);
  const int RANK = 1;

  // write nu grid
  int n_nu = wave_grid.size();
  float* tmp_array = new float[n_nu];
  hsize_t  dims_nu[RANK]={(hsize_t)n_nu};
  for (int j=0;j<n_nu;j++) tmp_array[j] = wave_grid.center(j);
  H5LTmake_dataset(file_id,"nu",RANK,dims_nu,H5T_NATIVE_FLOAT,tmp_array);
  delete[] tmp_array;

  // write mu grid
  tmp_array = new float[n_mu];
  hsize_t  dims_mu[RANK]={(hsize_t)n_mu};
  for (int j=0;j<n_mu;j++) tmp_array[j] = mu_grid[j];
  H5LTmake_dataset(file_id,"mu",RANK,dims_mu,H5T_NATIVE_FLOAT,tmp_array);
  delete[] tmp_array;

  // write time grid
  int n_t = time_grid.size();
  tmp_array = new float[n_t];
  hsize_t dims_t[RANK]={(hsize_t)n_t};
  for (int j=0;j<n_t;j++) tmp_array[j] = time_grid.center(j);
  H5LTmake_dataset(file_id,"time",RANK,dims_t,H5T_NATIVE_FLOAT,tmp_array);
  delete[] tmp_array;

  // write fluxes array
  if (n_mu == 1)
  {
    const int RANKF = 2;
    hsize_t  dims_flux[RANKF]={(hsize_t)n_t,(hsize_t)n_nu};
    H5LTmake_dataset(file_id,"Lnu",RANKF,dims_flux,H5T_NATIVE_DOUBLE,darray);
  }
  else
  {
    const int RANKF = 3;
    hsize_t  dims_flux[RANKF]={(hsize_t)n_t,(hsize_t)n_nu,(hsize_t)n_mu};
    H5LTmake_dataset(file_id,"Lnu",RANKF,dims_flux,H5T_NATIVE_DOUBLE,darray);
  }
  delete[] darray;

  H5Fclose (file_id);

}


void  spectrum_array::rescale(double r)
{
  for (size_t i=0;i<flux.size();i++) flux[i] *= r;
}



//--------------------------------------------------------------
// MPI average the spectrum contents
//--------------------------------------------------------------
// only process 0 gets the reduced spectrum to print
void spectrum_array::MPI_average()
{
#ifdef MPI_PARALLEL

  int receiving_ID = 0;
  int mpi_procs, myID;

  // average the flux (receive goes out of scope after section)
  {
    vector<double> receive;
    receive.resize(n_elements);
    for (int i=0;i<n_elements;++i) receive[i] = 0.0;

//    MPI_Reduce(&flux.front(), &receive.front(), n_elements, MPI_DOUBLE, MPI_SUM, receiving_ID, MPI_COMM_WORLD);
    MPI_Allreduce(&flux.front(), &receive.front(), n_elements, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    flux.swap(receive);
  }

  // average clicks (receive goes out of scope after section)
  {
    vector<int> receive;
    receive.resize(n_elements);
    for (int i=0;i<n_elements;++i) receive[i] = 0.0;
    MPI_Allreduce(&click.front(), &receive.front(), n_elements, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    click.swap(receive);
  }

  // only have the receiving ID do the division
  MPI_Comm_size( MPI_COMM_WORLD, &mpi_procs );
  MPI_Comm_rank( MPI_COMM_WORLD, &myID      );
  //if(myID == receiving_ID){
 //   #pragma omp parallel for
    for (int i=0;i<n_elements;i++)
      flux[i]  /= mpi_procs;

#endif
}
