#include <mpi.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <string.h>
#include <vector>
#include "spectrum_array.h"
#include "physical_constants.h"

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
  for (int i=0;i<click.size();i++) 
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
  FILE *out = fopen(name,"w");

  int n_times  = this->time_grid.size();
  int n_wave   = this->wave_grid.size();
  int n_mu     = this->mu_grid.size();
  int n_phi    = this->phi_grid.size();

  fprintf(out,"# %d %d %d %d\n",n_times,n_wave,n_mu,n_phi);

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

	  double norm = n_mu*n_phi;
	  if (n_wave > 1)  norm *= wave_grid.delta(j);
	  if (n_times > 1) norm *= time_grid.delta(i);
	  
	  fprintf(out,"%12.5e %10d\n", flux[id]/norm,click[id]);
	}
  fclose(out);
}


void  spectrum_array::rescale(double r)
{
  for (int i=0;i<flux.size();i++) flux[i] *= r;
}



//--------------------------------------------------------------
// MPI average the spectrum contents
//--------------------------------------------------------------
// only process 0 gets the reduced spectrum to print
void spectrum_array::MPI_average()
{
  int receiving_ID = 0;
  int mpi_procs, myID;

  // average the flux (receive goes out of scope after section)
  {
    vector<double> receive;
    receive.resize(n_elements);
    MPI_Reduce(&flux.front(), &receive.front(), n_elements, MPI_DOUBLE, MPI_SUM, receiving_ID, MPI_COMM_WORLD);
    flux.swap(receive);
  }

  // average clicks (receive goes out of scope after section)
  {
    vector<int> receive;
    receive.resize(n_elements);
    MPI_Reduce(&click.front(), &receive.front(), n_elements, MPI_INT, MPI_SUM, receiving_ID, MPI_COMM_WORLD);
    click.swap(receive);
  }
  
  // only have the receiving ID do the division
  MPI_Comm_size( MPI_COMM_WORLD, &mpi_procs );
  MPI_Comm_rank( MPI_COMM_WORLD, &myID      );
  if(myID == receiving_ID){
    #pragma omp parallel for
    for (int i=0;i<n_elements;i++)
    {
      flux[i]  /= mpi_procs;
      click[i] /= mpi_procs;
    }
  }
}
