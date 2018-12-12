#include "AtomicSpecies.h"
#include "physical_constants.h"
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_linalg.h>
#include <iostream>
#include "hdf5.h"
#include "hdf5_hl.h"

namespace pc = physical_constants;


//----------------------------------------------------------------
// Initialize an atom by reading in the data from the hdf5 file
// given by name "fname"
//----------------------------------------------------------------


int AtomicSpecies::initialize(std::string fname, int z, locate_array ng, int &levID)
{
  // set atomic number
  this->atomic_number = z;

  // copy over frequency grid
  nu_grid_.copy(ng);

  // open hdf5 file
  hid_t file_id = H5Fopen (fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  herr_t status;

  // get element group name
  char atomname[100];
  sprintf(atomname,"%d/",atomic_number);
  char dset[1000];

  // default, no fuzzlines
  fuzz_lines.n_lines = 0;

  // ----------------------------------------
  // read ions
  // ----------------------------------------
  status = H5LTget_attribute_int(file_id, atomname, "n_ions", &n_ions_);
  if (status != 0) return -1;
  ions_ = new AtomicIon[n_ions_];

  // read ionization potentials
  sprintf(dset,"%s%s",atomname,"ion_chi");
  double *ion_darr = new double[n_ions_];
  status = H5LTread_dataset_double(file_id,dset,ion_darr);
  if (status != 0) return -1;
  for (int i=0;i<n_ions_;++i) ions_[i].chi  = ion_darr[i];

  // read ground state levels
  sprintf(dset,"%s%s",atomname,"ion_ground");
  int *ion_iarr = new int[n_ions_];
  status = H5LTread_dataset_int(file_id,dset,ion_iarr);
  if (status != 0) return -1;
  for (int i=0;i<n_ions_;++i) ions_[i].ground  = ion_iarr[i];

  // initialize other quantities
  for (int i=0;i<n_ions_;++i)
  {
    ions_[i].part   = 0;
    ions_[i].frac   = 0;
    ions_[i].stage  = i;
  }

  // clean up
  delete[] ion_darr;
  delete[] ion_iarr;

  // ----------------------------------------
  // read levels
  // ----------------------------------------

  // get number of levels
  status = H5LTget_attribute_int(file_id, atomname, "n_levels",&n_levels_);
  if (status != 0) return -1;
  levels_ = new AtomicLevel[n_levels_];

  // make space to read data
  double *lev_darr = new double[n_levels_];
  int *lev_iarr = new int[n_levels_];

  // read level statistical weights
  sprintf(dset,"%s%s",atomname,"level_g");
  status = H5LTread_dataset_int(file_id,dset,lev_iarr);
  for (int i=0;i<n_levels_;++i) levels_[i].g = lev_iarr[i];

  // read level ionization stage
  sprintf(dset,"%s%s",atomname,"level_i");
  status = H5LTread_dataset_int(file_id,dset,lev_iarr);
  for (int i=0;i<n_levels_;++i) levels_[i].ion = lev_iarr[i];

  // read level excitation energy
  sprintf(dset,"%s%s",atomname,"level_E");
  status = H5LTread_dataset_double(file_id,dset,lev_darr);
  for (int i=0;i<n_levels_;++i) levels_[i].E = lev_darr[i];

  // set other parameters
  for (int i=0;i<n_levels_;++i)
  {
    levels_[i].n         = 0.0;
    levels_[i].E_ion     = ions_[levels_[i].ion].chi - levels_[i].E;

//    std::cout << levels_[i].ion << " " << i << " " << ions_[levels_[i].ion].chi << " " << levels[i].E << "\n";

    // set global ID (among all atoms being used)
    levels_[i].globalID  = levID;
    levID += 1;

    // set 0 g's to 1
    if (levels_[i].g == 0) levels_[i].g = 1;

    // find the level that this ionizes to (= -1 if none)
    levels_[i].ic = -1;
    for (int j=0;j<n_ions_;j++)
      if (ions_[j].stage == levels_[i].ion + 1)
      	levels_[i].ic  = ions_[j].ground;
  }

  // clean up
  delete[] lev_iarr;
  delete[] lev_darr;

  // ----------------------------------------
  // read lines
  // ----------------------------------------

  // get number of lines
  status = H5LTget_attribute_int(file_id, atomname, "n_lines",&n_lines_);
  if (status != 0) n_lines_ = 0;
  lines_ = new AtomicLine[n_lines_];

  if (n_lines_ > 0)
  {
    double *lin_darr = new double[n_lines_];
    int *lin_iarr    = new int[n_lines_];

    // read line lower levels
    sprintf(dset,"%s%s",atomname,"line_l");
    status = H5LTread_dataset_int(file_id,dset,lin_iarr);
    if (status != 0) return -1;
    for (int i=0;i<n_lines_;++i) lines_[i].ll = lin_iarr[i];

    // read line upper levels
    sprintf(dset,"%s%s",atomname,"line_u");
    status = H5LTread_dataset_int(file_id,dset,lin_iarr);
    for (int i=0;i<n_lines_;++i) lines_[i].lu = lin_iarr[i];

    // read line Einstein A
    sprintf(dset,"%s%s",atomname,"line_A");
    status = H5LTread_dataset_double(file_id,dset,lin_darr);
    for (int i=0;i<n_lines_;++i) lines_[i].A_ul = lin_darr[i];

    delete[] lin_darr;
    delete[] lin_iarr;
  }

  // set additional line properties
  for (int i=0;i<n_lines_;++i)
  {
    int ll = lines_[i].ll;
    int lu = lines_[i].lu;

    // get wavelength/frequency
    double delta_E = levels_[lu].E - levels_[ll].E;
    double nu      = delta_E*pc::ev_to_ergs/pc::h;
    double lam   = pc::c/nu*pc::cm_to_angs;
    lines_[i].nu  = nu;

    // set Einstein Coeficients
    int gl = levels_[ll].g;
    int gu = levels_[lu].g;
    double A = lines_[i].A_ul;
    lines_[i].B_ul = A*pc::c*pc::c/2.0/pc::h/nu/nu/nu;
    lines_[i].B_lu = lines_[i].B_ul*gu/gl;

    // set oscillator strength (see e.g., Rutten page 24)
    double lam_cm = lam*pc::angs_to_cm;
    lines_[i].f_lu = lam_cm*lam_cm*A*gu/gl/(8*pc::pi*pc::sigma_tot);

    // find index of bin in deal
    lines_[i].bin = nu_grid_.locate_within_bounds(nu);

    // default init tau and beta
    lines_[i].tau  = 0;
    lines_[i].beta = 1;
  }

  // ----------------------------------------
  // read photoionization cross-sections
  // if not available, just use hydrogenic approx
  // ----------------------------------------
  int npts     = 200;
  double fmax  = 10;
  for (int i=0;i<n_levels_;++i)
  {
     // set photoionization cross-section
     double E_ion = levels_[i].E_ion;
     double E_max = E_ion*fmax;
     double dE = (E_max - E_ion)/npts;
     levels_[i].s_photo.init(E_ion,E_max, dE);

     // effective excitation quantum number
     double E_ground = ions_[levels_[i].ion].chi;
     double n_eff = pow(1 - (E_ground - E_ion)/E_ground,-0.5);
     double s_fac = n_eff/(levels_[i].ion+1)/(levels_[i].ion + 1);

    double nu_t = E_ion*pc::ev_to_ergs/pc::h;
    int istart = 0; int istop = 0;
    for (int k=0;k<nu_grid_.size();k++)
    {
      if (nu_grid_.center(k) > nu_t)
        if (istart == 0) istart = k;

      if (nu_grid_.center(k) > 10*nu_t)
        if (istop == 0) istop = k;
    }
//    std::cout << "is " << i << " " << nu_t << " " << E_ion << "\n"; //istop - istart << " " << istop << "\n";

     //verner data
     // double V_Eth = 0.1360E+02;
     // double V_E0  = 0.4298E+00;
     // double V_s0  = 0.5475E+05;
     // double V_ya  = 0.3288E+02;
     // double V_P   = 0.2963E+01;
     // double V_yw  = 0.0000E+00;

     for (int j=0;j<npts;j++)
     {
        double E = levels_[i].s_photo.x[j];
        //double y = E/V_E0;
        //double Fr = ((y-1)*(y-1) + V_yw*V_yw)*
        double sigma = 6.3e-18*s_fac*pow(E/E_ion,-3.0); //2.5);
        levels_[i].s_photo.y[j] = sigma;
     }
  }

  H5Fclose(file_id);

  return 0;
}

//----------------------------------------------------------------
// setup atom for solving things in non-LTE
//----------------------------------------------------------------

int AtomicSpecies::set_use_nlte()
{
  use_nlte_ = true;

  //----------------------------------------------
  // allocate memory for arrays
  //----------------------------------------------
  rates_ = new double*[n_levels_];
  for (int i=0;i<n_levels_;++i) rates_[i] = new double[n_levels_];

  // matrix to solve
  M_nlte_ = gsl_matrix_calloc(n_levels_,n_levels_);
  gsl_matrix_set_zero(M_nlte_);

  // vector of level populations
  x_nlte_ = gsl_vector_calloc(n_levels_);
  gsl_vector_set_zero(x_nlte_);

  // right hand side vector
  b_nlte_ = gsl_vector_calloc(n_levels_);
  gsl_vector_set_zero(b_nlte_);

  // permuation vector, used internally for linear algebra solve
  p_nlte_ = gsl_permutation_alloc(n_levels_);
  gsl_permutation_init(p_nlte_);

  return 0;

}


int AtomicSpecies::read_fuzzfile(std::string fname)
{

  // open hdf5 file
  hid_t file_id = H5Fopen (fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  herr_t status;

  // get element group name
  char atomname[100];
  sprintf(atomname,"%d/",atomic_number);
  char dset[1000];

  int nl;
  status = H5LTget_attribute_int(file_id, atomname, "n_lines", &nl);
  if (status != 0) return -1;

  fuzz_lines.n_lines = nl;

  fuzz_lines.nu.resize(nl);
  fuzz_lines.gf.resize(nl);
  fuzz_lines.El.resize(nl);
  fuzz_lines.ion.resize(nl);
  fuzz_lines.bin.resize(nl);

  // read in arrays
  double *darr = new double[nl];
  int    *iarr = new int[nl];

  // read line frequency
  sprintf(dset,"%s%s",atomname,"nu");
  status = H5LTread_dataset_double(file_id,dset,darr);
  for (int i=0;i<nl;++i) fuzz_lines.nu[i] = darr[i];

  // read line gf
  sprintf(dset,"%s%s",atomname,"gf");
  status = H5LTread_dataset_double(file_id,dset,darr);
  for (int i=0;i<nl;++i) fuzz_lines.gf[i] = darr[i];

  // read line lower level excitation energy
  sprintf(dset,"%s%s",atomname,"El");
  status = H5LTread_dataset_double(file_id,dset,darr);
  for (int i=0;i<nl;++i) fuzz_lines.El[i] = darr[i];

  // read line ionization state
  sprintf(dset,"%s%s",atomname,"ion");
  status = H5LTread_dataset_int(file_id,dset,iarr);
  for (int i=0;i<nl;++i) fuzz_lines.ion[i] = iarr[i];

  // get frequency bin of line
  for (int i=0;i<nl;++i) {
    int ind = nu_grid_.locate(fuzz_lines.nu[i])-1;
    if (ind < 0) ind = 0;
    fuzz_lines.bin[i] = ind;
    //std::cout << fuzz_lines.nu[i] << " " << nu_grid_[fuzz_lines.bin[i]-1] << "\n";
  }

  delete[] darr;
  delete[] iarr;
  H5Fclose(file_id);

  return nl;

}
