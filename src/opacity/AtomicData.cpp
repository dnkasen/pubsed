#include "AtomicData.h"
#include "physical_constants.h"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "hdf5.h"
#include "hdf5_hl.h"

namespace pc = physical_constants;


//----------------------------------------------------------------
// Initialize an atom by reading in the data from the hdf5 file
// given by name "fname"
//----------------------------------------------------------------
AtomicData::AtomicData()
{
  // set all atoms to nothing
  for (int i=0;i<MAX_N_ATOMS;++i)
  {
    atomlist_[i].data_exists_ = false;
    atomlist_[i].n_ions_   = 0;
    atomlist_[i].n_levels_ = 0;
    atomlist_[i].n_lines_  = 0;
    atomlist_[i].fuzz_lines_.n_lines = 0;

    // default is to include all ion stages
    atomlist_[i].max_ion_stage_ = 9999;
  }
}

AtomicData::~AtomicData()
{
}


int AtomicData::initialize(std::string fname, locate_array ng)
{
  // copy over frequency grid
  nu_grid_.copy(ng);

  // copy over data file name
  atom_datafile_ = fname;

  // see if data file exists
  std::ifstream f(fname.c_str());
  if (!f.good())
  {
      std::cout << "# ERROR -- can't open atomic data file ";
      std::cout << fname << std::endl;
      return 1;
  }
  return 0;
}




//------------------------------------------------------------------------
// Print out general stats on all atomic data read and stored
//------------------------------------------------------------------------
void AtomicData::print()
{
  std::cout << "#-------------------------------------------------\n";
  std::cout << "# atomic data from: " << atom_datafile_ << "\n";
  std::cout << "#--------------------------------------------------\n";
  std::cout << "#  Z    n_ions  n_levels  n_lines  n_fuzz_lines\n";
  std::cout << "#-------------------------------------------------\n";
  for (int i=0;i<MAX_N_ATOMS;++i)
  {
    IndividualAtomData *atom = &(atomlist_[i]);
    if (atom->data_exists_ == false) continue;
    printf("# %3d   ",i);
    int nfuzz = atom->get_n_fuzz_lines();
    printf(" %4d %8d  %8d  %8d",atom->n_ions_,atom->n_levels_,atom->n_lines_,nfuzz);
    printf("\n");
  }
  std::cout << "#-------------------------------------------------\n";
}

//------------------------------------------------------------------------
// Print detailed stats for a certain element
//------------------------------------------------------------------------
void AtomicData::print_detailed(int z)
{

  IndividualAtomData *atom = &(atomlist_[z]);
  if (atom->data_exists_ == false)
  {
      std::cout << "# Can't print data for element " << z;
      std::cout << " ; doesn't exist\n";
      return;
  }

  std::cout << "#-------------------------------------------------\n";
  std::cout << "# atomic data from: " << atom_datafile_ << "\n";
  std::cout << "#--------------------------------------------------\n";
  std::cout << "#  Z    n_ions  n_levels  n_lines  n_fuzz_lines\n";
  std::cout << "#-------------------------------------------------\n";
  printf("# %3d   ",z);
  int nfuzz = atom->get_n_fuzz_lines();
  printf(" %4d %8d  %8d  %8d",atom->n_ions_,atom->n_levels_,atom->n_lines_,nfuzz);
  printf("\n");

 for (int i=0;i<atom->n_ions_;++i)
 {
    std::cout << "#-------------------------------------------------" << std::endl;
    std::cout << "ion stage       = " << atom->ions_[i].stage << std::endl;
    std::cout << "ionization chi  = " << atom->ions_[i].chi << " eV" << std::endl;
    std::cout << "ground level id = " << atom->ions_[i].ground << std::endl;
    std::cout << "#------------------------------------------------" << std::endl;;
    std::cout << std::endl;
  }
}

//------------------------------------------------------------------------
// Read all atomic data for species with atomic number Z
// Only include up to ionization stage max_ion
// (note ion = 0 is neutral)
//------------------------------------------------------------------------
int AtomicData::read_atomic_data(int z, int max_ion)
{
  // return if this is an invalid atom
  if ((z < 1)||(z >= MAX_N_ATOMS))
  {
    std::cout << "ERROR: Atomic number " << z << " is not allowed!" << std::endl;
    return 1;
  }

  atomlist_[z].max_ion_stage_ = max_ion;
  read_atomic_data(z);
}

//------------------------------------------------------------------------
// Read all atomic data for species with atomic number Z
//------------------------------------------------------------------------
int AtomicData::read_atomic_data(int z)
{
  // return if this is an invalid atom
  if ((z < 1)||(z >= MAX_N_ATOMS))
  {
    std::cout << "ERROR: Atomic number " << z << " is not allowed!" << std::endl;
    return 1;
  }

  // return if this atom has already been read
  if (atomlist_[z].data_exists_)
    return 0;

  // open hdf5 file
  herr_t status;
  status = H5Eset_auto1(NULL, NULL);
  hid_t file_id = H5Fopen (atom_datafile_.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  // get element group name
  char atomname[100];
  sprintf(atomname,"%d/",z);
  char dset[1000];

  IndividualAtomData *atom = &(atomlist_[z]);

  // check if this species is in the atom_datafile_
  // if not -- just fill up with single level
  status = H5Gget_objinfo (file_id, atomname, 0, NULL);
  if (status != 0)
  {
      std::cout << "# ERROR: No atomic data for species ";
      std::cout << z << " ; Filling with empty data\n";
      atom->n_ions_   = 1;
      atom->n_levels_ = 1;
      atom->n_lines_  = 0;
      atom->ions_ = new AtomicIon[1];
      atom->ions_[0].chi    = 999999;
      atom->ions_[0].ground = 0;
      atom->ions_[0].stage  = 0;
      atom->levels_ = new AtomicLevel[1];
      atom->levels_[0].E  = 0;
      atom->levels_[0].g  = 1;
      atom->levels_[0].ic = -1;
      atom->levels_[0].ion = 0;
      return 1;
  }

  // otherwise data exists, we'll read it in
  atomlist_[z].data_exists_ = true;

  //----------------------------------------
  // read ion data
  // ----------------------------------------
  int n_tot_ions;
  status = H5LTget_attribute_int(file_id, atomname, "n_ions", &n_tot_ions);
  if (status != 0) return -1;
  atom->n_ions_ = n_tot_ions;
  int add_last_stage = 0;
  if (atom->n_ions_  > atom->max_ion_stage_)
  {
    atom->n_ions_  = atom->max_ion_stage_;
    atom->ions_ = new AtomicIon[atom->n_ions_ + 1];
    add_last_stage = 1;
  }
  else
    atom->ions_ = new AtomicIon[atom->n_ions_];

  // allocate space to read in data
  int *ion_iarr = new int[n_tot_ions];
  double *ion_darr = new double[n_tot_ions];

  // read ionization potentials
  sprintf(dset,"%s%s",atomname,"ion_chi");
  status = H5LTread_dataset_double(file_id,dset,ion_darr);
  if (status != 0) return -1;
  for (int i=0;i<atom->n_ions_;++i) atom->ions_[i].chi  = ion_darr[i];

  // read ground state levels
  sprintf(dset,"%s%s",atomname,"ion_ground");
  status = H5LTread_dataset_int(file_id,dset,ion_iarr);
  if (status != 0) return -1;
  for (int i=0;i<atom->n_ions_;++i) atom->ions_[i].ground  = ion_iarr[i];

  for (int i=0;i<atom->n_ions_;++i)
    atom->ions_[i].stage  = i;


  // clean up
  delete[] ion_darr;
  delete[] ion_iarr;

  // ----------------------------------------
  // read levels
  // ----------------------------------------
  // get number of levels
  int tot_n_levels;
  status = H5LTget_attribute_int(file_id, atomname, "n_levels",&tot_n_levels);
  if (status != 0) return -1;

  // make space to read data
  double *lev_darr = new double[tot_n_levels];
  int *lev_iarr = new int[tot_n_levels];

  // read level ionization stages
  sprintf(dset,"%s%s",atomname,"level_i");
  status = H5LTread_dataset_int(file_id,dset,lev_iarr);

  // count number of levels to use
  int i;
  for (i=0;i<tot_n_levels;++i)
    if (lev_iarr[i] > atom->n_ions_-1) break;
  atom->n_levels_ = i;

  if (add_last_stage)
    atom->levels_ = new AtomicLevel[atom->n_levels_ + 1];
  else
    atom->levels_ = new AtomicLevel[atom->n_levels_];

  for (int i=0;i<atom->n_levels_;++i)
    atom->levels_[i].ion = lev_iarr[i];

  // read level statistical weights
  sprintf(dset,"%s%s",atomname,"level_g");
  status = H5LTread_dataset_int(file_id,dset,lev_iarr);
  for (int i=0;i<atom->n_levels_;++i) atom->levels_[i].g = lev_iarr[i];

  // read level excitation energy
  sprintf(dset,"%s%s",atomname,"level_E");
  status = H5LTread_dataset_double(file_id,dset,lev_darr);
  for (int i=0;i<atom->n_levels_;++i) atom->levels_[i].E = lev_darr[i];

  // set other parameters
  for (int i=0;i<atom->n_levels_;++i)
  {
    atom->levels_[i].E_ion = atom->ions_[atom->levels_[i].ion].chi - atom->levels_[i].E;

    // set 0 statistical weights to 1
    if (atom->levels_[i].g == 0) atom->levels_[i].g = 1;

    // find the level that this ionizes to (= -1 if none)
    atom->levels_[i].ic = -1;
    for (int j=0;j<atom->n_ions_;j++)
      if (atom->ions_[j].stage == atom->levels_[i].ion + 1)
      	atom->levels_[i].ic  = atom->ions_[j].ground;
  }

  // clean up
  delete[] lev_iarr;
  delete[] lev_darr;

  // ----------------------------------------
  // read lines
  // ----------------------------------------

  // get number of lines
  int n_tot_lines;
  status = H5LTget_attribute_int(file_id, atomname, "n_lines",&n_tot_lines);
  if (status != 0) n_tot_lines = 0;

  if (n_tot_lines > 0)
  {
    double *lin_darr = new double[n_tot_lines];
    int *lin_iarr    = new int[n_tot_lines];

    // read line lower level indices
    sprintf(dset,"%s%s",atomname,"line_l");
    status = H5LTread_dataset_int(file_id,dset,lin_iarr);
    if (status != 0) return -1;

    // count number of lines to use
    int i;
    for (i=0;i<n_tot_lines;++i)
      if (lin_iarr[i] >= atom->n_levels_) break;
    atom->n_lines_ = i;
    atom->lines_ = new AtomicLine[atom->n_lines_];

    // store lower level indices
    for (int i=0;i<atom->n_lines_;++i)
      atom->lines_[i].ll = lin_iarr[i];

    // read line upper level indices
    sprintf(dset,"%s%s",atomname,"line_u");
    status = H5LTread_dataset_int(file_id,dset,lin_iarr);
    for (int i=0;i<atom->n_lines_;++i) atom->lines_[i].lu = lin_iarr[i];

    // read line Einstein A
    sprintf(dset,"%s%s",atomname,"line_A");
    status = H5LTread_dataset_double(file_id,dset,lin_darr);
    for (int i=0;i<atom->n_lines_;++i) atom->lines_[i].A_ul = lin_darr[i];

    delete[] lin_darr;
    delete[] lin_iarr;
  }

  // set additional line properties
  for (int i=0;i<atom->n_lines_;++i)
  {
    int ll = atom->lines_[i].ll;
    int lu = atom->lines_[i].lu;

    //std::cout << i << " " << atom->lines_[i].ll << " " << atom->lines_[i].lu << "\n";

    // get wavelength/frequency
    double delta_E = atom->levels_[lu].E - atom->levels_[ll].E;
    double nu      = delta_E*pc::ev_to_ergs/pc::h;
    double lam   = pc::c/nu*pc::cm_to_angs;
    atom->lines_[i].nu  = nu;

    // set Einstein Coeficients
    int gl = atom->levels_[ll].g;
    int gu = atom->levels_[lu].g;
    double A = atom->lines_[i].A_ul;
    atom->lines_[i].B_ul = A*pc::c*pc::c/2.0/pc::h/nu/nu/nu;
    atom->lines_[i].B_lu = atom->lines_[i].B_ul*gu/gl;

    // set oscillator strength (see e.g., Rutten page 24)
    double lam_cm = lam*pc::angs_to_cm;
    atom->lines_[i].f_lu = lam_cm*lam_cm*A*gu/gl/(8*pc::pi*pc::sigma_tot);

    // find index of bin in deal
    atom->lines_[i].bin = nu_grid_.locate_within_bounds(nu);
  }

  if (add_last_stage)
  {
    int i = atom->n_ions_;
    int l = atom->n_levels_;
    atom->ions_[i].ground = l;
    atom->ions_[i].stage = i;
    atom->ions_[i].chi = 99999;
    atom->n_ions_ += 1;

    atom->levels_[l].ic = -1;
    atom->levels_[l].E   = 0;
    atom->levels_[l].g   = 1;
    atom->levels_[l].E_ion = 99999;
    atom->levels_[l].ion   = atom->n_ions_-1;
    atom->n_levels_     += 1;

    // find the level that this ionizes to (= -1 if none)
    for (int i=0;i<atom->n_levels_;++i)
    {
      atom->levels_[i].ic = -1;
      for (int j=0;j<atom->n_ions_;j++)
        if (atom->ions_[j].stage == atom->levels_[i].ion + 1)
          atom->levels_[i].ic  = atom->ions_[j].ground;
    }
  }

  // ----------------------------------------
  // read photoionization cross-sections
  // if not available, just use hydrogenic approx
  // ----------------------------------------
  int npts     = 100;
  double fmax  = 10;
  for (int i=0;i<atom->n_levels_;++i)
  {
     // set photoionization cross-section
     double E_ion = atom->levels_[i].E_ion;
     double E_max = E_ion*fmax;
     double dE = (E_max - E_ion)/npts;
     atom->levels_[i].s_photo.init(E_ion,E_max, dE);

     // effective excitation quantum number
     double E_ground = atom->ions_[atom->levels_[i].ion].chi;
     double n_eff = pow(1 - (E_ground - E_ion)/E_ground,-0.5);
     double s_fac = n_eff/(atom->levels_[i].ion+1)/(atom->levels_[i].ion + 1);

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
        double E = atom->levels_[i].s_photo.x[j];
        //double y = E/V_E0;
        //double Fr = ((y-1)*(y-1) + V_yw*V_yw)*
        double sigma = 6.3e-18*s_fac*pow(E/E_ion,-3.0); //2.5);
        atom->levels_[i].s_photo.y[j] = sigma;
     }
  }

   H5Fclose(file_id);
   return 0;
}


//------------------------------------------------------------------------
// Read a file constisting of (aka "fuzz") lines
// Do so for any atom that has data already defined
//------------------------------------------------------------------------
int AtomicData::read_fuzzfile_data(std::string fname)
{
  int n_lines = 0;
  for (int i=0;i<MAX_N_ATOMS;++i)
  {
    if (atomlist_[i].data_exists_)
      n_lines += read_fuzzfile_data_for_atom(fname,i);
  }
  return n_lines;
}


//------------------------------------------------------------------------
// Read a file constisting of (aka "fuzz") lines
// for a signel atom
//------------------------------------------------------------------------
int AtomicData::read_fuzzfile_data_for_atom(std::string fname, int Z)
{

  // open hdf5 file
  hid_t file_id = H5Fopen (fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  herr_t status;

  // get element group name
  char atomname[100];
  sprintf(atomname,"%d/",Z);
  char dset[1000];

  int nl;
  status = H5LTget_attribute_int(file_id, atomname, "n_lines", &nl);
  if (status != 0) return -1;

  atomlist_[Z].fuzz_lines_.n_lines = nl;
  atomlist_[Z].fuzz_lines_.nu.resize(nl);
  atomlist_[Z].fuzz_lines_.gf.resize(nl);
  atomlist_[Z].fuzz_lines_.El.resize(nl);
  atomlist_[Z].fuzz_lines_.ion.resize(nl);
  atomlist_[Z].fuzz_lines_.bin.resize(nl);

  // read in arrays
  double *darr = new double[nl];
  int    *iarr = new int[nl];

  // read line frequency
  sprintf(dset,"%s%s",atomname,"nu");
  status = H5LTread_dataset_double(file_id,dset,darr);
  for (int i=0;i<nl;++i) atomlist_[Z].fuzz_lines_.nu[i] = darr[i];

  // read line gf
  sprintf(dset,"%s%s",atomname,"gf");
  status = H5LTread_dataset_double(file_id,dset,darr);
  for (int i=0;i<nl;++i) atomlist_[Z].fuzz_lines_.gf[i] = darr[i];

  // read line lower level excitation energy
 sprintf(dset,"%s%s",atomname,"El");
 status = H5LTread_dataset_double(file_id,dset,darr);
 for (int i=0;i<nl;++i) atomlist_[Z].fuzz_lines_.El[i] = darr[i];

 // read line ionization state
 sprintf(dset,"%s%s",atomname,"ion");
 status = H5LTread_dataset_int(file_id,dset,iarr);
 for (int i=0;i<nl;++i) atomlist_[Z].fuzz_lines_.ion[i] = iarr[i];

 // get frequency bin of line
 for (int i=0;i<nl;++i)
 {
   int ind = nu_grid_.locate(atomlist_[Z].fuzz_lines_.nu[i])-1;
   if (ind < 0) ind = 0;
   atomlist_[Z].fuzz_lines_.bin[i] = ind;
 }

  delete[] darr;
  delete[] iarr;
  H5Fclose(file_id);

  return nl;
}
