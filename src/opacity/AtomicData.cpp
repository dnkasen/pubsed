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
    atomlist_[i].max_n_levels_  = 9999999;
  }
}

AtomicData::~AtomicData()
{
}

//------------------------------------------------------------------------
// Initialization
//------------------------------------------------------------------------

int AtomicData::initialize(std::string fname, locate_array ng)
{
  // copy over frequency grid
  nu_grid_.copy(ng);

  // copy over data file name
  atom_datafile_ = fname;

  // open hdf5 file to check file is OK to read
  herr_t status;
  status = H5Eset_auto1(NULL, NULL);
  hid_t file_id = H5Fopen (atom_datafile_.c_str(),  H5F_ACC_RDONLY, H5P_DEFAULT);

  // negative file id means open failed
  if (file_id < 0)
  {
    std::cerr << "# Error: Can't open atomic datafile: " << std::endl;
    std::cerr << "   \"" << atom_datafile_ <<"\"; " << std::endl;
    return 1;
  }

  // read file format version
  int version;
  status = H5LTread_dataset_int(file_id,"file_version"
    ,&version);
  // if version is not in file, assume original version
  if (status != 0)
    datafile_version_ = 0;
  else
    datafile_version_ = version;

  return 0;
}

//------------------------------------------------------------------------
// Return photo-ionization cross-section of level i
// at the center of frequency bin inu
// If the cross-section array does not exist, will calculate
// hydrogenic on the fly
//------------------------------------------------------------------------
double IndividualAtomData::get_lev_photo_cs(int i, int inu)
{

  // get index of this cross-section object
  int ics = levels_[i].cs;

  // if index ics < 0 the cross-section array doesn't exist
  // calculate hydrogenic cross-section on fly
  if (ics < 0)
  {
    double nu = nu_grid_->center(inu);
    double nu_t = levels_[i].nu_t;
    if (nu < nu_t) return 0;
    double fac = levels_[i].nu_t/nu;
    return levels_[i].sigma0*fac*fac*fac;
  }

  // if frequency bin is less than threshorld, cs = 0
  if (inu < photo_cs_[ics].i_start) return 0;

  // shift index (since cross-section arrays are only
  // stored for frequencies above threshold, which is at index i_start )
  int ishift = inu - photo_cs_[ics].i_start;
  return photo_cs_[ics].s[ishift];
}


//-----------------------------------------------
// Return the cross-section at energy E for ion state i
// E is in units of eV
// cross-section returned is in units of cm^2
// Returns 0 if no data for cross-section, or is invalid index
//-----------------------------------------------
double IndividualAtomData::get_nonthermal_ion_cross_section(int i, double E)
{
  if (i >= n_ions_) return 0;

  double cs = 0.0;
  for (int s=0;s<nt_ion_cs_[i].n_shells;++s)
  {
    double chi = nt_ion_cs_[i].chi[s];
    if (E < chi) continue;
    double u = E/chi;
    double v = nt_ion_cs_[i].A[s]*(1 - 1.0/u);
    v += nt_ion_cs_[i].B[s]*(1 - 1.0/u)*(1-1.0/u);
    v += nt_ion_cs_[i].C[s]*log(u);
    v += nt_ion_cs_[i].D[s]*log(u)/u;

    v = 1e-14*v/u/chi/chi;
    cs += v;
  }

  // if no-data available, use an approximation
  if (cs == 0)
  {
    double chi = get_ion_chi(i);
    double u = E/chi;
    double C = 2;
    cs = 1e-14/u/chi/chi*C*log(u);
  }

  return cs;
}


//-----------------------------------------------------
// Return the non-thermal bound-bound cross-section
// for line transition i
// for an electron of E (in units of eV)
// cross-section returned is in units of cm^2
// Uses Van-Regermorter approximation equation
// (see appendix of Botyanszki and Kasen 2018)
//---------------------------------------------------
double IndividualAtomData::get_nonthermal_bb_cross_section(int i, double E)
{
  double f_lu  = get_line_f(i);
  // energy difference in rydbergs
  double dE    = get_line_dE_eV(i)/pc::rydberg_ev;

  double effective_f_lu = f_lu;
  if (f_lu < 0.01) effective_f_lu = 0.01;

  // factor of 8 pi/sqrt(3)*pi*a_0**2 (a_0 = bohr radius)
  double fac = 4.0103404580362e-15;
  double gaunt = 0.2;
  double sigma = fac*(pc::rydberg_ev/E)/dE*effective_f_lu*gaunt;

  return sigma;
}


//-----------------------------------------------------
// get the collisional bound-bound rates for line i
// given a passed temperature T
// gives the excitation rate:  C_up
// and the de-excitaiton rate: C_down
// If no tabulated data exists, uses an analytic approx
//
// units here are cm^3/s -- must multiply by electron density
// to get rate per ion
//
// Rutten section 3.2.5 (page 52) points out that these van Regemorter
// rates are only valid for permitted dipole transitions, with f_lu in the
// range 10^-3 to 1. For forbidden lines with smaller f, he says the
// collisional transition rates "don't drop much below the values typical
// of permitted lines." That's not a very precise statement, but we can
// mock it up by not letting the f_lu factor drop below 10^-2
//------------------------------------------------------
void IndividualAtomData::get_collisional_bb_rates
(int i, double T, double &C_up, double &C_down)
{

  int ll       = get_line_l(i);
  int lu       = get_line_u(i);
  int gl       = get_lev_g(ll);
  int gu       = get_lev_g(lu);
  double El    = get_lev_E(ll);
  double Eu    = get_lev_E(lu);
  double dE    = (Eu - El)*pc::ev_to_ergs;
  double zeta  = dE/pc::k/T; // note dE is in ergs

  // get the omega
  double omega;
  // analytic approximation
  if (lines_[i].col_rate == NULL)
  {
    double f_lu  = get_line_f(i);
    int ion      = get_lev_ion(ll);

    double effective_f_lu = 0.;
    if (f_lu < 0.01) effective_f_lu = 0.01;
    else effective_f_lu = f_lu;

    if (ion == 0)
      omega =  2.16*pow(zeta,-1.68)*pow(T,-1.5)* effective_f_lu;
    else
      omega =  3.9*pow(zeta,-1.)*pow(T,-1.5)* effective_f_lu;
  }
  // tabulated data
  else
  {
    AtomicCollisionalBB_Rate  *col = lines_[i].col_rate;

    // locate this temperature
    if (T <= col->T[0])
      omega = col->O[0];
    else if (T >= col->T.back())
      omega = col->O.back();
    else
    {
      int ind = upper_bound(col->T.begin(), col->T.end(), T) - col->T.begin() - 1;
      if (ind < 0) ind = 0;

      int i1,i2;
      if (ind < col->T.size()-1)
      {
        i1    = ind;
        i2    = ind + 1;
      }
      else
      {
        i2    = ind;
        i1    = ind - 1;
      }
      // linear interpolation
      double slope = (col->O[i2]-col->O[i1])/(col->T[i2]-col->T[i1]);
      omega = col->O[ind] + slope*(T - col->T[i1]);
    }
  }
  
  C_up   = omega*exp(-zeta);
  // be careful about possible overflow
  if (zeta > 700) C_up = 0;
  C_down = omega*gl/gu;
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

  std::cout << "#---------- lines ---------------------------------" << std::endl;;
  for (int i=0;i<atom->n_lines_;++i)
  {
    std::cout << atom->lines_[i].nu << "\t" << atom->lines_[i].ll << "\t" << atom->lines_[i].lu << "\t"
      << atom->lines_[i].bin << std::endl;
  }
  std::cout << "#------------------------------------------------" << std::endl;;
  std::cout << std::endl;

  std::cout << "#---------- photoion_cs ---------------------------------" << std::endl;;
  for (int i=0;i<atom->photo_cs_.size();++i)
  {
    std::cout << "# photocs: " << atom->photo_cs_[i].id << std::endl;
    std::cout << "# ----------------------" << std::endl;
    for (int j=0;j<atom->photo_cs_[i].n_pts;++j)
    {
      int ishift = j + atom->photo_cs_[i].i_start;
      std::cout << nu_grid_.center(ishift) << "\t" << atom->photo_cs_[i].s[j] << std::endl;
    }
    std::cout << "#------------------------------------------------" << std::endl;;
  }
  std::cout << std::endl;
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

  // return if this atom has already been read
  if (atomlist_[z].data_exists_)
    return 0;

//  atomlist_[z].max_ion_stage_ = max_ion;
  atomlist_[z].nu_grid_ = &nu_grid_;

  if (datafile_version_ == 1)
    return read_newstyle_atomic_data(z);
  else
    return read_oldstyle_atomic_data(z);

}

//------------------------------------------------------------------------
// Read all atomic data for species with atomic number Z
// Maximum ionization state no specified so include allow
// ionization states in the atomic data file
//------------------------------------------------------------------------
int AtomicData::read_atomic_data(int z)
{
  return read_atomic_data(z,z);
}


//------------------------------------------------------------------------
// Read  atomic data for species with atomic number Z
// using a file in the new style of formatting
//------------------------------------------------------------------------
int AtomicData::read_newstyle_atomic_data(int z)
{
  // open hdf5 file
  herr_t status;
  status = H5Eset_auto1(NULL, NULL);
  hid_t file_id = H5Fopen (atom_datafile_.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  // get element group name
  char atomname[100];
  sprintf(atomname,"%d/",z);
  char dset[5000];

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
      atom->ions_.resize(1);
      atom->ions_[0].chi    = 999999;
      atom->ions_[0].ground = 0;
      atom->ions_[0].stage  = 0;
      atom->levels_.resize(1);
      atom->levels_[0].E  = 0;
      atom->levels_[0].g  = 1;
      atom->levels_[0].ic = -1;
      atom->levels_[0].ion = 0;
      return 1;
  }

  // otherwise data exists, we'll read it in
  atomlist_[z].data_exists_ = true;

  // determine number of ions
  int n_tot_ions = 0;
  for (int ion=0;ion < z;++ion)
  {
    char ionname[100];
    sprintf(ionname,"%d/%d/",z,ion);
    status = H5Gget_objinfo (file_id, ionname, 0, NULL);
    if (status == 0)
      n_tot_ions++;
  }
  atom->n_ions_ = n_tot_ions;
  if (atom->n_ions_  > atom->max_ion_stage_)
    atom->n_ions_  = atom->max_ion_stage_;
  atom->ions_.resize(atom->n_ions_);


  // ----------------------------------------
  // loop over ions and read data
  // ----------------------------------------
  int lev_count  = 0;
  for (int ion=0;ion < atom->n_ions_;++ion)
  {
    int ibase = lev_count;
    atom->ions_[ion].ground = ibase;
    atom->ions_[ion].stage  = ion;

    char ionname[100];
    sprintf(ionname,"%d/%d/",z,ion);
    status = H5Gget_objinfo (file_id, ionname, 0, NULL);
    if (status != 0)
    {
      std::cout << "# ERROR: Missing data for ion Z = " << z << ", I = " << ion << std::endl;
      return -1;
    }

    // get ionization potential
    sprintf(dset,"%s/ion_chi",ionname);
    double chi;
    status = H5LTread_dataset_double(file_id,dset,&chi);
    if (status != 0) return -1;
    atom->ions_[ion].chi = chi;

    // get number of levels
    int tot_n_levels;
    sprintf(dset,"%s/n_levels",ionname);
    status = H5LTread_dataset_int(file_id, dset,&tot_n_levels);
    if (status != 0) return -1;

    // read level statistical weights
    int    *g_iarr = new int[tot_n_levels];
    sprintf(dset,"%s%s",ionname,"level_g");
    status = H5LTread_dataset_int(file_id,dset,g_iarr);
    if (status != 0) return -1;

    // read photion_cs indices
    int    *cs_iarr = new int[tot_n_levels];
    sprintf(dset,"%s%s",ionname,"level_cs");
    status = H5LTread_dataset_int(file_id,dset,cs_iarr);
    if (status != 0) return -1;

    // read level excitation energy
    double *E_darr = new double[tot_n_levels];
    sprintf(dset,"%s%s",ionname,"level_E");
    status = H5LTread_dataset_double(file_id,dset,E_darr);
    if (status != 0) return -1;

    // Cap number of levels
    if (tot_n_levels > atom->max_n_levels_)
        tot_n_levels = atom->max_n_levels_;

    // add in levels
    for (int i=0;i<tot_n_levels;++i)
    {
      AtomicLevel lev;
      lev.ion = ion;
      lev.g = g_iarr[i];
      if (lev.g == 0) lev.g = 1;
      lev.E = E_darr[i];
      lev.E_ion = chi - lev.E;
      lev.ic = ibase + tot_n_levels;
      lev.cs = cs_iarr[i];

      // store hydrogenic photo cs stuff
      lev.nu_t = lev.E_ion*pc::ev_to_ergs/pc::h;
      lev.sigma0 = 6.3e-18;
            // effective excitation quantum number
      //double E_ground = atom->ions_[atom->levels_[i].ion].chi;
      //double n_eff = pow(1 - (E_ground - E_ion)/E_ground,-0.5);
      //double s_fac = n_eff/(atom->levels_[i].ion+1)/(atom->levels_[i].ion + 1);
      atom->levels_.push_back(lev);
    }

    // clean up
    delete[] g_iarr;
    delete[] cs_iarr;
    delete[] E_darr;


    // ---------------------------------------
    // Read and Setup line data
    // ---------------------------------------
    int n_tot_lines;
    sprintf(dset,"%s/n_lines",ionname);
    status = H5LTread_dataset_int(file_id, dset,&n_tot_lines);
    if (status != 0) n_tot_lines = 0;

    if (n_tot_lines > 0)
    {
      double *A_darr   = new double[n_tot_lines];
      int *lu_iarr     = new int[n_tot_lines];
      int *ll_iarr     = new int[n_tot_lines];
      int *col_iarr    = new int[n_tot_lines];

      // read line upper level indices
      sprintf(dset,"%s%s",ionname,"line_u");
      status = H5LTread_dataset_int(file_id,dset,lu_iarr);
      if (status != 0) return -1;

      // read line lower level indices
      sprintf(dset,"%s%s",ionname,"line_l");
      status = H5LTread_dataset_int(file_id,dset,ll_iarr);
      if (status != 0) return -1;

      // read line Einstein A
      sprintf(dset,"%s%s",ionname,"line_A");
      status = H5LTread_dataset_double(file_id,dset,A_darr);
      if (status != 0) return -1;

      // read line index of collisional rate
      sprintf(dset,"%s%s",ionname,"line_col_id");
      status = H5LTread_dataset_int(file_id,dset,col_iarr);
      if (status != 0) return -1;

      // add in the lines
      int n_lines_add = 0;
      for (int i=0;i<n_tot_lines;++i)
      {
        int ll = ll_iarr[i];
        int lu = lu_iarr[i];
        double A = A_darr[i];

        // skip lines that are to omitted levels
        if (lu >= tot_n_levels) continue;
        if (ll >= tot_n_levels) continue;
        n_lines_add++;

        // reindex levels to account for all ion stages
        lu   = lu + atom->ions_[ion].ground;
        ll   = ll + atom->ions_[ion].ground;

        // create a new atomic line to store this data
        AtomicLine lin;
        lin.lu   = lu;
        lin.ll   = ll;
        lin.A_ul = A;

        // get wavelength/frequency
        double delta_E = atom->levels_[lu].E - atom->levels_[ll].E;
        double nu      = delta_E*pc::ev_to_ergs/pc::h;
        double lam   = pc::c/nu*pc::cm_to_angs;
        lin.nu  = nu;

        // set Einstein Coeficients
        int gl = atom->levels_[ll].g;
        int gu = atom->levels_[lu].g;
        lin.B_ul = A*pc::c*pc::c/2.0/pc::h/nu/nu/nu;
        lin.B_lu = lin.B_ul*gu/gl;

        // set oscillator strength (see e.g., Rutten page 24)
        double lam_cm = lam*pc::angs_to_cm;
        lin.f_lu = lam_cm*lam_cm*A*gu/gl/(8*pc::pi*pc::sigma_tot);

        // find index of bin in frequency grid
        // This function will return -1 if out of grid bounds
        lin.bin = nu_grid_.locate_check_bounds(nu);

        // read collisional rate if it exists
        lin.col_rate = NULL;
        if (col_iarr[i] != -1)
        {
          lin.col_rate = new AtomicCollisionalBB_Rate;

          int n_pts;
          char cname[1000];
          sprintf(cname,"%s/collisional_data/col_%d",ionname,col_iarr[i]);
          sprintf(dset,"%s/n_pts",cname);
          status = H5LTread_dataset_int(file_id,dset,&n_pts);

          lin.col_rate->E = delta_E;
          double* tmp = new double[n_pts];
          sprintf(dset,"%s/T",cname);
          status = H5LTread_dataset_double(file_id,dset,tmp);
          lin.col_rate->T.assign(tmp,tmp+n_pts);
          sprintf(dset,"%s/Omega",cname);
          status = H5LTread_dataset_double(file_id,dset,tmp);
          lin.col_rate->O.assign(tmp,tmp+n_pts);
          delete[] tmp;
        }

        // if in bounds, add into vector
        if (lin.bin >= 0)
          atom->lines_.push_back(lin);

      }
      delete[] ll_iarr;
      delete[] lu_iarr;
      delete[] A_darr;
    }
    lev_count += tot_n_levels;

    // ---------------------------------------
    // Read and Setup photoion cross-section data
    // ---------------------------------------
    int n_photo_cs;
    sprintf(dset,"%s/photoion_data/n_photo_cs",ionname);
    status = H5LTread_dataset_int(file_id, dset,&n_photo_cs);
    if (status != 0) n_photo_cs = 0;

    for (int i=0;i<n_photo_cs;++i)
    {
      // check that this cross-section is used by some level
      int used = 0;
      for (int ilev=0;ilev<atom->levels_.size();++ilev)
        if (atom->levels_[ilev].cs  == i)
          used = 1;
      if (used == 0) continue;

      // set up to read in cross-section data
      sprintf(dset,"%s/photoion_data/cs_%d/n_pts",ionname,i);
      int n_pts;
      status = H5LTread_dataset_int(file_id, dset,&n_pts);
      xy_array cs_data(n_pts);

      // read the energy array
      double *darray = new double[n_pts];
      sprintf(dset,"%s/photoion_data/cs_%d/E",ionname,i);
      status = H5LTread_dataset_double(file_id, dset,darray);
      for (int j=0;j<n_pts;++j)
        cs_data.x[j] = darray[j]*pc::ev_to_ergs/pc::h;

      // read the cross-section array
      sprintf(dset,"%s/photoion_data/cs_%d/sigma",ionname,i);
      status = H5LTread_dataset_double(file_id, dset, darray);
      for (int j=0;j<n_pts;++j)
        cs_data.y[j] = darray[j];

      // edge data values
      double nu_edge = cs_data.x[n_pts-1];
      double s_edge = cs_data.y[n_pts-1];

      // create a new cross-section object to store this in
      AtomicPhotoCS this_cs;
      this_cs.id = i;

      // find the lowest frequency defined at
      double nu_t = cs_data.x[0];
      this_cs.i_start = nu_grid_.locate(nu_t);
      this_cs.n_pts = nu_grid_.size() - this_cs.i_start;

      this_cs.s.resize(this_cs.n_pts);
      for (int j=0;j<this_cs.n_pts;++j)
      {
        int ishift = j + this_cs.i_start;
        double nu = nu_grid_.center(ishift);
        if (nu < nu_edge) {
          this_cs.s[j] = cs_data.value_at_with_zero_edges(nu); }
        else {
          this_cs.s[j] = s_edge*pow(nu_edge/nu,3.0);  }
      }
      // store this cross-section and clean up
      atom->photo_cs_.push_back(this_cs);
      delete [] darray;
    }


  // finish loop for reading of ion's data
   }

   // add in extra fully ionized state and level
   AtomicIon aion;
   aion.chi = 99999;
   aion.ground = lev_count;
   aion.stage  = atom->ions_.size();
   atom->ions_.push_back(aion);

   AtomicLevel lev;
   lev.ic = -1;
   lev.E  = 0.0;
   lev.g  = 1;
   lev.E_ion = 99999;
   lev.ion  = atom->n_ions_;
   atom->levels_.push_back(lev);

   atom->n_levels_ = atom->levels_.size();
   atom->n_ions_   = atom->ions_.size();
   atom->n_lines_  = atom->lines_.size();

   // -------------------------------------
   // read non_thermal ion cross-section data
   // -------------------------------------

   atom->nt_ion_cs_.resize(atom->n_ions_);
   for (int i=0;i<atom->n_ions_;++i)
     atom->nt_ion_cs_[i].n_shells = 0;

   int n_data=0;
   status = H5LTread_dataset_int(file_id,"non_thermal_cs/n_data",&n_data);
   if (status == 0)
   {
     int* Z_data = new int[n_data];
     int* I_data = new int[n_data];
     double* A_data = new double[n_data];
     double* B_data = new double[n_data];;
     double* C_data = new double[n_data];
     double* D_data = new double[n_data];
     double* chi_data = new double[n_data];
     status = H5LTread_dataset_int(file_id,"non_thermal_cs/Z",Z_data);
     status = H5LTread_dataset_int(file_id,"non_thermal_cs/ion",I_data);
     status = H5LTread_dataset_double(file_id,"non_thermal_cs/A",A_data);
     status = H5LTread_dataset_double(file_id,"non_thermal_cs/B",B_data);
     status = H5LTread_dataset_double(file_id,"non_thermal_cs/C",C_data);
     status = H5LTread_dataset_double(file_id,"non_thermal_cs/D",D_data);
     status = H5LTread_dataset_double(file_id,"non_thermal_cs/chi",chi_data);

     for (int i=0;i<n_data;i++)
     {
       int ion = I_data[i];
       if (Z_data[i] != z) continue;
       if (ion >= atom->n_ions_) continue;
       atom->nt_ion_cs_[ion].A.push_back(A_data[i]);
       atom->nt_ion_cs_[ion].B.push_back(B_data[i]);
       atom->nt_ion_cs_[ion].C.push_back(C_data[i]);
       atom->nt_ion_cs_[ion].D.push_back(D_data[i]);
       atom->nt_ion_cs_[ion].chi.push_back(chi_data[i]);
       atom->nt_ion_cs_[ion].n_shells++;
     }
     delete[] A_data;
     delete[] B_data;
     delete[] C_data;
     delete[] D_data;
     delete[] chi_data;
     delete[] Z_data;
     delete[] I_data;
   }

   // all done
  H5Fclose(file_id);

/*  for (int i=0;i<atom->n_lines_;++i)
  {
    if (atom->lines_[i].col_rate == NULL)
      std::cout << i << "\t" << "NULL\n";
    else
    {
      std::cout << i  << " " << atom->lines_[i].col_rate->E << "\n";
      for (int j=0;j<atom->lines_[i].col_rate->T.size();++j)
        std::cout << atom->lines_[i].col_rate->T[j] << " " << atom->lines_[i].col_rate->O[j] << "\n";
    }
  }*/
  return 0;


}



//------------------------------------------------------------------------
// Read all atomic data for species with atomic number Z
//------------------------------------------------------------------------
int AtomicData::read_oldstyle_atomic_data(int z)
{

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
      atom->ions_.resize(1);
      atom->ions_[0].chi    = 999999;
      atom->ions_[0].ground = 0;
      atom->ions_[0].stage  = 0;
      atom->levels_.resize(1);
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
  if (status != 0) return 3;
  atom->n_ions_ = n_tot_ions;
  int add_last_stage = 0;
  if (atom->n_ions_  > atom->max_ion_stage_)
  {
    atom->n_ions_  = atom->max_ion_stage_;
    atom->ions_.resize(atom->n_ions_ + 1);
    add_last_stage = 1;
  }
  else
    atom->ions_.resize(atom->n_ions_);

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
  if (status != 0) return 2;

  // make space to read data
  double *lev_darr = new double[tot_n_levels];
  int *lev_iarr = new int[tot_n_levels];

  // read level ionization stages
  sprintf(dset,"%s%s",atomname,"level_i");
  status = H5LTread_dataset_int(file_id,dset,lev_iarr);
  if (status != 0) return 2;

  // count number of levels to use
  int i;
  for (i=0;i<tot_n_levels;++i)
    if (lev_iarr[i] > atom->n_ions_-1) break;
  atom->n_levels_ = i;

  if (add_last_stage)
    atom->levels_.resize(atom->n_levels_ + 1);
  else
    atom->levels_.resize(atom->n_levels_);

  for (int i=0;i<atom->n_levels_;++i)
    atom->levels_[i].ion = lev_iarr[i];

  // read level statistical weights
  sprintf(dset,"%s%s",atomname,"level_g");
  status = H5LTread_dataset_int(file_id,dset,lev_iarr);
  if (status != 0) return 2;
  for (int i=0;i<atom->n_levels_;++i) atom->levels_[i].g = lev_iarr[i];

  // read level excitation energy
  sprintf(dset,"%s%s",atomname,"level_E");
  status = H5LTread_dataset_double(file_id,dset,lev_darr);
  if (status != 0) return 2;
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

      // set photoionization cross-section to just hydrogenic
      atom->levels_[i].cs = -1;
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
    atom->lines_.resize(atom->n_lines_);

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
  // photoionization cross-sections
  // for oldstyle files are just hydrogenic
  // ----------------------------------------
  for (int i=0;i<atom->n_levels_;++i)
  {
     // set photoionization cross-section
     double E_ion = atom->levels_[i].E_ion;

     // effective excitation quantum number
     //double E_ground = atom->ions_[atom->levels_[i].ion].chi;
     //double n_eff = pow(1 - (E_ground - E_ion)/E_ground,-0.5);
     //double s_fac = n_eff/(atom->levels_[i].ion+1)/(atom->levels_[i].ion + 1);
    double nu_t = E_ion*pc::ev_to_ergs/pc::h;
    atom->levels_[i].nu_t   = nu_t;
    atom->levels_[i].sigma0 =  6.3e-18;

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
