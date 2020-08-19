#include <cstdlib>
#include <algorithm>
#include "physical_constants.h"
#include "GasState.h"
#include "brent.h"
#include "brent.cpp"
#include <iostream>
#include <fstream>
#include <iomanip>

#ifdef MPI_PARALLEL
#include "mpi.h"
#endif

namespace pc = physical_constants;

//----------------------------------------------------------------
// simple constructor
//----------------------------------------------------------------
GasState::GasState()
{
  nlte_turned_on_ = 0;
  e_gamma = 0;
  no_ground_recomb = 0;
  line_velocity_width_ = 0;
}

//----------------------------------------------------------------
// initialize the gas by specifying the atoms that will
// compose it, the freq. array, and atomic data pointer.
// Assumes the atomic data has already been read in and
// is passed as a pointer to the class
// inputs:
// AtomicData* adata:  pointer to atomic data holder class
// std::vector<int> e:  vector of atomic numbers
// std::vector<int> A:  vector of atomic weights (in atomic units)
// locate_array ng:  locate_array giving the freq. array
//---------------------------------------------------------------
void GasState::initialize
(AtomicData* adata, std::vector<int> e, std::vector<int> A, locate_array ng)
{
  // verbocity
#ifdef MPI_PARALLEL
  int my_rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
  verbose_ = (my_rank == 0);
#else
  verbose_ = 1;
#endif

  // copy the frequency grid
  nu_grid_.copy(ng);

  // set passed variables
  for (size_t i=0;i<e.size();++i) elem_Z.push_back(e[i]);
  for (size_t i=0;i<e.size();++i) elem_A.push_back(A[i]);
  mass_frac.resize(e.size());
  atomic_data_ = adata;


  atoms.resize(elem_Z.size());
  for (size_t i=0;i<atoms.size();++i)
  {
    int error = atomic_data_->read_atomic_data(elem_Z[i]);
    if ((error)&&(verbose_))
      std::cerr << "# ERROR: incomplete data (" << error << ") for atom Z=" << elem_Z[i] <<
        " in file " << atomdata_file_ << std::endl;
    atoms[i].initialize(elem_Z[i],atomic_data_);
  }
}

//----------------------------------------------------------------
// initialize the gas by specifying the atoms that will
// compose it and freq. array, and name of atomic data file
// will initialize the atomic data locally
// inputs:
// std::string atomfile: name of atom data file (in hdf5)
// std::vector<int> e:  vector of atomic numbers
// std::vector<int> A:  vector of atomic weights (in atomic units)
// locate_array ng:  locate_array giving the freq. array
//---------------------------------------------------------------
void GasState::initialize
(std::string af, std::vector<int> e, std::vector<int> A, locate_array ng)
{
  atomdata_file_ = af;
  // check if atomfile is there
  std::ifstream afile(atomdata_file_);
  if (!afile)
  {
    if (verbose_)
      std::cerr << "Can't open atom datafile " << atomdata_file_ << "; exiting" << std::endl;
    exit(1);
  }
  afile.close();

  // read in the atom data
  atomic_data_ = new AtomicData;
  atomic_data_->initialize(atomdata_file_,ng);
  initialize(atomic_data_,e,A,ng);
}



//-----------------------------------------------------------------
// Set the atoms to be solved in nlte
//-----------------------------------------------------------------
void GasState::set_atoms_in_nlte
(std::vector<int> useatoms)
{
  if (useatoms.size() == 0)
    nlte_turned_on_ = 0;
  else
    nlte_turned_on_ = 1;

  for (int i=0;i<useatoms.size();++i)
  {
    for (int j=0;j<atoms.size();++j)
    {
      if (elem_Z[j] == useatoms[i])
      {
        atoms[j].set_use_nlte();
	      atoms[j].use_collisions_nlte_ = use_collisions_nlte_;
      }
    }
  }
}

//-----------------------------------------------------------------
// Set mass fractions of each element in the gas
// this function will enforce that the mass fractions are
// normalize (i.e., add up to 1)
//
// input:
// std::vector<double> x: vector of mass fractions of each element
//-----------------------------------------------------------------
void GasState::set_mass_fractions(std::vector<double>& x)
{
  double norm = 0.0;
  for (size_t i=0;i<mass_frac.size();++i)
  {
    mass_frac[i] = x[i];
    norm += x[i];
  }

  // make sure it is normalized
  for (size_t i=0;i<mass_frac.size();++i) mass_frac[i] /= norm;

  // find mean weight of atoms/ions in gas (not counting free electrons)
  //  see e.g. Hansen Kawaler Trimble textbook 2nd edition page 18
  double inverse_mu_sum = 0.;
  for (size_t i=0;i<mass_frac.size();++i)
    inverse_mu_sum += mass_frac[i]/elem_A[i];

  this->mu_I = 1./inverse_mu_sum;
}


//-----------------------------------------------------------
// read fuzz lines from a file
// input:
// std::string fuzzfile: name of hdf5 file with fuzz data
// The lines will only be read for atoms that have
// already had data read in
//-----------------------------------------------------------
int GasState::read_fuzzfile(std::string fuzzfile)
{
  // check if fuzzfile exists
  FILE *fin = fopen(fuzzfile.c_str(),"r");
  if ((fin == NULL)&&(verbose_)&&(fuzzfile != ""))
    std::cerr << "# Warning: Can't open atomic data fuzzfile: "
    << fuzzfile << std::endl;
  if (fin == NULL) return 0;

  int n_lines = atomic_data_->read_fuzzfile_data(fuzzfile);
  return n_lines;
}


//-----------------------------------------------------------
// return the ionization state, i.e., electron density
// over total ion number density
//-----------------------------------------------------------
double GasState::get_ionization_state()
{
  double ni = dens_/(mu_I*pc::m_p);
  return n_elec_/ni;
}

//-----------------------------------------------------------
// Solve for the gas state (excitation/ionization)
// the level populations will be stored internally for
// Assumes no radiation field given, so in LTE
// returns: any error
//-----------------------------------------------------------
int GasState::solve_state()
{
  if (nlte_turned_on_)
  {
    std::cout << "ERROR; no radiation field Jnu passed to GasState NLTE solve\n";
    return 1;
  }
  std::vector<SedonaReal> J_nu;
  return solve_state(J_nu);
}

//-----------------------------------------------------------
// Solve for the gas state (excitation/ionization)
// the level populations will be stored internally for
// further calculations
// Returns: any error
//-----------------------------------------------------------
int GasState::solve_state(std::vector<SedonaReal>& J_nu)
{

  // set key properties of all atoms
  for (size_t i=0;i<atoms.size();++i)
  {
    atoms[i].n_dens_  = dens_*mass_frac[i]/(elem_A[i]*pc::m_p);
    atoms[i].e_gamma_ = e_gamma*mass_frac[i];
    atoms[i].no_ground_recomb_ = no_ground_recomb;
    atoms[i].gas_temp_ = temp_;
    // line widths
    double vd = sqrt(2*pc::k*temp_/pc::m_p/elem_A[i]);
    if (line_velocity_width_ > 0) vd = line_velocity_width_;
    atoms[i].line_beta_dop_ = vd/pc::c;

    // Calculate radiative rates, if using NLTE
    // and this atom to be treated in NLTE
    if (nlte_turned_on_ && atoms[i].use_nlte_)
        atoms[i].calculate_radiative_rates(J_nu);
  }

  // solve for the state of the gas.
  // This will solve the ion/excite state of each
  // atom and then iterate electron density until
  // charge conservation is fulfilled.

  gasMemFn f = &GasState::charge_conservation;
  solve_error_  = 0;

  charge_brent_args brent_args;

  brent_solver<GasState,charge_brent_args> solver;
  int n; // will store number of brent solver iterations

  // set bounds and tolerance for solving electron density
  double max_ne = 100*dens_/(mu_I*pc::m_p);
  double min_ne = 1e-10*dens_/(mu_I*pc::m_p);
  double tol    = 1e-3;

  int max_iters = 100;
  // call brent method to solve for n_e
  n_elec_ = solver.solve(*this, f, &brent_args, min_ne, max_ne,tol,max_iters,&n);

  if (n == -1) solve_error_ = 1;
  if (n == -2) solve_error_ = 2;


  return solve_error_;

}

//-----------------------------------------------------------
// This is the function that represents the non-linear equation
// for N_e.  It is used in Brents method below to solve
// for the root, thus determining N_e.  This equation is
// basically just the one for charge conservation.
//-----------------------------------------------------------
double GasState::charge_conservation(double ne, charge_brent_args* args) 
{
  // brent_args not actually used here
  
  // start with charge conservation function f set to zero
  double f  = 0;
  // loop over all atoms
  for (size_t i=0;i<atoms.size();++i)
  {
    // Solve the state of the atome with this value of
    // electron density ne.  Solution will either be in
    // LTE or NLTE depending on if atoms[i]->use_nlte_ is set
    if (nlte_turned_on_ && atoms[i].use_nlte_)
      atoms[i].solve_nlte(ne);
    else
      atoms[i].solve_lte(ne);

    // total electron donation from this atomic species
    f += dens_*mass_frac[i]/(elem_A[i]*pc::m_p)*atoms[i].get_net_ion_fraction();
  }

  // total ionization density minus electron density should equal zero
  f = f - ne;
  // return to non-linear solver to iterate this to zero
  return f;
}




//-----------------------------------------------------------
// return the fraction of atoms of index i that are in
// ionization state j.
//-----------------------------------------------------------
double GasState::get_ionization_fraction(int i, int j)
{
  if ((i < 0)||(i >= (int)atoms.size())) return -1;
  if ((j < 0)||(j >= atoms[i].n_ions_))   return -1;
  return atoms[i].ionization_fraction(j);
}



//-----------------------------------------------------------
// return the fraction of atoms of index i that are in
// ionization state j.
//-----------------------------------------------------------
double GasState::get_level_fraction(int i, int j)
{
  if ((i < 0)||(i >= (int)atoms.size())) return -1;
  if ((j < 0)||(j >= atoms[i].n_levels_)) return -1;
  return atoms[i].level_fraction(j);
}

//-----------------------------------------------------------
// return the fraction of atoms of index i that are in
// ionization state j.
//-----------------------------------------------------------
double GasState::get_level_departure(int i, int j)
{
  if ((i < 0)||(i >= (int)atoms.size())) return -1;
  if ((j < 0)||(j >= atoms[i].n_levels_)) return -1;
  return atoms[i].level_depature(j);
}

//-----------------------------------------------------------
// print out of the gas properties
//-----------------------------------------------------------
void GasState::print_properties()
{

  std::cout << "#-------------------------------------------------" << std::endl;
  std::cout << "# atomic data from file: " << std::endl;
  std::cout << "#   -> " << atomic_data_->get_input_filename() << std::endl;
  std::cout << "# File is sedona hdf5 format version = " << atomic_data_->get_version() << std::endl;
  std::cout << "# " << std::endl << "# Atomic species and data to be used are:" << std::endl;
  std::cout << "#--------------------------------------------------" << std::endl;
  std::cout << "#  Z    n_ions  n_levels  n_lines  n_fuzz_lines" << std::endl;
  std::cout << "#-------------------------------------------------" << std::endl;
  for (size_t i=0;i<atoms.size();++i)
  {
    printf("# %2d.%d ",elem_Z[i],elem_A[i]);
    printf(" %4d %8d  %8d  %8d",atoms[i].n_ions_,atoms[i].n_levels_,atoms[i].n_lines_,
        atoms[i].get_n_fuzz_lines());
    printf("\n");
  }
  std::cout << "#-------------------------------------------------" << std::endl;
  std::cout << "#" << std::endl;


  std::cout << "# opacity settings\n";
  std::cout << "#-------------------------------------------------\n";

  if ( (bulk_grey_opacity_ != 0) || (use_zone_specific_grey_opacity_ != 0) )
  {
    std::cout << "# bulk grey opacity              = " << bulk_grey_opacity_ << "\n";
    std::cout << "# use_zone_specific_grey_opacity = " << use_zone_specific_grey_opacity_ << "\n";
  }
  else
  {
    std::cout << "# use_e_scattering    = " << use_electron_scattering_opacity << "\n";
    std::cout << "# use_free_free       = " << use_free_free_opacity << "\n";
    std::cout << "# use_bound_free      = " << use_bound_free_opacity << "\n";
    std::cout << "# use_bound_bound     = " << use_bound_bound_opacity << "\n";
    std::cout << "# use_line_exp        = " << use_line_expansion_opacity << "\n";
    std::cout << "# use_fuzz_exp        = " << use_fuzz_expansion_opacity << "\n";

    std::cout << "# use_nlte            = " << nlte_turned_on_ << "\n";
    std::cout << "# use_collisions_nlte = " << use_collisions_nlte_ << "\n";
    if (nlte_turned_on_)
    {
      int n_in_nlte = 0;
      std::cout << "# These atoms to be treated in non-LTE: ";
      for (int i=0;i<atoms.size();i++)
        if (atoms[i].use_nlte_)
        {
          n_in_nlte += 1;
          std::cout << elem_Z[i] << " ";
        }
      if (n_in_nlte == 0)
        std::cout << "None";
      std::cout << "\n";

#ifdef USE_EIGEN
      std::cout << "# Will solve NLTE matrix equation with eigen\n";
#else
      std::cout << "# Will solve NLTE matrix equation with GSL\n";
#endif

    }
  }
  std::cout << "#---------------------------------------" << std::endl;
}

//-----------------------------------------------------------
// print out of the gas class
//-----------------------------------------------------------
void GasState::print()
{
  std::cout << "# dens = " << dens_ << "\n";
  std::cout << "# temp = " << temp_ << "\n";
  std::cout << "# mu_I = " << mu_I << "\n";
  for (size_t i=0;i<elem_Z.size();++i)
    printf("%4d %12.4e %12.4e\n",elem_Z[i],mass_frac[i],atoms[i].get_net_ion_fraction());
  for (size_t i=0;i<elem_Z.size();++i)
    atoms[i].print();
}


void GasState::write_levels(int iz)
{
  char zonefile[1000];
  sprintf(zonefile,"zone_%d.h5",iz);
  hid_t file_id = H5Fcreate(zonefile,H5F_ACC_TRUNC, H5P_DEFAULT,  H5P_DEFAULT);

  const int RANK = 1;

  for(size_t j=0;j<atoms.size();j++)
  {
    char afile[100];
    int this_Z = elem_Z[j];
    sprintf(afile,"Z_%d",this_Z);
    hid_t atom_id = H5Gcreate1(file_id,afile,0);


    double* tmp_ion = new double[elem_Z[j]+1];
    hsize_t dims_ion[RANK]={(hsize_t)elem_Z[j]+1};
    for(int k=0;k<elem_Z[j]+1;k++)
      tmp_ion[k] = get_ionization_fraction(j,k);
    H5LTmake_dataset(atom_id,"ion_fraction",RANK,dims_ion,H5T_NATIVE_DOUBLE,tmp_ion);

    int this_nl = atoms[j].n_levels_;
    double* tmp_level = new double[this_nl];
    hsize_t dims_level[RANK]={(hsize_t)this_nl};
    for(int k=0;k<this_nl;k++)
      tmp_level[k] = get_level_fraction(j,k);
    H5LTmake_dataset(atom_id,"level_fraction",RANK,dims_level,H5T_NATIVE_DOUBLE,tmp_level);

    for(int k=0;k<this_nl;k++)
      tmp_level[k] = get_level_departure(j,k);
    H5LTmake_dataset(atom_id,"level_departure",RANK,dims_level,H5T_NATIVE_DOUBLE,tmp_level);

    int this_nd = 1;
    double* tmp_ndens = new double[this_nd];
    hsize_t dims_ndens[RANK] = {(hsize_t)this_nd};
    tmp_ndens[0] = dens_ * mass_frac[j]/(elem_A[j]*pc::m_p);
    H5LTmake_dataset(atom_id,"n_dens",RANK,dims_ndens,H5T_NATIVE_DOUBLE,tmp_ndens);

    H5Gclose(atom_id);
    delete[] tmp_level;
    delete[] tmp_ion;
  }
  H5Fclose(file_id);

}

std::string format_with_commas(long int value)
{
    std::string numWithCommas = std::to_string(value);
    int insertPosition = numWithCommas.length() - 3;
    while (insertPosition > 0)
    {
        numWithCommas.insert(insertPosition, ",");
        insertPosition-=3;
    }
    return numWithCommas;
}

double GasState::print_memory_footprint()
{
    long int n_tot_lines  = 0;
    long int n_tot_levels = 0;
    for (int i=0;i<atoms.size();i++)
    {
        n_tot_lines  += atoms[i].n_lines_;
        n_tot_levels += atoms[i].n_levels_;
    }
    long int size_line = sizeof(AtomicLine);
    long int size_level = sizeof(AtomicLevel);

    long int tot_line  = size_line*n_tot_lines;
    long int tot_level = size_level*n_tot_levels;
    long int total = tot_line + tot_level;

    std::cout << "#-----------------------------------------------------|\n";
    std::cout << std::setw(10) << "#  data   |";
    std::cout << std::setw(15) << " number |";
    std::cout << std::setw( 12) << " each (B) |";
    std::cout << std::setw( 18) << " total (B) |\n";
    std::cout << "#-----------------------------------------------------|\n";
    std::cout << std::setw(10) << "# lines   |";
    std::cout << std::setw(15) << format_with_commas(n_tot_lines) + " |";
    std::cout << std::setw(12) << format_with_commas(size_line) + " |";
    std::cout << std::setw(18) << format_with_commas(tot_line) + " |\n";
    std::cout << std::setw(10) << "# levels  |";
    std::cout << std::setw(15) << format_with_commas(n_tot_levels) + " |";
    std::cout << std::setw(12) << format_with_commas(size_level) + " |";
    std::cout << std::setw(18) << format_with_commas(tot_level) + " |\n";
    std::cout << std::setw(10) << "# total   |";
    std::cout << std::setw(15) <<  " |";
    std::cout << std::setw(12) <<  " |";
    std::cout << std::setw(18) << format_with_commas(total) + " |\n";
    std::cout << "#-----------------------------------------------------|\n";

    return total;
}
