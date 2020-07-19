#ifndef _ATOMIC_DATA_H
#define _ATOMiC_DATA_H 1

#define MAX_N_ATOMS 120

#include <string>
#include <vector>
#include "sedona.h"
#include "xy_array.h"
#include "locate_array.h"

struct fuzz_line_structure
{
  int n_lines;
  std::vector<double> nu;
  std::vector<double> El;
  std::vector<double> gf;
  std::vector<int>   ion;
  std::vector<int>   bin;
};


struct AtomicIon
{
  int stage;          // ionization stage (0 = neutral, 1 = +, etc..)
  int ground;         // index of ground state level
  double chi;         // ionization energy, in eV
};

struct AtomicLine
{
  int lu,ll;           // index of upper/lower level
  double nu;           // rest frequency (Hz)
  double f_lu;         // oscillator strength
  double A_ul;         // Einstein A coeficient
  double B_ul;         // Einstein B coeficient
  double B_lu;         // Einstein B coeficient
  int    bin;          // index of the nu grid bin
};

struct AtomicLevel
{
  int   globalID;       // global id
  int  ion;             // ionization state (0 = neutral)
  int   ic;             // index of level we ionize to (-1 if none)
  int    g;             // statistical weight
  int   cs;             // index of photoionization cross-section
  double E;             // excitation energy above ground (in eV)
  double E_ion;         // energy to ionize (in eV)

  // photoionization cross-section as a function of wavelength
  xy_array s_photo;
  // recombination coefficient as a function of temperature
  xy_array a_rec;

};

struct AtomicPhotoCS
{
    int id;
    int n_pts;
    std::vector<double> E;
    std::vector<double> s;
};


//---------------------------------------------
// Holds the data for an individual atoms
//---------------------------------------------
class IndividualAtomData
{

public:

  // name of the data file that we read from
  std::string atom_datafile_;

  // bool to see if data was actually read
  bool data_exists_;

  std::vector<AtomicLevel>   levels_;      // array of level data
  std::vector<AtomicLine>    lines_;       // array of line data
  std::vector<AtomicIon>     ions_;        // array of ion data
  std::vector<AtomicPhotoCS> photo_cs_;    // array of cross-section

  int n_ions_;              // Number of ionic stages considered
  int n_levels_;            // number of energy levels
  int n_lines_;             // number of line transitions
  int max_ion_stage_;       // maximum ionization stage to use
  int max_n_levels_;        // maximum number of levels per ion stage

  fuzz_line_structure fuzz_lines_; // vector of fuzz lines

  double get_ion_chi(int i) {
    return ions_[i].chi;
  }
  int get_lev_g(int i) {
    return levels_[i].g;
  }
  double get_lev_E(int i) {
    return levels_[i].E;
  }
  int get_lev_ion(int i) {
    return levels_[i].ion;
  }
  int get_lev_ic(int i) {
    return levels_[i].ic;
  }
  double get_lev_Eion(int i) {
    return levels_[i].E_ion;
  }
  double get_lev_photo_cs(int i,double E_ev) {
    return levels_[i].s_photo.value_at_with_zero_edges(E_ev);
  }
  double get_line_nu(int i) {
    return lines_[i].nu;
  }
  double get_line_A(int i) {
    return lines_[i].A_ul;
  }
  double get_line_Bul(int i) {
    return lines_[i].B_ul;
  }
  double get_line_Blu(int i) {
    return lines_[i].B_lu;
  }
  double get_line_f(int i) {
    return lines_[i].f_lu;
  }
  int get_line_l(int i) {
    return lines_[i].ll;
  }
  int get_line_u(int i) {
    return lines_[i].lu;
  }
  int get_line_bin(int i) {
    return lines_[i].bin;
  }
  int get_n_fuzz_lines() {
    return fuzz_lines_.n_lines;
  }
  double get_fuzz_line_nu(int i) {
    return fuzz_lines_.nu[i];
  }
  double get_fuzz_line_gf(int i) {
    return fuzz_lines_.gf[i];
  }
  double get_fuzz_line_El(int i) {
    return fuzz_lines_.El[i];
  }
  int get_fuzz_line_ion(int i) {
    return fuzz_lines_.ion[i];
  }
  int get_fuzz_line_bin(int i) {
    return fuzz_lines_.bin[i];
  }

};



class AtomicData
{

public:

  // constructor
  AtomicData();
  // Destructor
  ~AtomicData();

  // name of the data file that we read from
  std::string atom_datafile_;
  // version of the atomic data file
  int datafile_version_;

  // array of structures for the atomic data
  // of all possible atoms
  IndividualAtomData atomlist_[MAX_N_ATOMS];

  locate_array nu_grid_;

  int initialize(std::string, locate_array ng);
  int read_atomic_data(std::string fname, int z);
  int read_atomic_data(int z,int max_ion);
  int read_atomic_data(int z);
  int read_oldstyle_atomic_data(int z);
  int read_newstyle_atomic_data(int z);

  void print();
  void print_detailed(int);

  int read_fuzzfile_data(std::string fname);
  int read_fuzzfile_data_for_atom(std::string fname, int);

  std::string get_input_filename()
  {
    return atom_datafile_;
  }

  IndividualAtomData* get_pointer_to_individual_atom(int z)
  {
    return &(atomlist_[z]);
  }


};


#endif
