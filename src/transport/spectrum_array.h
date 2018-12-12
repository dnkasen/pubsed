#ifndef _SPECTRUM_H
#define _SPECTRUM_H 1

#include <string>
#include "particle.h"
#include "locate_array.h"
using std::string;

// default values
#define DEFAULT_NAME "spectrum"

class spectrum_array {
 
private:

  // MPI
  int nproc, my_rank;
  // spectrum name
  char name[1000];

   // bin arrays
  locate_array time_grid;
  locate_array wave_grid;
  locate_array mu_grid;
  locate_array phi_grid;

  // counting arrays
  std::vector<double> flux;
  std::vector<int>    click;

  // Indexing
  int n_elements;
  int a1, a2, a3;
  int index(int,int,int,int);
    
public:

  // constructors
  spectrum_array();
  
  // Initialize
  void init(std::vector<double>,std::vector<double>,int,int);
  void set_name(std::string);
  
  // Count a packets
  void count(double t, double w, double E, double *D);

  //  void normalize();
  void rescale(double);
  void wipe();

  // MPI functions
  void MPI_average();

  // Print out
  void print();

  void writeCheckpointSpectrum(std::string fname, std::string spectrum_name);
  void readCheckpointSpectrum(std::string fname, std::string n);

  bool is_equal(spectrum_array sa, bool complain);
};

#endif
