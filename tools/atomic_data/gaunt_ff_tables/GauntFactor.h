#ifndef _GAUNT_FACTOR_H
#define _GAUNT_FACTOR_H 1

#include <string>
#include <vector>
#include "locate_array.h"
#include "hdf5.h"
#include "hdf5_hl.h"


class GauntFactor
{

public:



private:

void read_data()
{

  // get number of atoms
  int max_atoms;
  sprintf(dset,"max_atoms");
  status = H5LTread_dataset_int(file_id, dset,&max_atoms);
  if (status != 0) return -1;

  // read level excitation energy
  double *E_darr = new double[tot_n_levels];
  sprintf(dset,"%s%s",ionname,"level_E");
  status = H5LTread_dataset_double(file_id,dset,E_darr);
  if (status != 0) return -1;
}


};
