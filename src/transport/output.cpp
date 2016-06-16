#include "hdf5.h"
#include "hdf5_hl.h"
#include <string.h>
#include <iostream>
#include <sstream>

#include "transport.h"

//--------------------------------------------------------------
// Finalize the spectrum and output it
//--------------------------------------------------------------
void transport::output_spectrum(int it)
{
  std::stringstream ss;
  if (it < 0)  ss << "_0" << it;
  else ss << "_" << it;
  string base = ss.str();
  if (!this->steady_state) base = "";

  
  string specname = params_->getScalar<string>("spectrum_name");
  if (specname != "") 
    {
    optical_spectrum.set_name(specname + base + ".dat");
    optical_spectrum.MPI_average();
    if (verbose) optical_spectrum.print();
  }
  optical_spectrum.wipe();

  string gamname = params_->getScalar<string>("gamma_name");
  if (gamname != "") 
  {
    gamma_spectrum.set_name(gamname + base + ".dat");
    gamma_spectrum.MPI_average();
    if (verbose) gamma_spectrum.print();
  }
  gamma_spectrum.wipe();

}


void transport::output_spectrum()
{
  output_spectrum(0);
}



//------------------------------------------------------------
// Write detailed opacities and radiation field
// vs wavelength to an hdf file
//------------------------------------------------------------
void transport::write_opacities(int iw)
{
  // get file name
  char zonefile[1000];
  sprintf(zonefile,"grid_%05d.h5",iw);
  sprintf(zonefile,"grid.h5");  

  // open hdf5 file
  hid_t file_id = H5Fcreate( zonefile, H5F_ACC_TRUNC, H5P_DEFAULT,  H5P_DEFAULT);
  
  const int RANK = 1;
  hsize_t  dims_z[RANK]={grid->n_zones};
  float* zone_arr = new float[grid->n_zones];

  // print out zone coordinates
  double r[3];
  for (int i=0;i<grid->n_zones;i++) {
    grid->coordinates(i,r);  
    zone_arr[i] = r[0];  }
  H5LTmake_dataset(file_id,"x",RANK,dims_z,H5T_NATIVE_FLOAT,zone_arr);
  for (int i=0;i<grid->n_zones;i++) {
    grid->coordinates(i,r);  
    zone_arr[i] = r[1];  }
  H5LTmake_dataset(file_id,"y",RANK,dims_z,H5T_NATIVE_FLOAT,zone_arr);
  for (int i=0;i<grid->n_zones;i++) {
    grid->coordinates(i,r);  
    zone_arr[i] = r[2];  }
  H5LTmake_dataset(file_id,"z",RANK,dims_z,H5T_NATIVE_FLOAT,zone_arr);

  // print out zone scalars
  for (int i=0;i<grid->n_zones;i++)  zone_arr[i] = grid->z[i].rho;
  H5LTmake_dataset(file_id,"rho",RANK,dims_z,H5T_NATIVE_FLOAT,zone_arr);
  for (int i=0;i<grid->n_zones;i++)  zone_arr[i] = grid->z[i].T_gas;
  H5LTmake_dataset(file_id,"T_gas",RANK,dims_z,H5T_NATIVE_FLOAT,zone_arr);

  int n_nu = nu_grid.size();
  float* tmp_array = new float[n_nu];
  hsize_t  dims[RANK]={n_nu};

  for (int j=0;j<n_nu;j++) tmp_array[j] = nu_grid.center(j);
  H5LTmake_dataset(file_id,"nu",RANK,dims,H5T_NATIVE_FLOAT,tmp_array);

  // loop over zones for wavelength dependence opacities
  for (int i = 0; i < grid->n_zones; i++)
  {
    char zfile[100];
    sprintf(zfile,"zone_%d",i);
    hid_t zone_id =  H5Gcreate1( file_id, zfile, 0 );

    // write total opacity
    for (int j=0;j<n_nu;j++)
      tmp_array[j] = (scat_opacity_[i][j] + abs_opacity_[i][j])/grid->z[i].rho;
    H5LTmake_dataset(zone_id,"opacity",RANK,dims,H5T_NATIVE_FLOAT,tmp_array);

    // write absorption fraction
    for (int j=0;j<n_nu;j++) {
      tmp_array[j] = abs_opacity_[i][j]/(scat_opacity_[i][j] + abs_opacity_[i][j]);
      if (scat_opacity_[i][j] + abs_opacity_[i][j] == 0) tmp_array[j] = 0; }
    H5LTmake_dataset(zone_id,"epsilon",RANK,dims,H5T_NATIVE_FLOAT,tmp_array);

    // write emissivity 
    for (int j=0;j<n_nu;j++)  tmp_array[j] = emissivity_[i].get_value(j);
    H5LTmake_dataset(zone_id,"emissivity",RANK,dims,H5T_NATIVE_FLOAT,tmp_array);

    // write radiation field J
    for (int j=0;j<n_nu;j++)  tmp_array[j] = J_nu_[i][j];
    H5LTmake_dataset(zone_id,"Jnu",RANK,dims,H5T_NATIVE_FLOAT,tmp_array); 

    H5Gclose(zone_id);
  }
  
  H5Fclose (file_id);
  delete[] tmp_array;
  delete[] zone_arr;
 
}


//------------------------------------------------------------
// Write detailed levels of atoms to hdf5 file
//------------------------------------------------------------
void transport::write_levels(int iw)
{
  // get file name
  char zonefile[1000];
  sprintf(zonefile,"levels_%05d.h5",iw);

  // open hdf5 file
  hid_t file_id = H5Fcreate( zonefile, H5F_ACC_TRUNC, H5P_DEFAULT,  H5P_DEFAULT);
  
  const int RANK = 1;
  hsize_t  dims_z[RANK]={grid->n_zones};
  float* zone_arr = new float[grid->n_zones];

  // print out zone coordinates
  double r[3];
  for (int i=0;i<grid->n_zones;i++) {
    grid->coordinates(i,r);  
    zone_arr[i] = r[0];  }
  H5LTmake_dataset(file_id,"x",RANK,dims_z,H5T_NATIVE_FLOAT,zone_arr);
  for (int i=0;i<grid->n_zones;i++) {
    grid->coordinates(i,r);  
    zone_arr[i] = r[1];  }
  H5LTmake_dataset(file_id,"y",RANK,dims_z,H5T_NATIVE_FLOAT,zone_arr);
  for (int i=0;i<grid->n_zones;i++) {
    grid->coordinates(i,r);  
    zone_arr[i] = r[2];  }
  H5LTmake_dataset(file_id,"z",RANK,dims_z,H5T_NATIVE_FLOAT,zone_arr);

  // print out zone scalars
  for (int i=0;i<grid->n_zones;i++)  zone_arr[i] = grid->z[i].rho;
  H5LTmake_dataset(file_id,"rho",RANK,dims_z,H5T_NATIVE_FLOAT,zone_arr);
  for (int i=0;i<grid->n_zones;i++)  zone_arr[i] = grid->z[i].T_gas;
  H5LTmake_dataset(file_id,"T_gas",RANK,dims_z,H5T_NATIVE_FLOAT,zone_arr);

  int n_nu = nu_grid.size();
  float* tmp_array = new float[n_nu];
  hsize_t  dims[RANK]={n_nu};

  for (int j=0;j<n_nu;j++) tmp_array[j] = nu_grid.center(j);
  H5LTmake_dataset(file_id,"nu",RANK,dims,H5T_NATIVE_FLOAT,tmp_array);
	
  // loop over zones for wavelength dependence opacities
  for (int i = 0; i < grid->n_zones; i++)
  {
    // just recalculate state for now... I know...
    // set up the state of the gas in this zone
    gas.dens = grid->z[i].rho;
    gas.temp = grid->z[i].T_gas;
    gas.time = grid->t_now;
    gas.set_mass_fractions(grid->z[i].X_gas);
    // solve for the state 
    if (!gas.grey_opacity_) gas.solve_state(J_nu_[i]);

    char zfile[100];
    sprintf(zfile,"zone_%d",i);
    hid_t zone_id =  H5Gcreate1( file_id, zfile, 0 );

    for (int j=0;j<gas.atoms.size();j++)
    {
    	char afile[100];
    	int this_Z = gas.elem_Z[j];
	    sprintf(afile,"Z_%d",this_Z);
	    hid_t atom_id =  H5Gcreate1( zone_id, afile, 0 );

	    // write out ionization fractions for this atom
	   	float tmp_ion[100];
	   	hsize_t  dims_ion[RANK]={this_Z+1};
	   	for (int k=0;k<gas.elem_Z[j]+1;k++)
	   		tmp_ion[k] = gas.get_ionization_fraction(j,k);
	   	H5LTmake_dataset(atom_id,"ion_fraction",RANK,dims_ion,H5T_NATIVE_FLOAT,tmp_ion);

	   	// write out level populations for this atom
	   	int this_nl = gas.atoms[j].n_levels;
	  	float* tmp_level = new float[this_nl];
	   	hsize_t  dims_level[RANK]={this_nl};
	   	for (int k=0;k<this_nl;k++)
	   		tmp_level[k] = gas.get_level_fraction(j,k);
	   	H5LTmake_dataset(atom_id,"level_fraction",RANK,dims_level,H5T_NATIVE_FLOAT,tmp_level);

	   	// write out level departures for this atom
	   	for (int k=0;k<this_nl;k++)
	   		tmp_level[k] = gas.get_level_departure(j,k);
	   	H5LTmake_dataset(atom_id,"level_departure",RANK,dims_level,H5T_NATIVE_FLOAT,tmp_level);

      H5Gclose(atom_id);
  		delete[] tmp_level;
    }

  }
  
  H5Fclose (file_id);
  delete[] tmp_array;
  delete[] zone_arr;
  
}
