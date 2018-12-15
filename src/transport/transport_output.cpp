#include "hdf5.h"
#include "hdf5_hl.h"

#include <string.h>
#include <iostream>
#include <sstream>

#include "transport.h"
#include "physical_constants.h"

namespace pc = physical_constants;


//--------------------------------------------------------------
// Finalize the spectrum and output it
//--------------------------------------------------------------
void transport::output_spectrum(int it)
{

  std::stringstream ss;
  if (it < 0)  ss << "_final";
  else ss << "_" << it;
  string base = ss.str();

  string specname = params_->getScalar<string>("spectrum_name");
  if (specname != "")
    {
    optical_spectrum.set_name(specname + base);
    optical_spectrum.MPI_average();
    if (verbose) optical_spectrum.print();
  }

  string gamname = params_->getScalar<string>("gamma_name");
  if (gamname != "")
  {
    gamma_spectrum.set_name(gamname + base);
    gamma_spectrum.MPI_average();
    if (verbose) gamma_spectrum.print();
  }


}


void transport::wipe_spectra()
{

  gamma_spectrum.wipe();
  optical_spectrum.wipe();

}

void transport::output_spectrum()
{
  output_spectrum(0);
}

void transport::write_levels_to_plotfile(int iw)
{
  char pltfile[1000];
  sprintf(pltfile,"plt_%05d.h5",iw);

  const int RANK=1;


  hid_t file_id_dest = H5Fopen(pltfile,H5F_ACC_RDWR,H5P_DEFAULT);
  hid_t zone_id_dest = H5Gopen(file_id_dest,"zonedata",H5P_DEFAULT);
  for(int i=0;i<grid->n_zones;i++)
  {
    char zonefile[1000];
    sprintf(zonefile,"zone_%d.h5",i);

    hid_t file_id_src = H5Fopen(zonefile,H5F_ACC_RDONLY,H5P_DEFAULT);

    char zid[100];
    sprintf(zid,"%d",i);
    //hid_t zone_grp_src = H5Gopen(file_id_src,zid,H5P_DEFAULT);
    hid_t zone_grp_dest = H5Gopen(zone_id_dest,zid,H5P_DEFAULT);


    for(size_t j=0;j<gas_state_.atoms.size();j++)
    {
      char ag[100];
      int this_Z = gas_state_.elem_Z[j];
      sprintf(ag,"Z_%d",this_Z);
      hid_t atom_id_src = H5Gopen(file_id_src,ag,H5P_DEFAULT);
      hid_t atom_id_dest = H5Gcreate1(zone_grp_dest,ag,0);

      float* tmp_ion = new float[this_Z+1];
      hsize_t dims_ion[RANK]={(hsize_t)(this_Z+1)};
      
      H5LTread_dataset_float(atom_id_src,"ion_fraction",tmp_ion);
      H5LTmake_dataset(atom_id_dest,"ion_fraction",RANK,dims_ion,H5T_NATIVE_FLOAT,tmp_ion);

      int this_nl = gas_state_.atoms[j].n_levels_;
      float* tmp_level = new float[this_nl];
      hsize_t dims_level[RANK] = {(hsize_t)this_nl};
      H5LTread_dataset_float(atom_id_src,"level_fraction",tmp_level);
      H5LTmake_dataset(atom_id_dest,"level_fraction",RANK,dims_level,H5T_NATIVE_FLOAT,tmp_level);

      H5LTread_dataset_float(atom_id_src,"level_departure",tmp_level);
      H5LTmake_dataset(atom_id_dest,"level_departure",RANK,dims_level,H5T_NATIVE_FLOAT,tmp_level);

      H5Gclose(atom_id_dest);
      H5Gclose(atom_id_src);
      delete[] tmp_level;
      delete[] tmp_ion;
    }

    H5Gclose(zone_grp_dest);
    H5Fclose(file_id_src);

  }

  H5Gclose(zone_id_dest);
  H5Fclose(file_id_dest);
}


//------------------------------------------------------------
// Write detailed opacities and radiation field
// vs wavelength to an hdf file
// assumes that the pltfile has already been
// created
//------------------------------------------------------------
void transport::write_radiation_file(int iw)
{
  // get file name
  char zonefile[1000];
  sprintf(zonefile,"plt_%05d.h5",iw);

  // open hdf5 file
  hid_t file_id = H5Fopen( zonefile, H5F_ACC_RDWR, H5P_DEFAULT);
  const int RANK = 1;

  int n_nu = nu_grid.size();
  float* tmp_array = new float[n_nu];
  hsize_t  dims[RANK]={(hsize_t)n_nu};

  for (int j=0;j<n_nu;j++) tmp_array[j] = nu_grid.center(j);
  H5LTmake_dataset(file_id,"nu",RANK,dims,H5T_NATIVE_FLOAT,tmp_array);

  // write out mean opacities
  float* tz_array = new float[grid->n_zones];
  hsize_t  dims_z[RANK]={(hsize_t)grid->n_zones};
  for (int j=0;j<grid->n_zones;j++) tz_array[j] = planck_mean_opacity_[j];
  H5LTmake_dataset(file_id,"planck_mean",RANK,dims_z,H5T_NATIVE_FLOAT,tz_array);
  for (int j=0;j<grid->n_zones;j++) tz_array[j] = rosseland_mean_opacity_[j];
  H5LTmake_dataset(file_id,"rosseland_mean",RANK,dims_z,H5T_NATIVE_FLOAT,tz_array);
  delete[] tz_array;

  hid_t zone_dir = H5Gcreate1( file_id, "zonedata", 0 );

  // loop over zones for wavelength dependence opacities
  for (int i = 0; i < grid->n_zones; i++)
  {
    char zfile[100];
    sprintf(zfile,"%d",i);
    hid_t zone_id =  H5Gcreate1( zone_dir, zfile, 0 );

    // write total opacity
    if (omit_scattering_)
      for (int j=0;j<n_nu;j++)
        tmp_array[j] = (abs_opacity_[i][j])/grid->z[i].rho;
    else
      for (int j=0;j<n_nu;j++)
        tmp_array[j] = (scat_opacity_[i][j] + abs_opacity_[i][j])/grid->z[i].rho;
    H5LTmake_dataset(zone_id,"opacity",RANK,dims,H5T_NATIVE_FLOAT,tmp_array);

    // write absorption fraction
    for (int j=0;j<n_nu;j++)
    {
      double eps = 1;
      if (!omit_scattering_)
      {
        double topac = scat_opacity_[i][j] + abs_opacity_[i][j];
        if (topac == 0) eps = 1;
        else eps = abs_opacity_[i][j]/topac;
      }
      tmp_array[j] = eps;
    }
    H5LTmake_dataset(zone_id,"epsilon",RANK,dims,H5T_NATIVE_FLOAT,tmp_array);

    // write emissivity
    for (int j=0;j<n_nu;j++)  tmp_array[j] = emissivity_[i].get_value(j)/nu_grid.delta(j);
    H5LTmake_dataset(zone_id,"emissivity",RANK,dims,H5T_NATIVE_FLOAT,tmp_array);

    // write radiation field J
    for (int j=0;j<n_nu;j++)  {
      if (store_Jnu_) tmp_array[j] = J_nu_[i][j];
      else tmp_array[j] = 0; }
    H5LTmake_dataset(zone_id,"Jnu",RANK,dims,H5T_NATIVE_FLOAT,tmp_array);

    // if (write_levels)
    // {
    //   // just recalculate state for now... I know...
    //   // set up the state of the gas in this zone
    //   gas_state_.dens_ = grid->z[i].rho;
    //   gas_state_.temp_ = grid->z[i].T_gas;
    //   gas_state_.time_ = grid->t_now;
    //   gas_state_.set_mass_fractions(grid->z[i].X_gas);
    //   // solve for the state
    //   if (!gas_state_.grey_opacity_) gas_state_.solve_state(J_nu_[i]);

    //   for (size_t j=0;j<gas_state_.atoms.size();j++)
    //   {
    //     char afile[100];
    //     int this_Z = gas_state_.elem_Z[j];
    //     sprintf(afile,"Z_%d",this_Z);
    //     hid_t atom_id =  H5Gcreate1( zone_id, afile, 0 );

    //     // write out ionization fractions for this atom
    //     float tmp_ion[100];
    //     hsize_t  dims_ion[RANK]={(hsize_t)this_Z+1};
    //     for (int k=0;k<gas_state_.elem_Z[j]+1;k++)
    //       tmp_ion[k] = gas_state_.get_ionization_fraction(j,k);
    //     H5LTmake_dataset(atom_id,"ion_fraction",RANK,dims_ion,H5T_NATIVE_FLOAT,tmp_ion);

    //     // write out level populations for this atom
    //     int this_nl = gas_state_.atoms[j].n_levels_;
    //     float* tmp_level = new float[this_nl];
    //     hsize_t  dims_level[RANK]={(hsize_t)this_nl};
    //     for (int k=0;k<this_nl;k++)
    //       tmp_level[k] = gas_state_.get_level_fraction(j,k);
    //     H5LTmake_dataset(atom_id,"level_fraction",RANK,dims_level,H5T_NATIVE_FLOAT,tmp_level);

    //     // write out level departures for this atom
    //     for (int k=0;k<this_nl;k++)
    //       tmp_level[k] = gas_state_.get_level_departure(j,k);
    //     H5LTmake_dataset(atom_id,"level_departure",RANK,dims_level,H5T_NATIVE_FLOAT,tmp_level);

    //     H5Gclose(atom_id);
    //     delete[] tmp_level;
    //   }
    // }
    H5Gclose(zone_id);
  }
  H5Gclose(zone_dir);


  H5Fclose (file_id);
  delete[] tmp_array;
}
