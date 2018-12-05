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


//------------------------------------------------------------
// Write detailed opacities and radiation field
// vs wavelength to an hdf file
// assumes that the pltfile has already been
// created
//------------------------------------------------------------
void transport::write_radiation_file(int iw, int write_levels)
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

    if (write_levels)
    {
      // just recalculate state for now... I know...
      // set up the state of the gas in this zone
      gas_state_.dens_ = grid->z[i].rho;
      gas_state_.temp_ = grid->z[i].T_gas;
      gas_state_.time_ = grid->t_now;
      gas_state_.set_mass_fractions(grid->z[i].X_gas);
      // solve for the state
      if (!gas_state_.grey_opacity_) gas_state_.solve_state(J_nu_[i]);

      for (size_t j=0;j<gas_state_.atoms.size();j++)
      {
        char afile[100];
        int this_Z = gas_state_.elem_Z[j];
        sprintf(afile,"Z_%d",this_Z);
        hid_t atom_id =  H5Gcreate1( zone_id, afile, 0 );

        // write out ionization fractions for this atom
        float tmp_ion[100];
        hsize_t  dims_ion[RANK]={(hsize_t)this_Z+1};
        for (int k=0;k<gas_state_.elem_Z[j]+1;k++)
          tmp_ion[k] = gas_state_.get_ionization_fraction(j,k);
        H5LTmake_dataset(atom_id,"ion_fraction",RANK,dims_ion,H5T_NATIVE_FLOAT,tmp_ion);

        // write out level populations for this atom
        int this_nl = gas_state_.atoms[j].n_levels_;
        float* tmp_level = new float[this_nl];
        hsize_t  dims_level[RANK]={(hsize_t)this_nl};
        for (int k=0;k<this_nl;k++)
          tmp_level[k] = gas_state_.get_level_fraction(j,k);
        H5LTmake_dataset(atom_id,"level_fraction",RANK,dims_level,H5T_NATIVE_FLOAT,tmp_level);

        // write out level departures for this atom
        for (int k=0;k<this_nl;k++)
          tmp_level[k] = gas_state_.get_level_departure(j,k);
        H5LTmake_dataset(atom_id,"level_departure",RANK,dims_level,H5T_NATIVE_FLOAT,tmp_level);

        H5Gclose(atom_id);
        delete[] tmp_level;
      }
    }
    H5Gclose(zone_id);
  }
  H5Gclose(zone_dir);


  H5Fclose (file_id);
  delete[] tmp_array;
}

void transport::checkpoint_particles(char* fname) {
  // Figures out what every rank's offset is going to be in the big particle list
  int my_n_particles = n_particles();
  int my_offset, global_n_particles_total;
  hsize_t* global_n_particles = new hsize_t[MPI_nprocs];
  hsize_t* particle_offsets = new hsize_t[MPI_nprocs];
  MPI_Gather(&my_n_particles, 1, MPI_INT, global_n_particles, MPI_nprocs, MPI_INT, 0, MPI_COMM_WORLD);
  if (MPI_myID == 0) {
    for (int i = 0; i < MPI_nprocs; i++) {
      particle_offsets[i] = global_n_particles_total;
      global_n_particles_total += global_n_particles[i];
    }
  }
  // Rank 0 tells all of the other ranks what their offsets will be
  MPI_Scatter(particle_offsets, MPI_nprocs, MPI_INT, &my_offset, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&global_n_particles_total, 1, MPI_INT, 0, MPI_COMM_WORLD);
  // Rank 0 sets up all of the datasets in the H5 file
  if (MPI_myID == 0) {
    int ndim1 = 1;
    hsize_t dims1[1] =  {global_n_particles_total};
    int ndim3 = 2;
    hsize_t dims3[2] = {global_n_particles_total, 3};
    createGroup(fname, "particles");
    createDataset(fname, "particles", "type", ndim1, dims1, H5T_NATIVE_INT);
    createDataset(fname, "particles", "x", ndim3, dims3, H5T_NATIVE_FLOAT);
    createDataset(fname, "particles", "D", ndim3, dims3, H5T_NATIVE_FLOAT);
    createDataset(fname, "particles", "ind", ndim1, dims1, H5T_NATIVE_FLOAT);
    createDataset(fname, "particles", "t", ndim1, dims1, H5T_NATIVE_FLOAT);
    createDataset(fname, "particles", "e", ndim1, dims1, H5T_NATIVE_FLOAT);
    createDataset(fname, "particles", "nu", ndim1, dims1, H5T_NATIVE_FLOAT);
    createDataset(fname, "particles", "gamma", ndim1, dims1, H5T_NATIVE_FLOAT);
    createDataset(fname, "particles", "dshift", ndim1, dims1, H5T_NATIVE_FLOAT);
    createDataset(fname, "particles", "dvds", ndim1, dims1, H5T_NATIVE_FLOAT);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  for (int i = 0; i < MPI_nprocs; i++) {
    if (i == MPI_myID) {
      writeParticleProp(fname, "particles", "type", global_n_particles_total, my_offset);
      writeParticleProp(fname, "particles", "x", global_n_particles_total, my_offset);
      writeParticleProp(fname, "particles", "D", global_n_particles_total, my_offset);
      writeParticleProp(fname, "particles", "ind", global_n_particles_total, my_offset);
      writeParticleProp(fname, "particles", "t", global_n_particles_total, my_offset);
      writeParticleProp(fname, "particles", "e", global_n_particles_total, my_offset);
      writeParticleProp(fname, "particles", "nu", global_n_particles_total, my_offset);
      writeParticleProp(fname, "particles", "gamma", global_n_particles_total, my_offset);
      writeParticleProp(fname, "particles", "dshift", global_n_particles_total, my_offset);
      writeParticleProp(fname, "particles", "dvds", global_n_particles_total, my_offset);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  delete[] global_n_particles;
  delete[] particle_offsets;
}

// Writes out particle data, assuming that the particles group already exists in
// the hdf5 file named file. The 
void transport::writeParticleProp(char* fname, std::string fieldname, int total_particles, int offset) {
  int n_dims = 1;
  int n_particles_local = n_particles();
  int* buffer;
  if (fieldname == "type") {
    hid_t t = H5T_NATIVE_INT;
    int* buffer = new int[n_particles_local];
    for (int i = 0; i < n_particles_local; i++) {
      buffer[i] = particles[i].type;
    }
  }
  else if (fieldname == "x") {
    hid_t t = H5T_NATIVE_DOUBLE;
    n_dims = 2;
    double* buffer = new double[n_particles_local * 3];
    for (int i = 0; i < n_particles_local; i++) {
      buffer[i * 3] = particles[i].x[0];
      buffer[i * 3 + 1] = particles[i].x[1];
      buffer[i * 3 + 2] = particles[i].x[2];
    }
  }
  else if (fieldname == "D") {
    hid_t t = H5T_NATIVE_DOUBLE;
    n_dims = 2;
    double* buffer = new double[n_particles_local * 3];
    for (int i = 0; i < n_particles_local; i++) {
      buffer[i * 3] = particles[i].D[0];
      buffer[i * 3 + 1] = particles[i].D[1];
      buffer[i * 3 + 2] = particles[i].D[2];
    }
  }
  else if (fieldname == "ind") {
    hid_t t = H5T_NATIVE_INT;
    int* buffer = new int[n_particles_local];
    for (int i = 0; i < n_particles_local; i++) {
      buffer[i] = particles[i].ind;
    }
  }
  else if (fieldname == "t") {
    hid_t t = H5T_NATIVE_DOUBLE;
    double* buffer = new double[n_particles_local];
    for (int i = 0; i < n_particles_local; i++) {
      buffer[i] = particles[i].t;
    }
  }
  else if (fieldname == "e") {
    hid_t t = H5T_NATIVE_DOUBLE;
    double* buffer = new double[n_particles_local];
    for (int i = 0; i < n_particles_local; i++) {
      buffer[i] = particles[i].e;
    }
  }
  else if (fieldname == "nu") {
    hid_t t = H5T_NATIVE_DOUBLE;
    double* buffer = new double[n_particles_local];
    for (int i = 0; i < n_particles_local; i++) {
      buffer[i] = particles[i].nu;
    }
  }
  else if (fieldname == "gamma") {
    hid_t t = H5T_NATIVE_DOUBLE;
    double* buffer = new double[n_particles_local];
    for (int i = 0; i < n_particles_local; i++) {
      buffer[i] = particles[i].gamma;
    }
  }
  else if (fieldname == "dshift") {
    hid_t t = H5T_NATIVE_DOUBLE;
    double* buffer = new double[n_particles_local];
    for (int i = 0; i < n_particles_local; i++) {
      buffer[i] = particles[i].dshift;
    }
  }
  else if (fieldname == "dvds") {
    hid_t t = H5T_NATIVE_DOUBLE;
    double* buffer = new double[n_particles_local];
    for (int i = 0; i < n_particles_local; i++) {
      buffer[i] = particles[i].dvds;
    }
  }
  else {
    std::cerr << "Particle field " << fieldname << " does not exist. Terminating" << std::endl;
    exit(3);
  }

  if (n_dims == 1) {
    int[1] start = {offset};
    int[1] size = {n_particles_local};
    int[1] total_size = {total_particles};
  }
  else if (n_dims == 2) {
    int[2] start = {offset, 0};
    int[2] size = {n_particles_local, 3};
    int[2] total_size = {total_particles, 3};
  }
  else {
    std::cerr << "Dimension count " << n_dims << " not allowed." <<std::endl;
    exit(3);
  }

  writePatch(fname, "particles", fieldname, buffer, t, n_dims, start, size, total_size);
  if (buffer) {
    delete[] buffer;
  }
}
