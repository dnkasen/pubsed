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
  int suppress_txt = params_->getScalar<int>("spectrum_suppress_txt");
  if (specname != "")
    {
    optical_spectrum.set_name(specname + base);
    optical_spectrum.MPI_average();
    if (verbose) optical_spectrum.print(suppress_txt);
  }

  string gamname = params_->getScalar<string>("gamma_name");
  if (gamname != "")
  {
    gamma_spectrum.set_name(gamname + base);
    gamma_spectrum.MPI_average();
    if (verbose) gamma_spectrum.print(suppress_txt);
  }

  string partname = params_->getScalar<string>("spectrum_particle_list_name");
  if (partname != "")
  {
    if (verbose) {
      std::cout << "# writing escaped particle list" << std::endl;
      createFile(partname + base + ".h5");
    }
    // whenever escaped particles are written out, clear the list so we are less likely to run
    // into memory issues
    writeCheckpointParticles(particles_escaped, partname + base + ".h5", "particles_escaped");
    clearEscapedParticles();
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


    GasState* gas_state = &(gas_state_vec_[0]);
    for(size_t j=0;j<gas_state->atoms.size();j++)
    {
      char ag[100];
      int this_Z = gas_state->elem_Z[j];
      sprintf(ag,"Z_%d",this_Z);
      hid_t atom_id_src = H5Gopen(file_id_src,ag,H5P_DEFAULT);
      hid_t atom_id_dest = H5Gcreate1(zone_grp_dest,ag,0);

      double* tmp_ion = new double[this_Z+1];
      hsize_t dims_ion[RANK]={(hsize_t)(this_Z+1)};

      H5LTread_dataset_double(atom_id_src,"ion_fraction",tmp_ion);
      H5LTmake_dataset(atom_id_dest,"ion_fraction",RANK,dims_ion,H5T_NATIVE_DOUBLE,tmp_ion);

      int this_nl = gas_state->atoms[j].n_levels_;
      double* tmp_level = new double[this_nl];
      hsize_t dims_level[RANK] = {(hsize_t)this_nl};
      H5LTread_dataset_double(atom_id_src,"level_fraction",tmp_level);
      H5LTmake_dataset(atom_id_dest,"level_fraction",RANK,dims_level,H5T_NATIVE_DOUBLE,tmp_level);

      H5LTread_dataset_double(atom_id_src,"level_departure",tmp_level);
      H5LTmake_dataset(atom_id_dest,"level_departure",RANK,dims_level,H5T_NATIVE_DOUBLE,tmp_level);

      int this_nd = 1;
      double* tmp_ndens = new double[this_nd];
      hsize_t dims_ndens[RANK] = {(hsize_t)this_nd};
      H5LTread_dataset_double(atom_id_src,"n_dens",tmp_ndens);
      H5LTmake_dataset(atom_id_dest,"n_dens",RANK,dims_ndens,H5T_NATIVE_DOUBLE,tmp_ndens);

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

  int n_nu = nu_grid_.size();
  float* tmp_array = new float[n_nu];
  hsize_t  dims[RANK]={(hsize_t)n_nu};

  for (int j=0;j<n_nu;j++) tmp_array[j] = nu_grid_.center(j);
  H5LTmake_dataset(file_id,"nu",RANK,dims,H5T_NATIVE_FLOAT,tmp_array);

  // write out mean opacities
  float* tz_array = new float[grid->n_zones];
  hsize_t  dims_z[RANK]={(hsize_t)grid->n_zones};

  for (int j=0;j<grid->n_zones;j++) tz_array[j] = planck_mean_opacity_[j];
  H5LTmake_dataset(file_id,"planck_mean",RANK,dims_z,H5T_NATIVE_FLOAT,tz_array);


  for (int j=0;j<grid->n_zones;j++) tz_array[j] = rosseland_mean_opacity_[j];
  H5LTmake_dataset(file_id,"rosseland_mean",RANK,dims_z,H5T_NATIVE_FLOAT,tz_array);

  if (params_->getScalar<int>("opacity_use_nlte"))
  {
    for (int j=0;j<grid->n_zones;j++) tz_array[j] = bf_heating[j];
    H5LTmake_dataset(file_id,"bf_heating",RANK,dims_z,H5T_NATIVE_FLOAT,tz_array);
    for (int j=0;j<grid->n_zones;j++) tz_array[j] = bf_cooling[j];
    H5LTmake_dataset(file_id,"bf_cooling",RANK,dims_z,H5T_NATIVE_FLOAT,tz_array);
    for (int j=0;j<grid->n_zones;j++) tz_array[j] = ff_heating[j];
    H5LTmake_dataset(file_id,"ff_heating",RANK,dims_z,H5T_NATIVE_FLOAT,tz_array);
    for (int j=0;j<grid->n_zones;j++) tz_array[j] = ff_cooling[j];
    H5LTmake_dataset(file_id,"ff_cooling",RANK,dims_z,H5T_NATIVE_FLOAT,tz_array);
    for (int j=0;j<grid->n_zones;j++) tz_array[j] = coll_cooling[j];
    H5LTmake_dataset(file_id,"coll_cooling",RANK,dims_z,H5T_NATIVE_FLOAT,tz_array);
      //      for (int j=0;j<grid->n_zones;j++) tz_array[j] = grid->z[j].e_abs_compton;
      //      H5LTmake_dataset(file_id,"compton_heating_rate",RANK,dims_z,H5T_NATIVE_FLOAT,tz_array);

  }

  delete[] tz_array;


  if (use_ddmc_)
  {
    int* tzi_array = new int[grid->n_zones];
    for (int j=0;j<grid->n_zones;j++) tzi_array[j] = ddmc_use_in_zone_[j];
    H5LTmake_dataset(file_id,"use_ddmc",RANK,dims_z,H5T_NATIVE_INT,tzi_array);
    delete[] tzi_array;
  }

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
    for (int j=0;j<n_nu;j++)  tmp_array[j] = emissivity_[i].get_value(j)/nu_grid_.delta(j);
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
    //   gas_state->dens_ = grid->z[i].rho;
    //   gas_state->temp_ = grid->z[i].T_gas;
    //   gas_state->time_ = grid->t_now;
    //   gas_state->set_mass_fractions(grid->z[i].X_gas);
    //   // solve for the state
    //   if (!gas_state->grey_opacity_) gas_state->solve_state(J_nu_[i]);

    //   for (size_t j=0;j<gas_state->atoms.size();j++)
    //   {
    //     char afile[100];
    //     int this_Z = gas_state->elem_Z[j];
    //     sprintf(afile,"Z_%d",this_Z);
    //     hid_t atom_id =  H5Gcreate1( zone_id, afile, 0 );

    //     // write out ionization fractions for this atom
    //     float tmp_ion[100];
    //     hsize_t  dims_ion[RANK]={(hsize_t)this_Z+1};
    //     for (int k=0;k<gas_state->elem_Z[j]+1;k++)
    //       tmp_ion[k] = gas_state->get_ionization_fraction(j,k);
    //     H5LTmake_dataset(atom_id,"ion_fraction",RANK,dims_ion,H5T_NATIVE_FLOAT,tmp_ion);

    //     // write out level populations for this atom
    //     int this_nl = gas_state->atoms[j].n_levels_;
    //     float* tmp_level = new float[this_nl];
    //     hsize_t  dims_level[RANK]={(hsize_t)this_nl};
    //     for (int k=0;k<this_nl;k++)
    //       tmp_level[k] = gas_state->get_level_fraction(j,k);
    //     H5LTmake_dataset(atom_id,"level_fraction",RANK,dims_level,H5T_NATIVE_FLOAT,tmp_level);

    //     // write out level departures for this atom
    //     for (int k=0;k<this_nl;k++)
    //       tmp_level[k] = gas_state->get_level_departure(j,k);
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


void transport::writeCheckpointParticlesAll(std::string fname) {
  writeCheckpointParticles(particles, fname, "particles");
  writeCheckpointParticles(particles_escaped, fname, "particles_escaped");
}

void transport::writeCheckpointParticles(std::vector<particle>& particle_list, 
    std::string fname, std::string groupname) {
  // Figures out what every rank's offset is going to be in the big particle list
  int my_n_particles = particle_list.size();
  int my_offset;
  int global_n_particles_total = 0;
  int* global_n_particles = new int[MPI_nprocs];
  int* particle_offsets = new int[MPI_nprocs];
  MPI_Gather(&my_n_particles, 1, MPI_INT, global_n_particles, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (MPI_myID == 0) {
    for (int i = 0; i < MPI_nprocs; i++) {
      particle_offsets[i] = global_n_particles_total;
      global_n_particles_total += global_n_particles[i];
    }
  }
  // Rank 0 tells all of the other ranks what their offsets will be
  MPI_Scatter(particle_offsets, 1, MPI_INT, &my_offset, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&global_n_particles_total, 1, MPI_INT, 0, MPI_COMM_WORLD);
  // Rank 0 sets up all of the datasets in the H5 file
  if (MPI_myID == 0) {
    //createFile(fname);
    int ndim1 = 1;
    hsize_t dims1[1] =  {global_n_particles_total};
    hsize_t nranks_dim[1] = {MPI_nprocs};
    int ndim3 = 2;
    hsize_t dims3[2] = {global_n_particles_total, 3};
    createGroup(fname, groupname);
    createDataset(fname, groupname, "counts_by_rank", ndim1, nranks_dim, H5T_NATIVE_INT);
    writeSimple(fname, groupname, "counts_by_rank", global_n_particles, H5T_NATIVE_INT);

    createDataset(fname, groupname, "type", ndim1, dims1, H5T_NATIVE_INT);
    createDataset(fname, groupname, "x", ndim3, dims3, H5T_NATIVE_DOUBLE);
    createDataset(fname, groupname, "D", ndim3, dims3, H5T_NATIVE_DOUBLE);
    createDataset(fname, groupname, "x_interact", ndim3, dims3, H5T_NATIVE_DOUBLE);
    createDataset(fname, groupname, "ind", ndim1, dims1, H5T_NATIVE_INT);
    createDataset(fname, groupname, "t", ndim1, dims1, H5T_NATIVE_DOUBLE);
    createDataset(fname, groupname, "e", ndim1, dims1, H5T_NATIVE_DOUBLE);
    createDataset(fname, groupname, "nu", ndim1, dims1, H5T_NATIVE_DOUBLE);
    createDataset(fname, groupname, "gamma", ndim1, dims1, H5T_NATIVE_DOUBLE);
    createDataset(fname, groupname, "dshift", ndim1, dims1, H5T_NATIVE_DOUBLE);
    createDataset(fname, groupname, "dvds", ndim1, dims1, H5T_NATIVE_DOUBLE);
    createDataset(fname, groupname, "fate", ndim1, dims1, H5T_NATIVE_INT);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  for (int i = 0; i < MPI_nprocs; i++) {
    if (i == MPI_myID) {
      writeParticleProp(fname, "type", groupname, particle_list, global_n_particles_total, my_offset);
      writeParticleProp(fname, "x", groupname, particle_list, global_n_particles_total, my_offset);
      writeParticleProp(fname, "D", groupname, particle_list, global_n_particles_total, my_offset);
      writeParticleProp(fname, "x_interact", groupname, particle_list, global_n_particles_total, my_offset);
      writeParticleProp(fname, "ind", groupname, particle_list, global_n_particles_total, my_offset);
      writeParticleProp(fname, "t", groupname, particle_list, global_n_particles_total, my_offset);
      writeParticleProp(fname, "e", groupname, particle_list, global_n_particles_total, my_offset);
      writeParticleProp(fname, "nu", groupname, particle_list, global_n_particles_total, my_offset);
      writeParticleProp(fname, "gamma", groupname, particle_list, global_n_particles_total, my_offset);
      writeParticleProp(fname, "dshift", groupname, particle_list, global_n_particles_total, my_offset);
      writeParticleProp(fname, "dvds", groupname, particle_list, global_n_particles_total, my_offset);
      writeParticleProp(fname, "fate", groupname, particle_list, global_n_particles_total, my_offset);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  delete[] global_n_particles;
  delete[] particle_offsets;
}

// Writes out particle data, assuming that the particles group already exists in
// the hdf5 file named file. The
void transport::writeParticleProp(std::string fname, std::string fieldname, 
    std::string groupname, std::vector<particle>& particle_list, int total_particles, int offset) {
  int n_dims = 1;
  int n_particles_local = particle_list.size();
  int* buffer_i;
  double* buffer_d;
  hid_t t = H5T_NATIVE_DOUBLE;
  if (fieldname == "type") {
    t = H5T_NATIVE_INT;
    buffer_i = new int[n_particles_local];
    for (int i = 0; i < n_particles_local; i++) {
      buffer_i[i] = particle_list[i].type;
    }
  }
  else if (fieldname == "x") {
    n_dims = 2;
    buffer_d = new double[n_particles_local * 3];
    for (int i = 0; i < n_particles_local; i++) {
      buffer_d[i * 3] = particle_list[i].x[0];
      buffer_d[i * 3 + 1] = particle_list[i].x[1];
      buffer_d[i * 3 + 2] = particle_list[i].x[2];
    }
  }
  else if (fieldname == "D") {
    n_dims = 2;
    buffer_d = new double[n_particles_local * 3];
    for (int i = 0; i < n_particles_local; i++) {
      buffer_d[i * 3] = particle_list[i].D[0];
      buffer_d[i * 3 + 1] = particle_list[i].D[1];
      buffer_d[i * 3 + 2] = particle_list[i].D[2];
    }
  }
  else if (fieldname == "x_interact") {
    n_dims = 2;
    buffer_d = new double[n_particles_local * 3];
    for (int i = 0; i < n_particles_local; i++) {
      buffer_d[i * 3] = particle_list[i].x_interact[0];
      buffer_d[i * 3 + 1] = particle_list[i].x_interact[1];
      buffer_d[i * 3 + 2] = particle_list[i].x_interact[2];
    }
  }
  else if (fieldname == "ind") {
    t = H5T_NATIVE_INT;
    buffer_i = new int[n_particles_local];
    for (int i = 0; i < n_particles_local; i++) {
      buffer_i[i] = particle_list[i].ind;
    }
  }
  else if (fieldname == "t") {
    buffer_d = new double[n_particles_local];
    for (int i = 0; i < n_particles_local; i++) {
      buffer_d[i] = particle_list[i].t;
    }
  }
  else if (fieldname == "e") {
    buffer_d = new double[n_particles_local];
    for (int i = 0; i < n_particles_local; i++) {
      buffer_d[i] = particle_list[i].e;
    }
  }
  else if (fieldname == "nu") {
    buffer_d = new double[n_particles_local];
    for (int i = 0; i < n_particles_local; i++) {
      buffer_d[i] = particle_list[i].nu;
    }
  }
  else if (fieldname == "gamma") {
    buffer_d = new double[n_particles_local];
    for (int i = 0; i < n_particles_local; i++) {
      buffer_d[i] = particle_list[i].gamma;
    }
  }
  else if (fieldname == "dshift") {
    buffer_d = new double[n_particles_local];
    for (int i = 0; i < n_particles_local; i++) {
      buffer_d[i] = particle_list[i].dshift;
    }
  }
  else if (fieldname == "dvds") {
    buffer_d = new double[n_particles_local];
    for (int i = 0; i < n_particles_local; i++) {
      buffer_d[i] = particle_list[i].dvds;
    }
  }
  else if (fieldname == "fate") {
    t = H5T_NATIVE_INT;
    buffer_i = new int[n_particles_local];
    for (int i = 0; i < n_particles_local; i++) {
      buffer_i[i] = particle_list[i].fate;
    }
  }
  else {
    std::cerr << "Particle field " << fieldname << " does not exist. Terminating" << std::endl;
    exit(3);
  }

  // Because the code will only look at n_dims dimensions worth of measurements, we can always
  // fill in the whole 2-long arrays. If n_dims == 1, the code will just ignore the second
  // entries.
  int start[2] = {offset, 0};
  int size[2] = {n_particles_local, 3};
  int total_size[2] = {total_particles, 3};

  if (n_dims != 1 && n_dims != 2) {
    std::cerr << "Dimension count " << n_dims << " not allowed." <<std::endl;
    exit(3);
  }
  if (t == H5T_NATIVE_INT) {
    writePatch(fname, groupname, fieldname.c_str(), buffer_i, t, n_dims, start, size, total_size);
    delete[] buffer_i;
  }
  else if (t == H5T_NATIVE_DOUBLE) {
    writePatch(fname, groupname, fieldname.c_str(), buffer_d, t, n_dims, start, size, total_size);
    delete[] buffer_d;
  }
  else {
    std::cerr << "HDF5 type is wrong" <<std::endl;
  }
}

void transport::writeCheckpointSpectra(std::string fname) {
  optical_spectrum.writeCheckpointSpectrum(fname, "optical spectrum");
  gamma_spectrum.writeCheckpointSpectrum(fname, "gamma spectrum");
}

void transport::writeCheckpointRNG(std::string fname) {
  rangen.writeCheckpointRNG(fname);
}

void transport::readCheckpointParticles(std::vector<particle>& particle_list,
    std::string fname, std::string groupname, bool test, bool all_one_rank) {
  /* Get number of particles that are stored in the file */
  hsize_t global_n_particles_total, n_ranks_old;
  int my_n_particles, my_offset;
  int* global_n_particles;
  int* particle_offsets;
  // all_one_rank option reads all particles into this rank and does not divide
  // them up
  if (!all_one_rank) {
    for (int rank = 0; rank < MPI_nprocs; rank++) {
      if (MPI_myID == rank) {
        getH5dims(fname, groupname, "e", &global_n_particles_total);
        getH5dims(fname, groupname, "counts_by_rank", &n_ranks_old);
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }

    global_n_particles = new int[MPI_nprocs];
    if (n_ranks_old != MPI_nprocs) {
      if (verbose) {
        std::cerr << "Rank count pre- and post-restart don't match. Distributing" <<
          " particles evenly." << std::endl;
        std::cerr << n_ranks_old << " " << MPI_nprocs << std::endl;
        std::cerr << "Exiting." << std::endl;
        exit(10);
      }
      /* Make sure each particle gets "claimed" by one rank */
      my_n_particles = floor(global_n_particles_total / (1.0 * MPI_nprocs));
      int remainder = global_n_particles_total % MPI_nprocs;
      if (MPI_myID < remainder) {
        my_n_particles += 1;
      }
    }
    else {
      for (int rank = 0; rank < MPI_nprocs; rank++) {
        if (MPI_myID == rank) {
          readSimple(fname, groupname, "counts_by_rank", global_n_particles, H5T_NATIVE_INT);
          my_n_particles = global_n_particles[MPI_myID];
        }
        MPI_Barrier(MPI_COMM_WORLD);
      }
    }
    /* Figure out what each rank's offset will be in the particle array */
    particle_offsets = new int[MPI_nprocs];
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(&my_n_particles, 1, MPI_INT, global_n_particles, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (MPI_myID == 0) {
      particle_offsets[0] = 0;
      for (int i = 1; i < MPI_nprocs; i++) {
        particle_offsets[i] = particle_offsets[i - 1] + global_n_particles[i - 1];
      }
    }
    // Rank 0 tells all of the other ranks what their offsets will be
    MPI_Scatter(particle_offsets, 1, MPI_INT, &my_offset, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
  }
  else {
    getH5dims(fname, groupname, "e", &global_n_particles_total);
    std::cerr << "particle count " << global_n_particles_total << std::endl;
    my_n_particles = global_n_particles_total;
    my_offset = 0;
  }
  std::cerr << my_n_particles << std::endl;
  particle_list.resize(my_n_particles);


  /* Read in all of the quantities */
  for (int i = 0; i < MPI_nprocs; i++) {
    if (i == MPI_myID) {
      readParticleProp(fname, "type", groupname, particle_list, global_n_particles_total, my_offset);
      readParticleProp(fname, "x", groupname, particle_list, global_n_particles_total, my_offset);
      readParticleProp(fname, "D", groupname, particle_list, global_n_particles_total, my_offset);
      readParticleProp(fname, "x_interact", groupname, particle_list, global_n_particles_total, my_offset);
      readParticleProp(fname, "ind", groupname, particle_list, global_n_particles_total, my_offset);
      readParticleProp(fname, "t", groupname, particle_list, global_n_particles_total, my_offset);
      readParticleProp(fname, "e", groupname, particle_list, global_n_particles_total, my_offset);
      readParticleProp(fname, "nu", groupname, particle_list, global_n_particles_total, my_offset);
      readParticleProp(fname, "gamma", groupname, particle_list, global_n_particles_total, my_offset);
      readParticleProp(fname, "dshift", groupname, particle_list, global_n_particles_total, my_offset);
      readParticleProp(fname, "dvds", groupname, particle_list, global_n_particles_total, my_offset);
      readParticleProp(fname, "fate", groupname, particle_list, global_n_particles_total, my_offset);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  if (!all_one_rank) {
    delete[] global_n_particles;
    delete[] particle_offsets;
  }
}

void transport::readParticleProp(std::string fname, std::string fieldname,
    std::string groupname, std::vector<particle>& particle_list, int total_particles, int offset) {
  int n_dims = 1;
  int n_particles_local = particle_list.size();
  int* buffer_i;
  double* buffer_d;
  // Set up patch info
  int start[2] = {offset, 0};
  int size[2] = {n_particles_local, 3};
  int total_size[2] = {total_particles, 3};

  hid_t t = H5T_NATIVE_DOUBLE;
  if ((fieldname == "type") || (fieldname == "ind") || (fieldname == "fate")) {
    t = H5T_NATIVE_INT;
    buffer_i = new int[n_particles_local];
  }
  else if ((fieldname == "x") || (fieldname == "D") || (fieldname == "x_interact")) {
    n_dims = 2;
    buffer_d = new double[n_particles_local * 3];
  }
  else if ((fieldname == "t") || (fieldname == "e") || (fieldname == "nu") || (fieldname == "gamma") || (fieldname == "dshift") || (fieldname == "dvds")) {
    buffer_d = new double[n_particles_local];
  }
  else {
    std::cerr << "Particle field " << fieldname << " does not exist. Terminating" << std::endl;
    exit(3);
  }

  if (t == H5T_NATIVE_INT)
    readPatch(fname, groupname, fieldname.c_str(), buffer_i, t, n_dims, start, size, total_size);
  else if (t == H5T_NATIVE_DOUBLE)
    readPatch(fname, groupname, fieldname.c_str(), buffer_d, t, n_dims, start, size, total_size);
  else {
    std::cerr << "HDF5 type is wrong" << std::endl;
    exit(3);
  }

  if (fieldname == "type") {
    for (int i = 0; i < n_particles_local; i++) {
      particle_list[i].type = static_cast<PType>(buffer_i[i]);
    }
  }
  else if (fieldname == "x") {
    for (int i = 0; i < n_particles_local; i++) {
      particle_list[i].x[0] = buffer_d[i * 3];
      particle_list[i].x[1] = buffer_d[i * 3 + 1];
      particle_list[i].x[2] = buffer_d[i * 3 + 2];
    }
  }
  else if (fieldname == "D") {
    for (int i = 0; i < n_particles_local; i++) {
      particle_list[i].D[0] = buffer_d[i * 3];
      particle_list[i].D[1] = buffer_d[i * 3 + 1];
      particle_list[i].D[2] = buffer_d[i * 3 + 2];
    }
  }
  else if (fieldname == "x_interact") {
    for (int i = 0; i < n_particles_local; i++) {
      particle_list[i].x_interact[0] = buffer_d[i * 3];
      particle_list[i].x_interact[1] = buffer_d[i * 3 + 1];
      particle_list[i].x_interact[2] = buffer_d[i * 3 + 2];
    }
  }
  else if (fieldname == "ind") {
    for (int i = 0; i < n_particles_local; i++) {
      particle_list[i].ind = buffer_i[i];
    }
  }
  else if (fieldname == "t") {
    for (int i = 0; i < n_particles_local; i++) {
      particle_list[i].t = buffer_d[i];
    }
  }
  else if (fieldname == "e") {
    for (int i = 0; i < n_particles_local; i++) {
      particle_list[i].e = buffer_d[i];
    }
  }
  else if (fieldname == "nu") {
    for (int i = 0; i < n_particles_local; i++) {
      particle_list[i].nu = buffer_d[i];
    }
  }
  else if (fieldname == "gamma") {
    for (int i = 0; i < n_particles_local; i++) {
      particle_list[i].gamma = buffer_d[i];
    }
  }
  else if (fieldname == "dshift") {
    for (int i = 0; i < n_particles_local; i++) {
      particle_list[i].dshift = buffer_d[i];
    }
  }
  else if (fieldname == "dvds") {
    for (int i = 0; i < n_particles_local; i++) {
      particle_list[i].dvds = buffer_d[i];
    }
  }
  else if (fieldname == "fate") {
    for (int i = 0; i < n_particles_local; i++) {
      particle_list[i].fate = static_cast<ParticleFate>(buffer_i[i]);
    }
  }
  else {
    std::cerr << "Particle field " << fieldname << " does not exist. Terminating" << std::endl;
    exit(3);
  }

  if (t == H5T_NATIVE_INT)
    delete[] buffer_i;
  else if (t == H5T_NATIVE_DOUBLE)
    delete[] buffer_d;
  else {
    std::cerr << "HDF5 type is wrong" << std::endl;
    exit(3);
  }
}

void transport::readCheckpointSpectra(std::string fname, bool test) {
  optical_spectrum.readCheckpointSpectrum(fname, "optical spectrum");
  gamma_spectrum.readCheckpointSpectrum(fname, "gamma spectrum");
}
