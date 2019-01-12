#include <cstdlib>
#include <math.h>

#include "grid_general.h"
#include "physical_constants.h"
namespace pc = physical_constants;

//------------------------------------------------------------
// initialize the grid
//------------------------------------------------------------
void grid_general::init(ParameterReader* params)
{
  MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
  MPI_Comm_size( MPI_COMM_WORLD, &nproc );

  // If it's a restart, restart the grid. Otherwise read in the model file
  if (params->getScalar<int>("do_restart"))
    restartGrid(params);
  else
	  read_model_file(params);

	// complain if the grid is obviously not right
	if(z.size()==0)
	{
		std::cerr << "Error: there are no grid zones." << std::endl;
		exit(5);
		n_zones = z.size();
	}

}

void grid_general::write_hdf5_plotfile_zones
(hid_t file_id, hsize_t *dims_g, int ndims, double tt)
{
	// print out time
	hsize_t  dims_t[1]={1};
  	float time_a[1];
  	time_a[0] = tt;
  	H5LTmake_dataset(file_id,"time",1,dims_t,H5T_NATIVE_FLOAT,time_a);

	// print out zone arrays
  	float *arr = new float[n_zones];

	// print out rho
	for (int i=0;i<n_zones;++i) arr[i] = z[i].rho;
	H5LTmake_dataset(file_id,"rho",ndims,dims_g,H5T_NATIVE_FLOAT,arr);

 	// print out vel
	for (int i=0;i<n_zones;++i) arr[i] = z[i].v[0];
	H5LTmake_dataset(file_id,"velr",ndims,dims_g,H5T_NATIVE_FLOAT,arr);

	if (ndims > 1)
	{
		for (int i=0;i<n_zones;++i) arr[i] = z[i].v[2];
		H5LTmake_dataset(file_id,"velz",ndims,dims_g,H5T_NATIVE_FLOAT,arr);
	}

	// print out T_rad
	for (int i=0;i<n_zones;++i) arr[i] = pow(z[i].e_rad/pc::a,0.25);
	H5LTmake_dataset(file_id,"T_rad",ndims,dims_g,H5T_NATIVE_FLOAT,arr);

	// print out T_gas
	for (int i=0;i<n_zones;++i) arr[i] = z[i].T_gas;
	H5LTmake_dataset(file_id,"T_gas",ndims,dims_g,H5T_NATIVE_FLOAT,arr);

	// print out radioactive deposition
	for (int i=0;i<n_zones;++i) arr[i] = z[i].L_radio_dep;
	H5LTmake_dataset(file_id,"e_nuc_dep",ndims,dims_g,H5T_NATIVE_FLOAT,arr);

	// print out radioactive emission
	for (int i=0;i<n_zones;++i) arr[i] = z[i].L_radio_emit;
	H5LTmake_dataset(file_id,"e_nuc_emit",ndims,dims_g,H5T_NATIVE_FLOAT,arr);

  	delete [] arr;
}

void grid_general::writeCheckpointZones(std::string fname) {
  /* Mercifully, only rank 0 has to do any of this */
  if (my_rank == 0) {
    std::cerr << "creating group" << std::endl;
    createGroup(fname, "zones");
    std::cerr << "creating scalar prop" << std::endl;
    writeScalarZoneProp(fname, "v");
    writeScalarZoneProp(fname, "rho");
    writeScalarZoneProp(fname, "cs");
    writeScalarZoneProp(fname, "p_gas");
    writeScalarZoneProp(fname, "T_gas");
    writeScalarZoneProp(fname, "mu");
    writeScalarZoneProp(fname, "e_rad");
    writeScalarZoneProp(fname, "e_abs");
    writeScalarZoneProp(fname, "fx_rad");
    writeScalarZoneProp(fname, "fy_rad");
    writeScalarZoneProp(fname, "fz_rad");
    writeScalarZoneProp(fname, "eps_imc");
    writeScalarZoneProp(fname, "L_thermal");
    writeScalarZoneProp(fname, "L_radio_emit");
    writeScalarZoneProp(fname, "L_radio_dep");
    writeScalarZoneProp(fname, "G1");
    writeScalarZoneProp(fname, "G2");
    writeScalarZoneProp(fname, "G3");
    writeScalarZoneProp(fname, "P11");
    writeScalarZoneProp(fname, "P12");
    writeScalarZoneProp(fname, "P13");
    writeScalarZoneProp(fname, "P22");
    writeScalarZoneProp(fname, "P23");
    writeScalarZoneProp(fname, "P33");
    writeVectorZoneProp(fname, "X_gas");
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

void grid_general::writeScalarZoneProp(std::string fname, std::string fieldname) {
  int ndim1 = 1;
  hsize_t dims1[1] =  {n_zones};
  int ndim3 = 2;
  hsize_t dims3[2] = {n_zones, 3};
  real* buffer = new real[n_zones * 3];
  hid_t t;

  if (std::is_same<real, float>())
    t = H5T_NATIVE_FLOAT;
  else if (std::is_same<real, double>())
    t = H5T_NATIVE_DOUBLE;
  else {
    std::cerr << "real type not known. Cannot set up HDF5 data sets" << std::endl;
  }

  if (fieldname == "v") {
    createDataset(fname, "zones", "v", ndim3, dims3, t);
  }
  else {
    createDataset(fname, "zones", fieldname, ndim1, dims1, t);
  }

  if (fieldname == "v") {
    for (int i = 0; i < n_zones; i++) {
      buffer[i * 3] = z[i].v[0];
      buffer[i * 3 + 1] = z[i].v[1];
      buffer[i * 3 + 2] = z[i].v[2];
    }
  }
  else if (fieldname == "rho")
    for (int i = 0; i < n_zones; i++) buffer[i] = z[i].rho;
  else if (fieldname == "cs")
    for (int i = 0; i < n_zones; i++) buffer[i] = z[i].cs;
  else if (fieldname == "p_gas")
    for (int i = 0; i < n_zones; i++) buffer[i] = z[i].p_gas;
  else if (fieldname == "T_gas")
    for (int i = 0; i < n_zones; i++) buffer[i] = z[i].T_gas;
  else if (fieldname == "mu")
    for (int i = 0; i < n_zones; i++) buffer[i] = z[i].mu;
  else if (fieldname == "e_rad")
    for (int i = 0; i < n_zones; i++) buffer[i] = z[i].e_rad;
  else if (fieldname == "e_abs")
    for (int i = 0; i < n_zones; i++) buffer[i] = z[i].e_abs;
  else if (fieldname == "fx_rad")
    for (int i = 0; i < n_zones; i++) buffer[i] = z[i].fx_rad;
  else if (fieldname == "fy_rad")
    for (int i = 0; i < n_zones; i++) buffer[i] = z[i].fy_rad;
  else if (fieldname == "fz_rad")
    for (int i = 0; i < n_zones; i++) buffer[i] = z[i].fz_rad;
  else if (fieldname == "eps_imc")
    for (int i = 0; i < n_zones; i++) buffer[i] = z[i].eps_imc;
  else if (fieldname == "L_thermal")
    for (int i = 0; i < n_zones; i++) buffer[i] = z[i].L_thermal;
  else if (fieldname == "L_radio_emit")
    for (int i = 0; i < n_zones; i++) buffer[i] = z[i].L_radio_emit;
  else if (fieldname == "L_radio_dep")
    for (int i = 0; i < n_zones; i++) buffer[i] = z[i].L_radio_dep;
  else if (fieldname == "G1")
    for (int i = 0; i < n_zones; i++) buffer[i] = z[i].G1;
  else if (fieldname == "G2")
    for (int i = 0; i < n_zones; i++) buffer[i] = z[i].G2;
  else if (fieldname == "G3")
    for (int i = 0; i < n_zones; i++) buffer[i] = z[i].G3;
  else if (fieldname == "P11")
    for (int i = 0; i < n_zones; i++) buffer[i] = z[i].P11;
  else if (fieldname == "P12")
    for (int i = 0; i < n_zones; i++) buffer[i] = z[i].P12;
  else if (fieldname == "P13")
    for (int i = 0; i < n_zones; i++) buffer[i] = z[i].P13;
  else if (fieldname == "P22")
    for (int i = 0; i < n_zones; i++) buffer[i] = z[i].P22;
  else if (fieldname == "P23")
    for (int i = 0; i < n_zones; i++) buffer[i] = z[i].P23;
  else if (fieldname == "P33")
    for (int i = 0; i < n_zones; i++) buffer[i] = z[i].P33;
  else {
    std::cerr << "Field name " << fieldname << " not known." <<std::endl;
    exit(4);
  }
  writeSimple(fname, "zones", fieldname, buffer, t);
  delete [] buffer;
}

void grid_general::writeVectorZoneProp(std::string fname, std::string fieldname) {
  hid_t t;
  if (std::is_same<real, float>())
    t = H5T_NATIVE_FLOAT;
  else if (std::is_same<real, double>())
    t = H5T_NATIVE_DOUBLE;
  else {
    std::cerr << "real type not known. Cannot set up HDF5 data sets" << std::endl;
  }

  real* buffer = new real[n_zones * n_elems];
  if (fieldname == "X_gas") {
    int ndim = 2;
    hsize_t dims[2] = {n_zones, n_elems};
    createDataset(fname, "zones", fieldname, ndim, dims, t);
    for (int i = 0; i < n_zones; i++) {
      for (int j = 0; j < n_elems; j++) {
        buffer[i * n_elems + j] = z[i].X_gas[j];
      }
    }
  }
  else {
    std::cerr << "vector zone property " << fieldname << " unknown. Terminating." <<std::endl;
  }

  writeSimple(fname, "zones", fieldname, buffer, t);
}

void grid_general::readCheckpointZones(std::string fname, bool test) {
  /* To avoid doing extra communication, each rank will read in its own grid data. */
  for (int rank = 0; rank < nproc; rank++) {
    if (my_rank == rank) {
      /* Adjust vector lengths */
      hsize_t dims[2];
      getH5dims(fname, "zones", "X_gas", dims);
      n_zones = dims[0];
      n_elems = dims[1];
      z_new.resize(n_zones);

      readScalarZoneProp(fname, "v");
      readScalarZoneProp(fname, "rho");
      readScalarZoneProp(fname, "cs");
      readScalarZoneProp(fname, "p_gas");
      readScalarZoneProp(fname, "T_gas");
      readScalarZoneProp(fname, "mu");
      readScalarZoneProp(fname, "e_rad");
      readScalarZoneProp(fname, "e_abs");
      readScalarZoneProp(fname, "fx_rad");
      readScalarZoneProp(fname, "fy_rad");
      readScalarZoneProp(fname, "fz_rad");
      readScalarZoneProp(fname, "eps_imc");
      readScalarZoneProp(fname, "L_thermal");
      readScalarZoneProp(fname, "L_radio_emit");
      readScalarZoneProp(fname, "L_radio_dep");
      readScalarZoneProp(fname, "G1");
      readScalarZoneProp(fname, "G2");
      readScalarZoneProp(fname, "G3");
      readScalarZoneProp(fname, "P11");
      readScalarZoneProp(fname, "P12");
      readScalarZoneProp(fname, "P13");
      readScalarZoneProp(fname, "P22");
      readScalarZoneProp(fname, "P23");
      readScalarZoneProp(fname, "P33");
      readVectorZoneProp(fname, "X_gas");

      if (not test) {
        z = z_new;
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
}

void grid_general::readScalarZoneProp(std::string fname, std::string fieldname) {
  real* buffer = new real[n_zones * 3];
  hid_t t;
  if (std::is_same<real, float>())
    t = H5T_NATIVE_FLOAT;
  else if (std::is_same<real, double>())
    t = H5T_NATIVE_DOUBLE;
  else {
    std::cerr << "real type not known. Cannot set up HDF5 data sets" << std::endl;
  }
  /* Fill buffer */
  readSimple(fname, "zones", fieldname, buffer, t);
  if (fieldname == "v") {
    for (int i = 0; i < n_zones; i++) {
      z_new[i].v[0] = buffer[i * 3];
      z_new[i].v[1] = buffer[i * 3 + 1];
      z_new[i].v[2] = buffer[i * 3 + 2];
    }
  }
  else if (fieldname == "rho")
    for (int i = 0; i < n_zones; i++) z_new[i].rho = buffer[i];
  else if (fieldname == "cs")
    for (int i = 0; i < n_zones; i++) z_new[i].cs = buffer[i];
  else if (fieldname == "p_gas")
    for (int i = 0; i < n_zones; i++) z_new[i].p_gas = buffer[i];
  else if (fieldname == "T_gas")
    for (int i = 0; i < n_zones; i++) z_new[i].T_gas= buffer[i];
  else if (fieldname == "mu")
    for (int i = 0; i < n_zones; i++) z_new[i].mu = buffer[i];
  else if (fieldname == "e_rad")
    for (int i = 0; i < n_zones; i++) z_new[i].e_rad = buffer[i];
  else if (fieldname == "e_abs")
    for (int i = 0; i < n_zones; i++) z_new[i].e_abs = buffer[i];
  else if (fieldname == "fx_rad")
    for (int i = 0; i < n_zones; i++) z_new[i].fx_rad = buffer[i];
  else if (fieldname == "fy_rad")
    for (int i = 0; i < n_zones; i++) z_new[i].fy_rad = buffer[i];
  else if (fieldname == "fz_rad")
    for (int i = 0; i < n_zones; i++) z_new[i].fz_rad = buffer[i];
  else if (fieldname == "eps_imc")
    for (int i = 0; i < n_zones; i++) z_new[i].eps_imc = buffer[i];
  else if (fieldname == "L_thermal")
    for (int i = 0; i < n_zones; i++) z_new[i].L_thermal = buffer[i];
  else if (fieldname == "L_radio_emit")
    for (int i = 0; i < n_zones; i++) z_new[i].L_radio_emit = buffer[i];
  else if (fieldname == "L_radio_dep")
    for (int i = 0; i < n_zones; i++) z_new[i].L_radio_dep = buffer[i];
  else if (fieldname == "G1")
    for (int i = 0; i < n_zones; i++) z_new[i].G1 = buffer[i];
  else if (fieldname == "G2")
    for (int i = 0; i < n_zones; i++) z_new[i].G2 = buffer[i];
  else if (fieldname == "G3")
    for (int i = 0; i < n_zones; i++) z_new[i].G3 = buffer[i];
  else if (fieldname == "P11")
    for (int i = 0; i < n_zones; i++) z_new[i].P11 = buffer[i];
  else if (fieldname == "P12")
    for (int i = 0; i < n_zones; i++) z_new[i].P12 = buffer[i];
  else if (fieldname == "P13")
    for (int i = 0; i < n_zones; i++) z_new[i].P13 = buffer[i];
  else if (fieldname == "P22")
    for (int i = 0; i < n_zones; i++) z_new[i].P22 = buffer[i];
  else if (fieldname == "P23")
    for (int i = 0; i < n_zones; i++) z_new[i].P23 = buffer[i];
  else if (fieldname == "P33")
    for (int i = 0; i < n_zones; i++) z_new[i].P33 = buffer[i];
  else {
    std::cerr << "Unknown zone field name " << fieldname << std::endl;
    exit(3);
  }

  delete[] buffer;
}

void grid_general::readVectorZoneProp(std::string fname, std::string fieldname) {
  hid_t t;
  if (std::is_same<real, float>())
    t = H5T_NATIVE_FLOAT;
  else if (std::is_same<real, double>())
    t = H5T_NATIVE_DOUBLE;
  else {
    std::cerr << "real type not known. Cannot read HDF5 data sets" << std::endl;
  }

  real* buffer = new real[n_zones * n_elems];
  readSimple(fname, "zones", fieldname, buffer, t);

  if (fieldname == "X_gas") {
    for (int i = 0; i < n_zones; i++) {
      z_new[i].X_gas.resize(n_elems);
      for (int j = 0; j < n_elems; j++) {
        z_new[i].X_gas[j] = buffer[i * n_elems + j];
      }
    }
  }
  else {
    std::cerr << "vector zone property " << fieldname << " unknown. Terminating." <<std::endl;
  }

  delete[] buffer;
}

void grid_general::writeCheckpointGeneralGrid(std::string fname) {
  /* Mercifully, only rank 0 has to do any of this, but calling function will handle telling
   * which ranks to do what */
  createGroup(fname, "grid");
  hsize_t single_val = 1;
  hsize_t elem_dim = n_elems;

  createDataset(fname, "grid", "t_now", 1, &single_val, H5T_NATIVE_DOUBLE);
  writeSimple(fname, "grid", "t_now", &t_now, H5T_NATIVE_DOUBLE);
  createDataset(fname, "grid", "n_zones", 1, &single_val, H5T_NATIVE_INT);
  writeSimple(fname, "grid", "n_zones", &n_zones, H5T_NATIVE_INT);

  createDataset(fname, "grid", "n_elems", 1, &single_val, H5T_NATIVE_INT);
  writeSimple(fname, "grid", "n_elems", &n_elems, H5T_NATIVE_INT);

  writeVector(fname, "grid", "elems_Z", elems_Z, H5T_NATIVE_INT);
  writeVector(fname, "grid", "elems_A", elems_A, H5T_NATIVE_INT);

}

void grid_general::readCheckpointGeneralGrid(std::string fname, bool test) {
  std::cerr << "t_now" << std::endl;
  readSimple(fname, "grid", "t_now", &t_now_new, H5T_NATIVE_DOUBLE);
  std::cerr << "nzones" << std::endl;
  readSimple(fname, "grid", "n_zones", &n_zones_new, H5T_NATIVE_INT);
  std::cerr << "nelems" << std::endl;
  readSimple(fname, "grid", "n_elems", &n_elems_new, H5T_NATIVE_INT);
  std::cerr << "elems z" << std::endl;
  readVector(fname, "grid", "elems_Z", elems_Z_new, H5T_NATIVE_INT);
  std::cerr << "elems A" << std::endl;
  readVector(fname, "grid", "elems_A", elems_A_new, H5T_NATIVE_INT);
  if (not test) {
    t_now = t_now_new;
    n_zones = n_zones_new;
    n_elems = n_elems_new;
    elems_A = elems_A_new;
    elems_Z = elems_Z_new;
  }
}

void grid_general::write_integrated_quantities(int iw, double tt)
{
	// open output file
  	FILE* fout = NULL;
  	if (iw == 0) fout = fopen(TIME_PLOTFILE_NAME,"w");
  	else  fout = fopen(TIME_PLOTFILE_NAME,"a");
  	if (fout == NULL) return;

  	// write header
  	if (iw == 0)
      fprintf(fout,"#    time(sec)     E_radiation      L_nuc_dep      L_nuc_emit      Mass\n");

	// integrated qunaitites
  	double L_dep = 0, L_emit = 0, E_rad = 0, mass = 0;
  	for (int i=0;i<n_zones;++i)
 	{
    	double vol = zone_volume(i);
    	L_dep  += z[i].L_radio_dep*vol;
    	L_emit += z[i].L_radio_emit*vol;
    	E_rad  += z[i].e_rad*vol;
    	mass   += z[i].rho*vol;
 	 }
    fprintf(fout,"%15.6e %15.6e %15.6e %15.6e %15.6e\n",tt,E_rad,L_dep,L_emit,mass);
 	fclose(fout);
}
