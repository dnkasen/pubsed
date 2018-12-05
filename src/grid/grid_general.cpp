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
	// read the model file or fill in custom model
	read_model_file(params);

	// complain if the grid is obviously not right
	if(z.size()==0)
	{
		std::cerr << "Error: there are no grid zones." << std::endl;
		exit(5);
		n_zones = z.size();
	}

  writeCheckpointZones("zones.h5");
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
  int my_rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
  /* Mercifully, only rank 0 has to do any of this */
  if (my_rank == 0) {
    createFile(fname);
    createGroup(fname, "zones");
    writeScalarZoneProp(fname, "v");
    writeScalarZoneProp(fname, "rho");
    writeScalarZoneProp(fname, "cs");
    writeScalarZoneProp(fname, "p_gas");
    writeScalarZoneProp(fname, "T_gas");
    writeScalarZoneProp(fname, "mu");
    writeScalarZoneProp(fname, "erad");
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
  if (fieldname == "rho")
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
