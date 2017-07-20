#include <cstdlib>
#include <mpi.h>
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
		std::cout << "Error: there are no grid zones." << std::endl;
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

	// Close the file 
	H5Fclose(file_id);
  	
  	delete [] arr;
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
   	 fprintf(fout,"#    time(sec)     E_radiation      L_nuc_dep      L_nuc_emit\n");

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
  	fprintf(fout,"%15.6e %15.6e %15.6e %15.6e\n",tt,E_rad,L_dep,L_emit);
 	fclose(fout);
}
