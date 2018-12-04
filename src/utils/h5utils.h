#ifndef _H5UTILS_H
#define _H5UTILS_H 1

#include <hdf5.h>

// Creates and closes a new HDF5 file called fname
void createFile( char * fname );

// Creates a group called gname in the existing HDF5 file fname
void createGroup( char * fname , char * gname );

// Creates a dataset called dname in the existing HDF5 group and file named gname
// and fname. The created dataset will be have dimension fdims[0] by fdims[1] by ... fdims[dim-1].
// dim is the number of dimensions. type is the HDF5 data type name, e.g. H5T_NATIVE_DOUBLE.
void createDataset( char * fname , char * gname , char * dname , int dim , hsize_t * fdims , hid_t type );

// Writes to the extant dataset file/group/dset. data is a pointer to a buffer
// containing the data to be written. type is the HDF5 data type name, e.g. H5T_NATIVE_DOUBLE
void writeSimple( char * file , char * group , char * dset , void * data , hid_t type );

// Writes to a chunk of the extant dataset file/group/dset. data is a pointer
// to a buffer containing the data to be written. type is the HDF5 data type name,
// e.g. H5T_NATIVE_DOUBLE. dim is the number of dimensions in the dataset.
// glo_size is an array of the size of the whole dataset in each dimension.
// loc_size is an array of the size of the patch of data to write in each dmension.
// start is the offset of the patch in each dimension.
void writePatch( char * file , char * group , char * dset , void * data , hid_t type , int dim , int * start , int * loc_size , int * glo_size);

// Opens and returns a handle to hdf5 file fname
hid_t openH5File( char * fname );

// Opens and returns a handle to the group called group in the file h5file
hid_t openH5Group(hid_t h5file, char* group);

// Opens and returns a handle to the dataset called dset in the group h5group
hid_t openH5Dset(hid_t h5group, char* dset);

// Closes file with handle h5file
void closeH5File(hid_t h5file);

// Closes group with handle h5group
void closeH5Group(hid_t h5group);

// Closes dataset with handle h5dset
void closeH5Dset(hid_t h5dset);

// Gets dimensions of the dataset dset in h5grp. Dimensions are returned in
// the buffer dims
void getH5dims( hid_t h5grp , char * dset , hsize_t * dims );

// Reads in the entire content of the dataset at fname/group/dset. Data is read
// into the buffer data
void readSimple( char* fname, char* group, char * dset , void * data , hid_t type );

// Reads in a patch of data from the dataset fname/group/dset to the buffer data.
// dim is the number of dimensions in the dataset. The patch to read has size
// loc_size, and starts at start. glo_size is the size of the whole dataset
void readPatch( char* fname, char* group, char * dset , void * data , hid_t type , int dim , int * start , int * loc_size , int * glo_size);
#endif
