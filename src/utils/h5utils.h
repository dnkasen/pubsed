#ifndef _H5UTILS_H
#define _H5UTILS_H 1

#include <hdf5.h>

void createFile( char * fname );
void createGroup( char * fname , char * gname );
void createDataset( char * fname , char * gname , char * dname , int dim , hsize_t * fdims , hid_t type );
void writeSimple( char * file , char * group , char * dset , void * data , hid_t type );
void writePatch( char * file , char * group , char * dset , void * data , hid_t type , int dim , int * start , int * loc_size , int * glo_size);

hid_t openH5File( char * fname );
hid_t openH5Group(hid_t h5file, char* group);
hid_t openH5Dset(hid_t h5group, char* dset);
void closeH5File(hid_t h5file);
void closeH5Group(hid_t h5group);
void closeH5Dset(hid_t h5dset);
void getH5dims( hid_t h5grp , char * dset , hsize_t * dims );
void readSimple( char* fname, char* group, char * dset , void * data , hid_t type );
void readPatch( char* fname, char* group, char * dset , void * data , hid_t type , int dim , int * start , int * loc_size , int * glo_size);
#endif
