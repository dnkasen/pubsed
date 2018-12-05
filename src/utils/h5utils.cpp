/* HDF5 helper functions, from Paul Duffell */

#include <hdf5.h>
#include "h5utils.h"

void createFile( char * fname ){
  hid_t h5file = H5Fcreate( fname , H5F_ACC_TRUNC , H5P_DEFAULT , H5P_DEFAULT );
  H5Fclose( h5file );
}

void createGroup( char * fname , char * gname ){
  hid_t h5file = H5Fopen( fname , H5F_ACC_RDWR , H5P_DEFAULT );
  hid_t h5group = H5Gcreate1(h5file, gname, 0);
  H5Gclose( h5group );
  H5Fclose( h5file );
}

void createDataset( char * fname , char * gname , char * dname , int dim , hsize_t * fdims , hid_t type ){
  hid_t h5file  = H5Fopen( fname , H5F_ACC_RDWR , H5P_DEFAULT );
  hid_t h5group = H5Gopen1( h5file , gname );
  hid_t fspace  = H5Screate_simple(dim,fdims,NULL);
  hid_t h5dset  = H5Dcreate1( h5group , dname , type , fspace , H5P_DEFAULT );
  H5Sclose( fspace );
  H5Dclose( h5dset );
  H5Gclose( h5group );
  H5Fclose( h5file );
}

void writeSimple( char * file , char * group , char * dset , void * data , hid_t type ){
  hid_t h5fil = H5Fopen( file , H5F_ACC_RDWR , H5P_DEFAULT );
  hid_t h5grp = H5Gopen1( h5fil , group );
  hid_t h5dst = H5Dopen1( h5grp , dset );

  H5Dwrite( h5dst , type , H5S_ALL , H5S_ALL , H5P_DEFAULT , data );

  H5Dclose( h5dst );
  H5Gclose( h5grp );
  H5Fclose( h5fil );
}

void writePatch( char * file , char * group , char * dset , void * data , hid_t type , int dim , int * start , int * loc_size , int * glo_size){
  hid_t h5fil = H5Fopen( file , H5F_ACC_RDWR , H5P_DEFAULT );
  hid_t h5grp = H5Gopen1( h5fil , group );
  hid_t h5dst = H5Dopen1( h5grp , dset );

  hsize_t mdims[dim];
  hsize_t fdims[dim];

  hsize_t fstart[dim];
  hsize_t fstride[dim];
  hsize_t fcount[dim];
  hsize_t fblock[dim];

  int d;
  for( d=0 ; d<dim ; ++d ){
    mdims[d] = loc_size[d];
    fdims[d] = glo_size[d];

    fstart[d]  = start[d];
    fstride[d] = 1;
    fcount[d]  = loc_size[d];
    fblock[d]  = 1;
  }
  hid_t mspace = H5Screate_simple(dim,mdims,NULL);
  hid_t fspace = H5Screate_simple(dim,fdims,NULL);

  H5Sselect_hyperslab( fspace , H5S_SELECT_SET , fstart , fstride , fcount , fblock );

  H5Dwrite( h5dst , type , mspace , fspace , H5P_DEFAULT , data );

  H5Sclose( mspace );
  H5Sclose( fspace );
  H5Dclose( h5dst );
  H5Gclose( h5grp );
  H5Fclose( h5fil );
}

hid_t openH5File( char * fname ){
  hid_t h5fil = H5Fopen( fname , H5F_ACC_RDWR , H5P_DEFAULT );
  return h5fil;
}

hid_t openH5Group(hid_t h5file, char* group) {
  hid_t h5group = H5Gopen1( h5file , group );
  return h5group;
}

hid_t openH5Dset(hid_t h5group, char* dset) {
  hid_t h5dset = H5Dopen1 (h5group, dset);
  return h5dset;
}

void closeH5File(hid_t h5file) {
  H5Fclose(h5file);
}

void closeH5Group(hid_t h5group) {
  H5Gclose(h5group);
}

void closeH5Dset(hid_t h5dset) {
  H5Dclose(h5dset);
}

void getH5dims( hid_t h5grp , char * dset , hsize_t * dims ){
  hid_t h5dst = H5Dopen1( h5grp , dset );
  hid_t h5spc = H5Dget_space( h5dst );

  H5Sget_simple_extent_dims( h5spc , dims , NULL);

  H5Sclose( h5spc );
  H5Dclose( h5dst );
}

void readSimple( char* fname, char* group, char * dset , void * data , hid_t type ){
  hid_t h5file = H5Fopen( fname , H5F_ACC_RDWR , H5P_DEFAULT );
  hid_t h5group = H5Gopen1( h5file , group );
  hid_t h5dst = H5Dopen1( h5group , dset );

  H5Dread( h5dst , type , H5S_ALL , H5S_ALL , H5P_DEFAULT , data );

  H5Dclose( h5dst );
  H5Gclose( h5group );
  H5Fclose( h5file);
}

void readPatch(char* fname, char* group, char * dset , void * data , hid_t type , int dim , int * start , int * loc_size , int * glo_size){
  hid_t h5file = H5Fopen( fname , H5F_ACC_RDWR , H5P_DEFAULT );
  hid_t h5group = H5Gopen1( h5file , group );
  hid_t h5dst = H5Dopen1( h5group , dset );

  hsize_t mdims[dim];
  hsize_t fdims[dim];

  hsize_t fstart[dim];
  hsize_t fstride[dim];
  hsize_t fcount[dim];
  hsize_t fblock[dim];

  int d;
  for( d=0 ; d<dim ; ++d ){
    mdims[d] = loc_size[d];
    fdims[d] = glo_size[d];

    fstart[d]  = start[d];
    fstride[d] = 1;
    fcount[d]  = loc_size[d];
    fblock[d]  = 1;
  }
  hid_t mspace = H5Screate_simple(dim,mdims,NULL);
  hid_t fspace = H5Screate_simple(dim,fdims,NULL);

  H5Sselect_hyperslab( fspace , H5S_SELECT_SET , fstart , fstride , fcount , fblock );

  H5Dread( h5dst , type , mspace , fspace , H5P_DEFAULT , data );

  H5Sclose( mspace );
  H5Sclose( fspace );
  H5Dclose( h5dst );
  H5Gclose( h5group );
  H5Fclose( h5file);
}

