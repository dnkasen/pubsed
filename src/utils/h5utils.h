#ifndef _H5UTILS_H
#define _H5UTILS_H 1

#include <vector>
#include <iostream>
#include <hdf5.h>

// Returns variable length string HDF5 type
hid_t H5VarString();
// Creates and closes a new HDF5 file called fname
void createFile(std::string fname);

// Creates a group called gname in the existing HDF5 file fname
void createGroup(std::string fname, std::string gname);
void createGroup(hid_t parent, std::string gname);

// Creates a dataset called dname in the existing HDF5 group and file named gname
// and fname. The created dataset will be have dimension fdims[0] by fdims[1] by ... fdims[dim-1].
// dim is the number of dimensions. type is the HDF5 data type name, e.g. H5T_NATIVE_DOUBLE.
void createDataset(std::string fname, std::string gname, std::string dname, int dim, hsize_t* fdims, hid_t type);
void createDataset(hid_t h5group, std::string dname, int dim, hsize_t* fdims, hid_t type);

// Writes to the extant dataset file/group/dset. data is a pointer to a buffer
// containing the data to be written. type is the HDF5 data type name, e.g. H5T_NATIVE_DOUBLE
void writeSimple(std::string file, std::string group, std::string dset, void* data, hid_t type);
void writeSimple(hid_t h5grp, std::string dset, void* data, hid_t type);

// Writes to a chunk of the extant dataset file/group/dset. data is a pointer
// to a buffer containing the data to be written. type is the HDF5 data type name,
// e.g. H5T_NATIVE_DOUBLE. dim is the number of dimensions in the dataset.
// glo_size is an array of the size of the whole dataset in each dimension.
// loc_size is an array of the size of the patch of data to write in each dmension.
// start is the offset of the patch in each dimension.
void writePatch(std::string file, std::string group, std::string dset, void* data, hid_t type, int dim, int* start, int* loc_size, int* glo_size);

template <typename T>
void writeVector(std::string file, std::string group, std::string dset, std::vector<T>& vec, hid_t t) {
  T* buffer = vec.data();
  int ndims = 1;
  hsize_t dims = vec.size();
  createDataset(file, group, dset, ndims, &dims, t);
  writeSimple(file, group, dset, buffer, t);
}

template <typename T>
void writeVector(hid_t h5group, std::string dset, std::vector<T>& vec, hid_t t) {
  T* buffer = vec.data();
  int ndims = 1;
  hsize_t dims = vec.size();
  createDataset(h5group, dset, ndims, &dims, t);
  writeSimple(h5group, dset, buffer, t);
}
// Opens and returns a handle to hdf5 file fname
hid_t openH5File(std::string fname);

// Opens and returns a handle to the group called group in the file h5file
hid_t openH5Group(hid_t h5file, std::string group);

// Opens and returns a handle to the dataset called dset in the group h5group
hid_t openH5Dset(hid_t h5group, std::string dset);

// Closes file with handle h5file
void closeH5File(hid_t h5file);

// Closes group with handle h5group
void closeH5Group(hid_t h5group);

// Closes dataset with handle h5dset
void closeH5Dset(hid_t h5dset);

// Gets dimensions of the dataset dset in h5grp. Dimensions are returned in
// the buffer dims
void getH5dims(hid_t h5grp, std::string dset, hsize_t* dims);
void getH5dims(std::string fname, std::string gname, std::string dset, hsize_t* dims);

// Reads in the entire content of the dataset at fname/group/dset or h5grp/dset. Data are read
// into the buffer data
void readSimple(hid_t h5grp, std::string dset, void* data, hid_t type);
void readSimple(std::string fname, std::string group, std::string dset, void* data, hid_t type);

template <typename T>
void readVector(std::string file, std::string group, std::string dset, std::vector<T>& vec, hid_t t) {
  hsize_t dim;
  getH5dims(file, group, dset, &dim);
  std::cerr << dim << std::endl;

  T* buffer = new T[dim];

  readSimple(file, group, dset, buffer, t);
  vec.resize(dim);
  for (int i = 0; i < dim; i++) {
    vec[i] = buffer[i];
  }
  delete[] buffer;
}

template <typename T>
void readVector(hid_t h5parent, std::string dset, std::vector<T>& vec, hid_t t) {
  hsize_t dim;
  getH5dims(h5parent, dset, &dim);

  T* buffer = new T[dim];

  readSimple(h5parent, dset, buffer, t);
  vec.resize(dim);
  for (int i = 0; i < dim; i++) {
    vec[i] = buffer[i];
  }
  delete[] buffer;
}

// Reads in a patch of data from the dataset fname/group/dset to the buffer data.
// dim is the number of dimensions in the dataset. The patch to read has size
// loc_size, and starts at start. glo_size is the size of the whole dataset
void readPatch(hid_t h5grp, std::string dset, void* data, hid_t type, int dim, int* start, int* loc_size, int* glo_size);
void readPatch(std::string fname, std::string group, std::string dset, void* data, hid_t type, int dim, int* start, int* loc_size, int* glo_size);
#endif
