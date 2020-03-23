
#include <algorithm>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include "locate_array.h"

using namespace std;




//---------------------------------------------------------
// Just allocation the memory for this
//---------------------------------------------------------
void locate_array::init(const int n)
{
  x_.assign(n,0);
}

//---------------------------------------------------------
// Initialize with start, stop and delta
// if start==stop make it a catch-all
//---------------------------------------------------------
void locate_array::init(const double start, const double stop, const double del)
{
  if(start>=stop){
    locate_type_ = none;
    x_.resize(1);
    min_ = -numeric_limits<double>::infinity();
    x_[0] = numeric_limits<double>::infinity();
  }

  else{
    locate_type_ = do_lin;
    del_ = del;
    int n = ceil( (stop-start)/del );
    n = max(n,1);
    x_.resize(n);
    do_log_interpolate_ = 0;

    min_ = start;
    for (int i=0; i<n-1; i++) x_[i] = start + (i+1)*del;
    x_[n-1] = stop;
  }
}


//---------------------------------------------------------
// Initialize with start, stop and delta
// logarithmic spacing
// if start==stop make it a catch-all
//---------------------------------------------------------
void locate_array::log_init(const double start, const double stop, const double del)
{
  if(start>=stop){
    locate_type_ = none;
    x_.resize(1);
    min_ = -numeric_limits<double>::infinity();
    x_[0] = numeric_limits<double>::infinity();
  }

  else{
    locate_type_ = do_log;
    del_ = del;
    do_log_interpolate_ = 0;

    std::vector<double> arr;
    min_ =start;
    double v = start + start*del;
    while (v  < stop)
    {
      arr.push_back(v);
      v = v + v*del;
    }
    arr.push_back(stop);
    int n = arr.size();
    x_.resize(n);
    for (int i=0;i<n;i++) x_[i] = arr[i];
  }
}

//---------------------------------------------------------
// Initialize with start, stop and n_pts
// if start==stop make it a catch-all
// if n==0 make it a catch-all
//---------------------------------------------------------
void locate_array::init(const double start, const double stop, const int n)
{
  if(start==stop || n==0){
    locate_type_ = none;
    x_.resize(1);
    min_ = -numeric_limits<double>::infinity();
    x_[0] = numeric_limits<double>::infinity();
  }

  else{
    locate_type_ = do_lin;
    del_ = (stop - start)/(double)n;
    x_.resize(n);
    do_log_interpolate_ = 0;

    min_ = start;
    #pragma omp parallel for
    for (int i=0; i<n-1; i++) x_[i] = start + (i+1)*del_;
    x_[n-1] = stop;
  }
}

//---------------------------------------------------------
// Initialize with passed vector
//---------------------------------------------------------
void locate_array::init(const std::vector<double> &a, const double minval)
{
  locate_type_ = flex;
  min_ = minval;
  do_log_interpolate_ = 0;
  x_.assign(a.begin(), a.end());
}


//---------------------------------------------------------
// Initialize with pointer to an array
//---------------------------------------------------------
void locate_array::init(const double* a, const int n, const double minval)
{
  locate_type_ = flex;
  min_ = minval;
  do_log_interpolate_ = 0;
  x_.resize(n);
  for (int i = 0; i < n; i++) x_[i] = a[i];
}

//---------------------------------------------------------
// copy from another
//---------------------------------------------------------
void locate_array::copy(locate_array l)
{
  min_ = l.min_;
  del_ = l.del_;
  do_log_interpolate_ = l.do_log_interpolate_;
  locate_type_ = l.locate_type_;
  x_.resize(l.size());
  for (int i=0;i<l.size();i++) x_[i] = l.x_[i];
}

bool locate_array::is_equal(locate_array l, bool complain) {
  bool equal = true;
  if (min_ != l.min_) {
    if (complain) std::cerr << "locate array minima are different" << std::endl;
    equal = false;
  }
  if (do_log_interpolate_ != l.do_log_interpolate_) {
    if (complain) std::cerr << "locate array do_log_interpolate are different" << do_log_interpolate_ << " " << l.do_log_interpolate_ << std::endl;
    equal = false;
  }
  if (del_ != l.del_) {
    if (complain) std::cerr << "locate array dels are different" << std::endl;
    equal = false;
  }
  if (locate_type_ != l.locate_type_) {
    if (complain) std::cerr << "locate array types are different" << std::endl;
    equal = false;
  }
  if (x_.size() != l.x_.size()) {
    if (complain) std::cerr << "locate array array sizes are different" << std::endl;
    equal = false;
  }
  for (int i = 0; i < x_.size(); i++) {
    if (x_[i] != l.x_[i]) {
      if (complain) std::cerr << "locate array array elements are different" << std::endl;
      equal = false;
    }
  }
  if (!equal && complain) {
    print();
  }
  return equal;
}

    

//---------------------------------------------------------
// locate (return closest index above the value)
// if off left side of boundary, returns 0
// if off right side of boundary, returns size
// Warning: The returned index will be out of bounds
// if off right hand side.
// Bins are closed to the left and open to the right
//---------------------------------------------------------
int locate_array::locate(const double xval) const
{
  // First handle some trivial cases
  if (size() == 1) return 0;
  if (xval >= maxval()) return size();
  if (xval < minval()) return 0;
  int ind;
  if (locate_type_ == flex) {
    ind = upper_bound(x_.begin(), x_.end(), xval) - x_.begin();
  }
  else if (locate_type_ == do_lin) {
    ind = floor((xval - min_) / del_);
  }
  else if (locate_type_ == do_log) {
    ind = floor(log(xval/min_) / log(1 + del_));
  }
  else if (locate_type_ == none) {
    ind = 0;
  }
  else {
    std::cerr << "locate_type not recognized, falling back to flex" << std::endl;
    ind = upper_bound(x_.begin(), x_.end(), xval) - x_.begin();
  }
  if (locate_type_ != none && ind != 0 && ind != size() && (left(ind) > xval or right(ind) <= xval)) {
    std::cerr << "Calculated index incorrect. " << "type: " << locate_type_ << " ind: " << ind
      << " left: " << left(ind) << " val: " << xval << " right: " << right(ind) << " min del: "
      <<  min_ << " " <<  del_ << std::endl;
    ind = upper_bound(x_.begin(), x_.end(), xval) - x_.begin();
  }
  return ind;
}

//---------------------------------------------------------
// locate (return closest index below the value)
// if off left side of boundary, returns 0
// if off right side of boundary, returns size-1
//---------------------------------------------------------
int locate_array::locate_within_bounds(const double xval) const
{
  int ind = locate(xval);
  if (ind == x_.size()) return x_.size()-1;
  return ind;
}



//---------------------------------------------------------
// sample uniformally in zone
//---------------------------------------------------------
double locate_array::sample(const int i, const double rand) const
{
  if (x_.size() == 1) return 0;
  if (i == 0) return min_    + (x_[0] - min_   )*rand;
  else return        x_[i-1] + (x_[i] - x_[i-1])*rand;
}

//---------------------------------------------------------
// simple printout
//---------------------------------------------------------
void locate_array::print() const
{
  printf("# Print Locate Array; n_elements = %lu\n",x_.size());
  printf("min %12.4e\n",min_);
  printf("del %12.4e\n",del_);
  printf("type %4d\n",locate_type_);
  for (int i=0;i<x_.size();i++)
    printf("%4d %12.4e\n",i,x_[i]);
}
  
void locate_array::writeCheckpoint(std::string fname, std::string gname, std::string dset) {
  hsize_t single_val = 1;
  hid_t h5file = openH5File(fname);
  hid_t h5group = openH5Group(h5file, gname);
  createGroup(h5group, dset);
  hid_t h5_locatearray_group = openH5Group(h5group, dset);
  
  writeVector(h5_locatearray_group, "x", x_, H5T_NATIVE_DOUBLE);

  createDataset(h5_locatearray_group, "min", 1, &single_val, H5T_NATIVE_DOUBLE); 
  writeSimple(h5_locatearray_group, "min", &min_, H5T_NATIVE_DOUBLE);
  createDataset(h5_locatearray_group, "do_log_interpolate", 1, &single_val, H5T_NATIVE_INT);
  writeSimple(h5_locatearray_group, "do_log_interpolate", &do_log_interpolate_, H5T_NATIVE_INT);
  createDataset(h5_locatearray_group, "del", 1, &single_val, H5T_NATIVE_DOUBLE); 
  writeSimple(h5_locatearray_group, "del", &del_, H5T_NATIVE_DOUBLE);
  createDataset(h5_locatearray_group, "locate_type", 1, &single_val, H5T_NATIVE_INT); 
  writeSimple(h5_locatearray_group, "locate_type", &locate_type_, H5T_NATIVE_INT);
  
  closeH5Group(h5_locatearray_group);
  closeH5Group(h5group);
  closeH5File(h5file);

}

void locate_array::readCheckpoint(std::string fname, std::string gname, std::string dset) {
  hid_t h5file = openH5File(fname);
  hid_t h5group = openH5Group(h5file, gname);
  hid_t h5_locatearray_group = openH5Group(h5group, dset);

  readVector(h5_locatearray_group, "x", x_, H5T_NATIVE_DOUBLE);
  readSimple(h5_locatearray_group, "min", &min_, H5T_NATIVE_DOUBLE);
  readSimple(h5_locatearray_group, "do_log_interpolate", &do_log_interpolate_, H5T_NATIVE_INT);
  readSimple(h5_locatearray_group, "del", &del_, H5T_NATIVE_DOUBLE);
  readSimple(h5_locatearray_group, "locate_type", &locate_type_, H5T_NATIVE_INT);

  closeH5Group(h5_locatearray_group);
  closeH5Group(h5group);
  closeH5File(h5file);
}

void locate_array::swap(locate_array new_array){
  // swap the vectors
  x_.swap(new_array.x_);

  // swap the minimum values
  double min_tmp = min_;
  min_ = new_array.min_;
  new_array.min_ = min_tmp;

  // swap the do_log_interpolate parameters
  int tmp = do_log_interpolate_;
  do_log_interpolate_ = new_array.do_log_interpolate_;
  new_array.do_log_interpolate_ = tmp;

  // swap the dels
  double del_tmp = del_;
  del_ = new_array.del_;
  new_array.del_ = del_tmp;
  
  // swap locate array types
  LAType locate_tmp = locate_type_;
  locate_type_ = new_array.locate_type_;
  new_array.locate_type_ = locate_tmp;
}

void locate_array::scale(double e) {
  for (int i = 0; i < size(); i++) {
    x_[i] *= e;
  }
  if (locate_type_ == do_lin) {
    del_ *= e;
  }
  min_ *= e;
}

