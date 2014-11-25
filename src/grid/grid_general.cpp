#include <mpi.h>
#include <math.h>

#include "grid_general.h"



//------------------------------------------------------------
// initialize the grid
//------------------------------------------------------------
void grid_general::init(ParameterReader* params)
{
  // read the model file or fill in custom model
  read_model_file(params);

  // complain if the grid is obviously not right
  if(z.size()==0){
    std::cout << "Error: there are no grid zones." << std::endl;
    exit(5);

    n_zones = z.size();
  }
}



//------------------------------------------------------------
// Combine the radiation tallies in all zones
// from all processors using MPI allreduce
//------------------------------------------------------------
void grid_general::reduce_radiation()
{
  // largest block to reduce
  int bsize    = 10000;
  int n_blocks = floor(n_zones/bsize);

  // reduce in blocks
  for (int i=0;i<n_blocks;i++) 
    reduce_radiation_block(bsize,i*bsize);
  
  // get remainder
  int remainder = n_zones - n_blocks*bsize;
  if (remainder > 0) reduce_radiation_block(remainder,n_blocks*bsize);
}


//------------------------------------------------------------
// Combine the radiation tallies in blocks
// from all processors  using MPI allreduce
//------------------------------------------------------------
void grid_general::reduce_radiation_block(int bsize, int start)
{
  int j,size;
  MPI_Comm_size( MPI_COMM_WORLD, &size );

  double *src = new double[bsize];
  double *dst = new double[bsize];
  
  // reduce (average) e_rad
  for (j=0;j<bsize;j++) 
  {
    src[j] = z[start + j].e_rad/size;
    dst[j] = 0;
  }
  MPI_Allreduce(src,dst,bsize,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  for (j=0;j<bsize;j++) z[start + j].e_rad = dst[j];


  // reduce (average) e_abs
  for (j=0;j<bsize;j++) 
  {
    src[j] = z[start + j].e_abs/size;
    dst[j] = 0;
  }
  MPI_Allreduce(src,dst,bsize,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  for (j=0;j<bsize;j++) z[start + j].e_abs = dst[j];
    
  // need to put in other quantities...

  delete src;
  delete dst;

}


void grid_general::reduce_T_gas()
{
  // reduce gas temperature 
  double *src_ptr = new double[n_zones];
  double *dst_ptr = new double[n_zones];

  for (int i=0;i<n_zones;i++) {src_ptr[i] = z[i].T_gas; dst_ptr[i] = 0.0;}
  MPI_Allreduce(src_ptr,dst_ptr,n_zones,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  for (int i=0;i<n_zones;i++) z[i].T_gas = dst_ptr[i];

  delete src_ptr;
  delete dst_ptr;

}
