#include <math.h>
#include <cassert>
#include "transport.h"
#include "physical_constants.h"

namespace pc = physical_constants;


//------------------------------------------------------------
// Clear radiation quantities
//------------------------------------------------------------
void transport::wipe_radiation()
{
  for (int i=0;i<grid->n_zones;i++) 
  {
    grid->z[i].e_rad  = 0;
    grid->z[i].e_abs  = 0;
    grid->z[i].L_radio_dep = 0;
    if (store_Jnu_)
      for (int j=0;j<nu_grid.size();j++)
        J_nu_[i][j] = 0;
    grid->z[i].L_radio_emit = 0;
    //grid->z[i].fx_rad = 0;
    //grid->z[i].fy_rad = 0;
    //grid->z[i].fz_rad = 0;
  }
}

//------------------------------------------------------------
// Combine the opacity calculations in all zones
// from all processors using MPI 
//------------------------------------------------------------
void transport::reduce_opacities()
{
  //=************************************************
  // do zone vectors
  //=************************************************
  if (MPI_nprocs == 1) return;

  // dimensions
  int nw = nu_grid.size();
  int nz = grid->n_zones;

  // maximum size of transfer blocks
  int max_blocksize = Max_MPI_Blocksize;
  if (nw > Max_MPI_Blocksize) {
    std::cout << "Error, frequency grid is bigger than MPI_Max_Blocksize\n";
    exit(1);
  }

  // number of zones that fit into a transfer block
  int nz_per_block      = floor(1.0*max_blocksize/nw);
  // actual blocksize and # of blocks
  int blocksize         = nz_per_block*nw;
  int n_blocks          = floor(1.0*nw*nz/blocksize);
  // size of last block to pick up remainder
  int last_blocksize    = nw*nz - n_blocks*blocksize;
  int last_nz_per_block = nz - nz_per_block*n_blocks;

  // sanity check
  if (blocksize > Max_MPI_Blocksize)
  {
    std::cout << "Error, Blocksize greater than MPI_Max_Blocksize\n";
    exit(1);
  }

  int cnt;
  //-----------------------------
  // loop over blocks
  //-----------------------------
  for (int i=0;i<n_blocks+1;i++)
  {
    int this_nz        = nz_per_block;
    int this_blocksize = blocksize;
    if (i == n_blocks) {
      this_nz = last_nz_per_block;
      this_blocksize = last_blocksize; }

    //-----------------------------
    // absorptive opacity
    //-----------------------------
    cnt = 0;
    for (int j=0;j<this_nz;j++)
    {
      int iz = i*nz_per_block + j; 
      for (int k=0;k<nw;k++)
      {
        src_MPI_block[cnt] = abs_opacity_[iz][k];
        dst_MPI_block[cnt] = 0.0;
        cnt++;
      }
    }
    MPI_Allreduce(src_MPI_block,dst_MPI_block,this_blocksize,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    cnt = 0;
    for (int j=0;j<this_nz;j++)
    {
      int iz = i*nz_per_block + j; 
      for (int k=0;k<nw;k++)
      {
        abs_opacity_[iz][k] = (OpacityType)dst_MPI_block[cnt];
        cnt++;
      }
    }

    //-----------------------------
    // scattering opacity
    //-----------------------------
    if (!omit_scattering_)
    {
      cnt = 0;
      for (int j=0;j<this_nz;j++)
      {
        int iz = i*nz_per_block + j; 
        for (int k=0;k<nw;k++)
        {
          src_MPI_block[cnt] = scat_opacity_[iz][k];
          dst_MPI_block[cnt] = 0.0;
          cnt++;
        }
      }
      MPI_Allreduce(src_MPI_block,dst_MPI_block,this_blocksize,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      cnt = 0;
      for (int j=0;j<this_nz;j++)
      {
        int iz = i*nz_per_block + j; 
        for (int k=0;k<nw;k++)
        {  
          scat_opacity_[iz][k] = (OpacityType)dst_MPI_block[cnt];
          cnt++;
        }
      }
    }
    //-----------------------------
    // emissivity
    //-----------------------------
    cnt = 0;
    for (int j=0;j<this_nz;j++)
    {
      int iz = i*nz_per_block + j; 
      for (int k=0;k<nw;k++)
      {
        src_MPI_block[cnt] = emissivity_[iz].get(k);
        dst_MPI_block[cnt] = 0.0;
        cnt++;
      }
    }
    MPI_Allreduce(src_MPI_block,dst_MPI_block,this_blocksize,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    cnt = 0;
    for (int j=0;j<this_nz;j++)
    {
      int iz = i*nz_per_block + j; 
      for (int k=0;k<nw;k++)
      {
        emissivity_[iz].set(k,(OpacityType)dst_MPI_block[cnt]);
        cnt++;
      }
    }
  }

  //=************************************************
  // do zone scalars
  //=************************************************
  for (int i=0;i<nz;i++)
  {
    src_MPI_zones[i] = compton_opac[i];
    dst_MPI_zones[i] = 0.0;
  }
  MPI_Allreduce(src_MPI_zones,dst_MPI_zones,nz,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  for (int i=0;i<nz;i++) compton_opac[i] = dst_MPI_zones[i];

  for (int i=0;i<nz;i++)
  {
    src_MPI_zones[i] = photoion_opac[i];
    dst_MPI_zones[i] = 0.0;
  }
  MPI_Allreduce(src_MPI_zones,dst_MPI_zones,nz,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  for (int i=0;i<nz;i++) photoion_opac[i] = dst_MPI_zones[i];
}

//------------------------------------------------------------
// Combine the solved for temperature in zones
// from all processors using MPI 
//------------------------------------------------------------
 void transport::reduce_Tgas()
 {
  if (MPI_nprocs == 1) return;

  //=************************************************
  // do zone scalar
  //=************************************************
  int nz = grid->n_zones;
  for (int i=0;i<nz;i++)
  {
    src_MPI_zones[i] = 0;
    dst_MPI_zones[i] = 0.0;
  }
  for (int i=my_zone_start_;i<my_zone_stop_;i++)
    src_MPI_zones[i] = grid->z[i].T_gas;

  MPI_Allreduce(src_MPI_zones,dst_MPI_zones,nz,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  for (int i=0;i<nz;i++) grid->z[i].T_gas = dst_MPI_zones[i];
 }

//------------------------------------------------------------
// Combine the radiation tallies in all zones
// from all processors using MPI 
//------------------------------------------------------------
 void transport::reduce_radiation(double dt)
{
  if (MPI_nprocs > 1)
  {
  // eventually do a smarter reduction
    int ng = nu_grid.size();
    int nz = grid->n_zones;
    int blocksize = 100000; // transfer around 10 million # per round
    int nz_per_block = floor(blocksize/ng);
    if (nz_per_block > nz) nz_per_block = nz;
    if (nz_per_block < 1)  nz_per_block  = 1;

    // new block size
    if (store_Jnu_)
    {
      blocksize = ng;
      double *src = new double[blocksize];
      double *dst = new double[blocksize];
      for (int i=0;i<nz;i++)
      {
        for (int j=0;j<blocksize;j++)
        {
        src[j] = J_nu_[i][j];
        dst[j] = 0.0;
        }
        MPI_Allreduce(src,dst,blocksize,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        for (int j=0;j<blocksize;j++)
          J_nu_[i][j] = dst[j]/MPI_nprocs;
      }
      delete[] src;
      delete[] dst;
    } 

     //=************************************************
    // do zone scalars
    //=************************************************
    for (int i=0;i<nz;i++)
    {
      src_MPI_zones[i] = grid->z[i].e_abs;
      dst_MPI_zones[i] = 0.0;
    }
    MPI_Allreduce(src_MPI_zones,dst_MPI_zones,nz,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    for (int i=0;i<nz;i++) grid->z[i].e_abs = dst_MPI_zones[i]/MPI_nprocs;
    for (int i=0;i<nz;i++)
    {
      src_MPI_zones[i] = grid->z[i].L_radio_dep;
      dst_MPI_zones[i] = 0.0;
    }
    MPI_Allreduce(src_MPI_zones,dst_MPI_zones,nz,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    for (int i=0;i<nz;i++) grid->z[i].L_radio_dep = dst_MPI_zones[i]/MPI_nprocs;
  }

  //=************************************************
  // properly normalize the radiative quantities
  //=************************************************
  for (int i=0;i<grid->n_zones;i++) 
  {
    double vol = grid->zone_volume(i);
    grid->z[i].e_abs   /= vol*dt;
    grid->z[i].L_radio_dep /= vol*dt;
    grid->z[i].L_radio_emit /= vol;

    //grid->z[i].fx_rad  /= vol*pc::c*dt; 
    //grid->z[i].fy_rad  /= vol*pc::c*dt;
    //grid->z[i].fz_rad  /= vol*pc::c*dt;

    if ((nu_grid.size() == 1)||(!store_Jnu_))
    {
      grid->z[i].e_rad = J_nu_[i][0]/(vol*dt*pc::c);
    }
    else 
    {
      double esum = 0;
      for (int j=0;j<nu_grid.size();j++) 
      {
        J_nu_[i][j] /= vol*dt*4*pc::pi*nu_grid.delta(j);
        esum += J_nu_[i][j]*nu_grid.delta(j)*4*pc::pi/pc::c; 
      }
      grid->z[i].e_rad = esum;
    }
  }
}





//   vector<real> send, receive;
//   int my_begin, my_end, size;

//   //-- EACH PROCESSOR GETS THE REDUCTION INFORMATION IT NEEDS
//   for(int proc=0; proc<MPI_nprocs; proc++){

//     // set the begin and end indices so a process covers range [begin,end)
//     my_begin = ( proc==0 ? 0 : my_zone_end[proc-1] );
//     my_end = my_zone_end[proc];

//     // set the computation size and create the send/receive vectors
//     size = my_end - my_begin;
//     send.resize(size);
//     receive.resize(size);

//     // reduce e_rad
//     for(int i=my_begin; i<my_end; i++) send[i-my_begin] = grid->z[i].e_abs;
//     MPI_Reduce(&send.front(), &receive.front(), size, MPI_real, MPI_SUM, proc, MPI_COMM_WORLD);
//     for(int i=my_begin; i<my_end; i++) grid->z[i].e_abs = receive[i-my_begin] / (real)MPI_nprocs;
//     }

//     // TODO - need to put in other quantities...
//   }
