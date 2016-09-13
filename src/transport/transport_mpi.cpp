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

  // eventually do a smarter reduction
  int ng = nu_grid.size();
  int nz = grid->n_zones;
  int blocksize = 100000; // transfer around 10 million # per round
  int nz_per_block = floor(blocksize/ng);
  if (nz_per_block > nz) nz_per_block = nz;
  if (nz_per_block < 1)  nz_per_block  = 1;

  // new block size
  blocksize = ng;
  double *src = new double[blocksize];
  double *dst = new double[blocksize];
  for (int i=0;i<nz;i++)
  {
    //-----------------------------
    // absorptive opacity
    //-----------------------------
    for (int j=0;j<blocksize;j++)
    {
      src[j] = abs_opacity_[i][j];
      dst[j] = 0.0;
    }
    MPI_Allreduce(src,dst,blocksize,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    for (int j=0;j<blocksize;j++)
      abs_opacity_[i][j] = dst[j];

    //-----------------------------
    // scattering opacity
    //-----------------------------
    for (int j=0;j<blocksize;j++)
    {
      src[j] = scat_opacity_[i][j];
      dst[j] = 0.0;
    }
    MPI_Allreduce(src,dst,blocksize,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    for (int j=0;j<blocksize;j++)
      scat_opacity_[i][j] = dst[j];

    //-----------------------------
    // emissivity
    //-----------------------------
    for (int j=0;j<blocksize;j++)
    {
      src[j] = emissivity_[i].get_value(j);
      dst[j] = 0.0;
    }
    MPI_Allreduce(src,dst,blocksize,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    for (int j=0;j<blocksize;j++)
      emissivity_[i].set_value(j,dst[j]);

    // normalize emissivity cdf
    emissivity_[i].normalize();
  }
  delete[] src;
  delete[] dst;

  //=************************************************
  // do zone scalars
  //=************************************************
  src = new double[nz];
  dst = new double[nz];
  for (int i=0;i<nz;i++)
  {
    src[i] = compton_opac[i];
    dst[i] = 0.0;
  }
  MPI_Allreduce(src,dst,nz,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  for (int i=0;i<nz;i++) compton_opac[i] = dst[i];

  for (int i=0;i<nz;i++)
  {
    src[i] = photoion_opac[i];
    dst[i] = 0.0;
  }
  MPI_Allreduce(src,dst,nz,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  for (int i=0;i<nz;i++) photoion_opac[i] = dst[i];

  delete[] src;
  delete[] dst;
}

//------------------------------------------------------------
// Combine the radiation tallies in all zones
// from all processors using MPI 
//------------------------------------------------------------
 void transport::reduce_radiation(double dt)
{
  // eventually do a smarter reduction
  int ng = nu_grid.size();
  int nz = grid->n_zones;
  int blocksize = 100000; // transfer around 10 million # per round
  int nz_per_block = floor(blocksize/ng);
  if (nz_per_block > nz) nz_per_block = nz;
  if (nz_per_block < 1)  nz_per_block  = 1;

  // new block size
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

   //=************************************************
  // do zone scalars
  //=************************************************
  src = new double[nz];
  dst = new double[nz];
  for (int i=0;i<nz;i++)
  {
    src[i] = grid->z[i].e_abs;
    dst[i] = 0.0;
  }
  MPI_Allreduce(src,dst,nz,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  for (int i=0;i<nz;i++) grid->z[i].e_abs = dst[i]/MPI_nprocs;
  for (int i=0;i<nz;i++)
  {
    src[i] = grid->z[i].L_radio_dep;
    dst[i] = 0.0;
  }
  MPI_Allreduce(src,dst,nz,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  for (int i=0;i<nz;i++) grid->z[i].L_radio_dep = dst[i]/MPI_nprocs;


  delete[] src;
  delete[] dst;

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

    if (nu_grid.size() == 1)
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
