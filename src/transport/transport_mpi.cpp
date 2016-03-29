#include <math.h>
#include <cassert>
#include "transport.h"
#include "physical_constants.h"

namespace pc = physical_constants;

//------------------------------------------------------------
// Combine the radiation tallies in all zones
// from all processors using MPI 
//------------------------------------------------------------
 void transport::reduce_radiation(double dt)
{
  // properly normalize the radiative quantities
  for (int i=0;i<grid->n_zones;i++) 
  {
    double vol = grid->zone_volume(i);
    grid->z[i].e_rad   /= vol*pc::c*dt;
    grid->z[i].e_abs   /= vol*dt; 
    //grid->z[i].fx_rad  /= vol*pc::c*dt; 
    //grid->z[i].fy_rad  /= vol*pc::c*dt;
    //grid->z[i].fz_rad  /= vol*pc::c*dt;
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
}
