#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <limits>

#include "hydro_1D_movingmesh.h"
#include "physical_constants.h"
namespace pc = physical_constants;


void hydro_1D_movingmesh::init(ParameterReader *params, grid_general *g)
{
  grid = g;
  theDomain = new Hydro1DMovingMeshDomain[1];

  theDomain->t      = grid->t_now;
  theDomain->t_init = grid->t_now;

  // also setup parameter file

  // final time and outputs -- try to eliminate these
  //theDomain->t_fin   = theDomain->theParList.t_max;
  //theDomain->N_rpt = theDomain->theParList.NumRepts;
  //theDomain->N_snp = theDomain->theParList.NumSnaps;
  //theDomain->N_chk = theDomain->theParList.NumChecks;
  // theDomain.count_steps = 0;
  //theDomain->final_step = 0;
  //theDomain->nrpt=-1;
  //theDomain->nsnp=-1;
  //theDomain->nchk=-1;
  theDomain->rank = 0;
  theDomain->size = 1;

  // point mass
  //theDomain->point_mass = theDomain->theParList.grav_pointmass;

  theDomain->Ng   = HYDRO_1D_MOVINGMESH_NUM_G;
  theDomain->Nr   = grid->n_zones;
  theDomain->theCells = new Hydro1DMovingMeshCell[theDomain->Nr];



  // setup the cells
  for( int i=0 ; i<theDomain->Nr ; ++i )
  {
    // setup
    //theDomain.theCells[i].riph = rp;  // outer radius of zone
    //theDomain.theCells[i].dr   = rp - rm;

    Hydro1DMovingMeshCell *c = &(theDomain->theCells[i]);

    c->prim[RHO] = grid->z[i].rho;  //rho;    // rho
    c->prim[PPP] = 0.0; //P_min;  // pressure
    c->prim[VRR] = 0.0;   // radial velocity
    c->prim[XXX] = 0.0;   // passive scalar
    c->prim[AAA] = 0.0;   // Rayliegh Taylor

    double dV = grid->zone_volume(i);
    prim2cons( c->prim , c->cons , dV );
    cons2prim( c->cons , c->prim , dV );
  }

  // set from parameter file
  //GAMMA_LAW = theDomain->theParList.Adiabatic_Index;
  //RHO_FLOOR = theDomain->theParList.Density_Floor;
  //PRE_FLOOR = theDomain->theParList.Pressure_Floor;
  //USE_RT = theDomain->theParList.rt_flag;   // Rayleigh taylor

  //
  //riemann_solver = theDomain->theParList.Riemann_Solver;
  //rt_flag = theDomain->theParList.rt_flag;
  //gamma_law = theDomain->theParList.Adiabatic_Index;

}


//------------------------------------------------------------
// Convert primitive hydro variables to conserverative
//------------------------------------------------------------
void hydro_1D_movingmesh::prim2cons
( double * prim , double * cons , double dV )
{
   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double vr  = prim[VRR];
   double v2 = vr*vr;
   double gam = GAMMA_LAW;
   double rhoe = Pp/(gam-1.);

   cons[DDD] = rho*dV;
   cons[SRR] = rho*vr*dV;
   cons[TAU] = (.5*rho*v2 + rhoe)*dV;

   for(int q=XXX ; q<HYDRO_1D_MOVINGMESH_NUM_Q; ++q )
      cons[q] = cons[DDD]*prim[q];

}

//------------------------------------------------------------
// Convert conservativ hydro variables to primitive
//------------------------------------------------------------
void hydro_1D_movingmesh::cons2prim
( double * cons , double * prim , double dV )
{

   double rho = cons[DDD]/dV;
   double Sr  = cons[SRR]/dV;
   double E   = cons[TAU]/dV;

   double vr = Sr/rho;
   double v2 = vr*vr;
   double rhoe = E - .5*rho*v2;
   double gam = GAMMA_LAW;
   double Pp = (gam-1.)*rhoe;

   if( rho<RHO_FLOOR ) rho=RHO_FLOOR;
   if( Pp < PRE_FLOOR*rho ) Pp = PRE_FLOOR*rho;

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[VRR] = vr;

   int q;
   for( q=XXX ; q<HYDRO_1D_MOVINGMESH_NUM_Q; ++q ){
      prim[q] = cons[q]/cons[DDD];
   }

}

//------------------------------------------------------------
// Calclulate hydro time step
//------------------------------------------------------------
double hydro_1D_movingmesh::get_time_step()
{
  set_wcell();
  double dt = getmindt();
  return dt;
}

//------------------------------------------------------------
//
//------------------------------------------------------------
void hydro_1D_movingmesh::set_wcell()
{

   struct Hydro1DMovingMeshCell * theCells = theDomain->theCells;
   int mesh_motion = theDomain->theParList.Mesh_Motion;
   int Nr = theDomain->Nr;

   int i;
   for( i=0 ; i<Nr-1 ; ++i ){
      struct Hydro1DMovingMeshCell * cL = theCells+i;
      double w = 0.0;
      if( mesh_motion ){
         struct Hydro1DMovingMeshCell * cR = theCells+i+1;
         double wL = get_vr( cL->prim );
         double wR = get_vr( cR->prim );
         w = .5*(wL + wR);
         if( i==0 && theDomain->rank==0 ) w = 0.5*wR*(cL->riph)/(cR->riph-.5*cR->dr);//*(cR->riph - .5*cR->dr)/(cL->riph);//0.0;//2./3.*wR;
         //if( i==0 && theDomain->rank==0 ) w = 0.0;
      }
      cL->wiph = w;
   }
}


//------------------------------------------------------------
//
//------------------------------------------------------------
double hydro_1D_movingmesh::getmindt()
{
     struct Hydro1DMovingMeshCell* theCells = theDomain->theCells;
     int Nr = theDomain->Nr;

     double dt = 1e100;
     int i;
     for( i=1 ; i<Nr-1 ; ++i ){
        int im = i-1;
        struct Hydro1DMovingMeshCell* c = theCells+i;
        double dr = c->dr;
        double r = c->riph-.5*dr;
        double wm = theCells[im].wiph;
        double wp = c->wiph;
        double w = .5*(wm+wp);
        double dt_temp = mindt( c->prim , w , r , dr );
        if( dt > dt_temp ) dt = dt_temp;
     }
     dt *= theDomain->theParList.CFL;

     return( dt );
  }

//------------------------------------------------------------
//
//------------------------------------------------------------
double hydro_1D_movingmesh::mindt
( double * prim , double w , double r , double dr )
{
   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double vr  = prim[VRR];
   double gam = GAMMA_LAW;

   double cs = sqrt(fabs(gam*Pp/rho));
   double eta = get_eta( prim , NULL , r );

   double maxvr = cs + fabs( vr - w );
   double dt = dr/maxvr;
   double dt_eta = dr*dr/eta;
   if( dt > dt_eta && USE_RT ) dt = dt_eta;

   return( dt );
}

//------------------------------------------------------------
//
//------------------------------------------------------------
double hydro_1D_movingmesh::get_eta
( double * prim , double * grad_prim , double r )
{
   double C = 0.06*1.7;//0.03;
   double gam = GAMMA_LAW;

   double cs = sqrt( gam*fabs(prim[PPP]/prim[RHO]) );

   double alpha = prim[AAA];
   if( alpha < 0.0 ) alpha = 0.0;

   double u_eddy = cs*sqrt( alpha );
   double lambda = r*sqrt(alpha);
   double eta = C*u_eddy*lambda;

   return( eta );

}


void hydro_1D_movingmesh::step(double dt)
{


}
