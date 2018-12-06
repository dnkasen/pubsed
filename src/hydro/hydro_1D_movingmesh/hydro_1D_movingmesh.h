#ifndef _HYDRO_1D_MOVINGMESH_H
#define _HYDRO_1D_MOVINGMESH_H

#include <vector>
#include <fstream>
#include <iostream>
#include "zone.h"
#include "hydro_general.h"

#define HYDRO_1D_MOVINGMESH_NUM_Q 8
#define HYDRO_1D_MOVINGMESH_NUM_G 2


enum{RHO,PPP,VRR,XXX,AAA};
enum{DDD,TAU,SRR};


//------------------------------------
// struct that holds hydro data for a
// single zone A.K.A. cell
//------------------------------------
struct Hydro1DMovingMeshCell
{
   double prim[HYDRO_1D_MOVINGMESH_NUM_Q];
   double cons[HYDRO_1D_MOVINGMESH_NUM_Q];
   double RKcons[HYDRO_1D_MOVINGMESH_NUM_Q];
   double grad[HYDRO_1D_MOVINGMESH_NUM_Q];
   double riph;
   double dr;
   double miph;
   double dm;
   double wiph;
   double pot;
};

//------------------------------------
// struct for holding hydro parameters
//------------------------------------
struct Hydro1DMovingMeshCellParams
{

   int Num_R;
   double t_min, t_max;
   double rmin,rmax;
   int NumRepts, NumSnaps, NumChecks;
   int Out_LogTime;

   int LogZoning;
   double LogRadius;
   int Mesh_Motion, Riemann_Solver;
   double MaxShort, MaxLong;
   int Absorb_BC, Initial_Regrid, rt_flag;

   int grav_flag, grow_flag;
   double grav_G, grav_pointmass;

   double CFL, PLM;
   double Density_Floor, Pressure_Floor;
   double Adiabatic_Index;

   double rt_A,rt_B,rt_C,rt_D;

};


//------------------------------------
// struct for holding everything
// you need for the hydro
//------------------------------------
struct Hydro1DMovingMeshDomain
{

   struct Hydro1DMovingMeshCell * theCells;

   int Nr;  // number of cells
   int Ng;  // number of ghost zones

   double point_mass;

   time_t Wallt_init; // timer
   int rank,size;

   struct Hydro1DMovingMeshCellParams theParList;

   double t;
   int count_steps;
   double t_init, t_fin;
   int nrpt;
   int N_rpt;
   int nsnp;
   int N_snp;
   int nchk;
   int N_chk;

   int final_step;

};

//------------------------------------
// class for doing hydro
//------------------------------------
class hydro_1D_movingmesh : public hydro_general
{

private:

  struct Hydro1DMovingMeshDomain *theDomain;

  double GAMMA_LAW;
  double RHO_FLOOR;
  double PRE_FLOOR;
  double USE_RT;

  void cons2prim( double * cons , double * prim , double dV );
  void prim2cons( double * prim , double * cons , double dV );
  void set_wcell();
  double mindt( double * prim , double w , double r , double dr );

  double get_vr( double * prim )
  {  return( prim[VRR] );}
  double getmindt();
  double get_eta( double * prim , double * grad_prim , double r );

public:

  void   init(ParameterReader*, grid_general*);
  double get_time_step();
  void   step(double dt);

};

#endif
