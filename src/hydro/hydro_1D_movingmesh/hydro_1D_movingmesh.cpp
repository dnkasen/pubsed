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
    // NEED TO setup
    //theDomain.theCells[i].riph = rp;  // outer radius of zone
    //theDomain.theCells[i].dr   = rp - rm;

    // NEED to setup
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

  // SET save boundary condition values
  boundary_prim_[RHO] = 0.0;
  boundary_prim_[PPP] = 0.0;
  boundary_prim_[VRR] = 0.0;
  boundary_prim_[XXX] = 0.0;
  boundary_prim_[AAA] = 0.0;


  // NEED to setup
  //GAMMA_LAW = theDomain->theParList.Adiabatic_Index;
  //RHO_FLOOR = theDomain->theParList.Density_Floor;
  //PRE_FLOOR = theDomain->theParList.Pressure_Floor;
  //USE_RT = theDomain->theParList.rt_flag;   // Rayleigh taylor

  // NEED to setup
  //riemann_solver = theDomain->theParList.Riemann_Solver;
  //rt_flag = theDomain->theParList.rt_flag;
  //gamma_law = theDomain->theParList.Adiabatic_Index;
  riemann_solver = 1;
  rt_flag = 0;
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
// Take a single hydro time step
//------------------------------------------------------------
void hydro_1D_movingmesh::step(double dt )
{
  // NEED TO PUT data in place

   struct Hydro1DMovingMeshCell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;

   int i;

   for( i=0 ; i<Nr ; ++i ){
      struct Hydro1DMovingMeshCell * c = theCells+i;
      memcpy( c->RKcons , c->cons , HYDRO_1D_MOVINGMESH_NUM_Q*sizeof(double) );
   }

   take_onestep(0.0 ,     dt , 1 , 0 );
   take_onestep(0.5 , 0.5*dt , 0 , 1 );

//   onestep( theDomain , 0.0 ,     dt , 1 , 1 );

   theDomain->t += dt;
   theDomain->count_steps += 1;

   // NEED TO PUT data back and adjust grid faces

}

void hydro_1D_movingmesh::take_onestep
(double RK, double dt, int first_step, int last_step)
{
    adjust_RK_cons( RK );
    radial_flux(dt);
    add_source(dt);

    if( first_step ) move_cells(RK , dt );
    calc_dr();
    calc_prim();

    //if( last_step ){
    //   AMR( theDomain );
  //  }
    boundary();


}




//------------------------------------------------------------
// Set's velocity of cells
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
//
//------------------------------------------------------------
double hydro_1D_movingmesh::getmindt()
{
  struct Hydro1DMovingMeshCell* theCells = theDomain->theCells;
  int Nr = theDomain->Nr;

  double dt = 1e100;
  for( int i=1 ; i<Nr-1 ; ++i )
  {
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


//------------------------------------------------------------
//
//------------------------------------------------------------
void hydro_1D_movingmesh::adjust_RK_cons(double RK )
{
   struct Hydro1DMovingMeshCell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;

   int i,q;
   for( i=0 ; i<Nr ; ++i ){
      Hydro1DMovingMeshCell * c = theCells+i;
      for( q=0 ; q<HYDRO_1D_MOVINGMESH_NUM_Q ; ++q ){
         c->cons[q] = (1.-RK)*c->cons[q] + RK*c->RKcons[q];
      }
   }
}

//------------------------------------------------------------
//
//------------------------------------------------------------
void hydro_1D_movingmesh::radial_flux(double dt )
{
   struct Hydro1DMovingMeshCell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   plm();

   for( int i=0 ; i<Nr-1 ; ++i )
   {
      struct Hydro1DMovingMeshCell * cL = theCells+i;
      struct Hydro1DMovingMeshCell * cR = theCells+i+1;
      double r = cL->riph;
      double dA = get_dA(r);
      riemann( cL , cR , r , dA*dt );
   }
}


//------------------------------------------------------------
// Geometrical functions
//------------------------------------------------------------
double hydro_1D_movingmesh::get_dA( double r )
{
   return( 4.*M_PI*r*r );
}

double hydro_1D_movingmesh::get_dV( double rp , double rm )
{
   double dr  = rp-rm;
   double r2    = (rp*rp+rm*rm+rp*rm)/3.;
   return( 4.*M_PI*r2*dr );
}

double hydro_1D_movingmesh::get_moment_arm( double rp , double rm )
{
   double r3 = (rp*rp*rp + rp*rp*rm + rp*rm*rm + rm*rm*rm)/4.;
   double r2 = (rp*rp + rp*rm + rm*rm)/3.;
   double r = r3/r2;
   return( r );
}


void hydro_1D_movingmesh::move_cells(double RK , double dt)
{
   struct Hydro1DMovingMeshCell *theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   for(int i=0 ; i<Nr ; ++i )
   {
      Hydro1DMovingMeshCell * c = theCells+i;
      c->riph += c->wiph*dt;
   }

}


//------------------------------------------------------------
// Do Riemann solve
//------------------------------------------------------------
void hydro_1D_movingmesh::riemann
( struct Hydro1DMovingMeshCell * cL , struct Hydro1DMovingMeshCell * cR, double r , double dAdt )
{
   double primL[HYDRO_1D_MOVINGMESH_NUM_Q];
   double primR[HYDRO_1D_MOVINGMESH_NUM_Q];

   double drL = .5*cL->dr;
   double drR = .5*cR->dr;
   int q;

   for(q=0 ; q<HYDRO_1D_MOVINGMESH_NUM_Q ; ++q )
   {
      primL[q] = cL->prim[q] + cL->grad[q]*drL;
      primR[q] = cR->prim[q] - cR->grad[q]*drR;
   }

   double Sl,Sr,Ss;

   vel( primL , primR , &Sl , &Sr , &Ss );

   double Fl[HYDRO_1D_MOVINGMESH_NUM_Q];
   double Fr[HYDRO_1D_MOVINGMESH_NUM_Q];
   double Ul[HYDRO_1D_MOVINGMESH_NUM_Q];
   double Ur[HYDRO_1D_MOVINGMESH_NUM_Q];

   double Flux[HYDRO_1D_MOVINGMESH_NUM_Q];

   double w = cL->wiph;

   if( w < Sl )
   {
      flux( primL , Fl );
      prim2cons( primL , Ul , 1.0 );

      for( q=0 ; q<HYDRO_1D_MOVINGMESH_NUM_Q ; ++q ){
         Flux[q] = Fl[q] - w*Ul[q];
      }
   }
   else if( w > Sr )
   {
      flux( primR , Fr );
      prim2cons( primR , Ur , 1.0 );

      for( q=0 ; q<HYDRO_1D_MOVINGMESH_NUM_Q ; ++q ){
         Flux[q] = Fr[q] - w*Ur[q];
      }
   }
   else
   {
      if( riemann_solver == _HLL_ ){
         double Fstar;
         double Ustar;
         double aL =  Sr;
         double aR = -Sl;

         prim2cons( primL , Ul , 1.0 );
         prim2cons( primR , Ur , 1.0 );
         flux( primL , Fl );
         flux( primR , Fr );

         for( q=0 ; q<HYDRO_1D_MOVINGMESH_NUM_Q ; ++q )
         {
            Fstar = ( aL*Fl[q] + aR*Fr[q] + aL*aR*( Ul[q] - Ur[q] ) )/( aL + aR );
            Ustar = ( aR*Ul[q] + aL*Ur[q] + Fl[q] - Fr[q] )/( aL + aR );

            Flux[q] = Fstar - w*Ustar;
         }
    }
    else
    {
         double Ustar[HYDRO_1D_MOVINGMESH_NUM_Q];
         double Uk[HYDRO_1D_MOVINGMESH_NUM_Q];
         double Fk[HYDRO_1D_MOVINGMESH_NUM_Q];
         if( w < Ss ){
            prim2cons( primL , Uk , 1.0 );
            getUstar( primL , Ustar , Sl , Ss );
            flux( primL , Fk );

            for( q=0 ; q<HYDRO_1D_MOVINGMESH_NUM_Q ; ++q ){
               Flux[q] = Fk[q] + Sl*( Ustar[q] - Uk[q] ) - w*Ustar[q];
            }
         }else{
            prim2cons( primR , Uk , 1.0 );
            getUstar( primR , Ustar , Sr , Ss );
            flux( primR , Fk );

            for( q=0 ; q<HYDRO_1D_MOVINGMESH_NUM_Q ; ++q ){
               Flux[q] = Fk[q] + Sr*( Ustar[q] - Uk[q] ) - w*Ustar[q];
            }
         }
      }
   }

   if( rt_flag )
   {
      double prim[HYDRO_1D_MOVINGMESH_NUM_Q];
      double consL[HYDRO_1D_MOVINGMESH_NUM_Q];
      double consR[HYDRO_1D_MOVINGMESH_NUM_Q];
      prim2cons( cL->prim , consL , 1.0 );
      prim2cons( cR->prim , consR , 1.0 );
      double gprim[HYDRO_1D_MOVINGMESH_NUM_Q];
      double gcons[HYDRO_1D_MOVINGMESH_NUM_Q];
      for( q=0 ; q<HYDRO_1D_MOVINGMESH_NUM_Q ; ++q ){
         prim[q] = .5*(primL[q]+primR[q]);
         gprim[q] = (cR->prim[q] - cL->prim[q])/(drL+drR);
         gcons[q] = (consR[q] - consL[q])/(drL+drR);
      }
//If new model, gcons[q] = P^1/gamma*grad( cons/P^1/gamma ).
/*
      double Pgam  = pow( prim[PPP] , 1./gamma_law );
      double PgamL = pow( cL->prim[PPP] , 1./gamma_law );
      double PgamR = pow( cR->prim[PPP] , 1./gamma_law );
      for( q=0 ; q<NUM_Q ; ++q ){
         gcons[q] = Pgam*( consR[q]/PgamR - consL[q]/PgamL )/(drL+drR);
      }
*/
////////
      double eta = get_eta( prim , gprim , r );
      for( q=0 ; q<HYDRO_1D_MOVINGMESH_NUM_Q ; ++q ){
         Flux[q] += -eta*gcons[q];
      }
   }

   for( q=0 ; q<HYDRO_1D_MOVINGMESH_NUM_Q ; ++q ){
      cL->cons[q] -= Flux[q]*dAdt;
      cR->cons[q] += Flux[q]*dAdt;
   }

}


//------------------------------------------------------------
//
//------------------------------------------------------------
void hydro_1D_movingmesh::flux( double * prim , double * flux )
{

   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double vr  = prim[VRR];
   double v2  = vr*vr;
   double gam = GAMMA_LAW;
   double rhoe = Pp/(gam-1.);

   flux[DDD] = rho*vr;
   flux[SRR] = rho*vr*vr + Pp;
   flux[TAU] = (.5*rho*v2 + rhoe + Pp)*vr;

   int q;
   for( q=XXX ; q<HYDRO_1D_MOVINGMESH_NUM_Q ; ++q ){
      flux[q] = flux[DDD]*prim[q];
   }
}


//------------------------------------------------------------
//
//------------------------------------------------------------
void hydro_1D_movingmesh::plm()
{
   struct Hydro1DMovingMeshCell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   double PLM = theDomain->theParList.PLM;
   int i,q;
   for( i=0 ; i<Nr ; ++i ){
      int im = i-1;
      int ip = i+1;
      if( i==0 ) im = 0;
      if( i==Nr-1 ) ip = Nr-1;
      struct Hydro1DMovingMeshCell * c  = theCells+i;
      struct Hydro1DMovingMeshCell * cL = theCells+im;
      struct Hydro1DMovingMeshCell * cR = theCells+ip;
      double drL = cL->dr;
      double drC = c->dr;
      double drR = cR->dr;
      for( q=0 ; q<HYDRO_1D_MOVINGMESH_NUM_Q ; ++q ){
         double pL = cL->prim[q];
         double pC = c->prim[q];
         double pR = cR->prim[q];
         double sL = pC - pL;
         sL /= .5*( drC + drL );
         double sR = pR - pC;
         sR /= .5*( drR + drC );
         double sC = pR - pL;
         sC /= .5*( drL + drR ) + drC;
         c->grad[q] = minmod( PLM*sL , sC , PLM*sR );
      }
   }
}

//------------------------------------------------------------
// Calculate wave speeds
//------------------------------------------------------------
void hydro_1D_movingmesh::vel
(double * prim1 , double * prim2 , double * Sl , double * Sr , double * Ss)
{
   double gam = GAMMA_LAW;

   double P1   = prim1[PPP];
   double rho1 = prim1[RHO];
   double vn1  = prim1[VRR];

   double cs1 = sqrt(fabs(gam*P1/rho1));

   double P2   = prim2[PPP];
   double rho2 = prim2[RHO];
   double vn2  = prim2[VRR];

   double cs2 = sqrt(fabs(gam*P2/rho2));

   *Ss = ( P2 - P1 + rho1*vn1*(-cs1) - rho2*vn2*cs2 )/( rho1*(-cs1) - rho2*cs2 );

   *Sr =  cs1 + vn1;
   *Sl = -cs1 + vn1;

   if( *Sr <  cs2 + vn2 ) *Sr =  cs2 + vn2;
   if( *Sl > -cs2 + vn2 ) *Sl = -cs2 + vn2;

}


//------------------------------------------------------------
//
//------------------------------------------------------------
void hydro_1D_movingmesh::getUstar
(double * prim , double * Ustar , double Sk , double Ss )
{

   double rho = prim[RHO];
   double vr  = prim[VRR];
   double Pp  = prim[PPP];
   double v2  = vr*vr;

   double gam = GAMMA_LAW;

   double rhoe = Pp/(gam-1.);

   double rhostar = rho*(Sk - vr)/(Sk - Ss);
   double Pstar = Pp*(Ss - vr)/(Sk - Ss);
   double Us = rhoe*(Sk - vr)/(Sk - Ss);

   Ustar[DDD] = rhostar;
   Ustar[SRR] = rhostar*( Ss );
   Ustar[TAU] = .5*rhostar*v2 + Us + rhostar*Ss*(Ss - vr) + Pstar;

   int q;
   for( q=XXX ; q<HYDRO_1D_MOVINGMESH_NUM_Q ; ++q ){
      Ustar[q] = prim[q]*Ustar[DDD];
   }

}


//------------------------------------------------------------
//
//------------------------------------------------------------
double hydro_1D_movingmesh::minmod
( double a , double b , double c )
{
   double m = a;
   if( a*b < 0.0 ) m = 0.0;
   if( fabs(b) < fabs(m) ) m = b;
   if( b*c < 0.0 ) m = 0.0;
   if( fabs(c) < fabs(m) ) m = c;
   return(m);
}


//------------------------------------------------------------
//
//------------------------------------------------------------
void hydro_1D_movingmesh::add_source(double dt )
{

   struct Hydro1DMovingMeshCell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   double t = theDomain->t;
   double grad[HYDRO_1D_MOVINGMESH_NUM_Q];

   int i,q;
   for( i=0 ; i<Nr ; ++i ){
      struct Hydro1DMovingMeshCell * c = theCells+i;
      double rp = c->riph;
      double rm = rp-c->dr;
      double r = get_moment_arm(rp,rm);
      double dV = get_dV(rp,rm);
      source( c->prim , c->cons , rp , rm , dV*dt );
      //debug
      //source_nozz( c->prim , c->cons , rp , rm , t , dV*dt );
      int inside = i>0 && i<Nr-1;
      for( q=0 ; q<HYDRO_1D_MOVINGMESH_NUM_Q ; ++q ){
         if( inside ){
            struct Hydro1DMovingMeshCell * cp = theCells+i+1;
            struct Hydro1DMovingMeshCell * cm = theCells+i-1;
            double dR = .5*cp->dr + c->dr + .5*cm->dr;
            grad[q] = (cp->prim[q]-cm->prim[q])/dR;
         }else{
            grad[q] = 0.0;
         }
      }
      source_alpha( c->prim , c->cons , grad , r , dV*dt );
   }

   //debug
//   int gravity_flag = theDomain->theParList.grav_flag;
  // if( gravity_flag ) gravity_addsrc( theDomain , dt );

}



//------------------------------------------------------------
//
//------------------------------------------------------------
void hydro_1D_movingmesh::source
( double * prim , double * cons , double rp , double rm , double dVdt )
{
   double Pp  = prim[PPP];
   double r  = .5*(rp+rm);
   double r2 = (rp*rp+rm*rm+rp*rm)/3.;
   cons[SRR] += 2.*Pp*(r/r2)*dVdt;
}

//------------------------------------------------------------
//
//------------------------------------------------------------
void hydro_1D_movingmesh::source_alpha
( double * prim , double * cons , double * grad_prim , double r , double dVdt )
{

   double A = 2e-5/1.7; //2e-5;//1e-4;
   double B = 1.2;//0.9;
   double D = 0.0;

   double gam = GAMMA_LAW;
   double alpha = prim[AAA];

   double Pp = prim[PPP];
   double rho = prim[RHO];
   double P1 = grad_prim[PPP];
   double rho1 = grad_prim[RHO];

   double g2 = -P1*rho1;
   if( g2 < 0.0 ) g2 = 0.0;
   double cs = sqrt(gam*fabs(Pp/rho));

   cons[AAA] += ( (A+B*alpha)*sqrt(g2) - D*rho*alpha*cs/r )*dVdt;
   if( cons[AAA] < 0. ) cons[AAA] = 0.;

}


//------------------------------------------------------------
//
//------------------------------------------------------------
void hydro_1D_movingmesh::calc_dr( )
{
   struct Hydro1DMovingMeshCell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;

   int i;
   for( i=1 ; i<Nr ; ++i )
   {
      int im = i-1;
      double rm = theCells[im].riph;
      double rp = theCells[i ].riph;
      double dr = rp-rm;
      theCells[i].dr = dr;
   }
   if( theDomain->rank==0 ) theCells[0].dr = theCells[0].riph;
}

//------------------------------------------------------------
//
//------------------------------------------------------------
void hydro_1D_movingmesh::calc_prim()
{
   struct Hydro1DMovingMeshCell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;

   int i;
   for( i=0 ; i<Nr ; ++i )
   {
      struct Hydro1DMovingMeshCell *c = theCells+i;
      double rp = c->riph;
      double rm = rp-c->dr;
      double dV = get_dV( rp , rm );
      cons2prim( c->cons , c->prim , dV );
   }
}


//------------------------------------------------------------
//
//------------------------------------------------------------
void hydro_1D_movingmesh::boundary()
{

   struct Hydro1DMovingMeshCell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int rank = theDomain->rank;
   int size = theDomain->size;


   if( rank == size-1 )
   {
      struct Hydro1DMovingMeshCell * cB = theCells+Nr-1;
      double rp = cB->riph;
      double rm = rp-cB->dr;
      double r = get_moment_arm(rp,rm);
      //set boundary condition
      cB->prim[RHO] = boundary_prim_[RHO];
      cB->prim[PPP] = boundary_prim_[PPP];
      cB->prim[VRR] = boundary_prim_[VRR];
      cB->prim[XXX] = boundary_prim_[XXX];
      cB->prim[AAA] = boundary_prim_[AAA];
   }
/*
   if( rank == 0 ){
      struct cell * cB = theCells+0;
      struct cell * cP = theCells+1;
      cB->prim[RHO] = cP->prim[RHO];
      cB->prim[PPP] = cP->prim[PPP];
   }
*/

   int ABSORB_R0 = theDomain->theParList.Absorb_BC;
   if( ABSORB_R0 && rank==0 )
   {
      struct Hydro1DMovingMeshCell * c3 = theCells+1;
      struct Hydro1DMovingMeshCell * c4 = theCells;

      //struct cell * cB = theCells+0;
      //struct cell * cP = theCells+1;
      //if( cP->prim[VRR] < 0.0 ){
         //cB->prim[RHO] = cP->prim[RHO];
         //cB->prim[PPP] = cP->prim[PPP];
         //cB->prim[VRR] = cP->prim[VRR];
      //}
      //if( c3->prim[VRR] < 0.0 ){
      for( int q=0 ; q<HYDRO_1D_MOVINGMESH_NUM_Q ; ++q ){
        c4->prim[q] = c3->prim[q]; }
      //c4->prim[AAA] = 0.0;
      //if( c4->prim[VRR] > 0.0 ) c4->prim[VRR] *= .5;
      //c4->prim[VRR] = 0.0;
      //}
   }

}
