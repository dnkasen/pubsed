#include "physical_constants.h"
#include "nlte_gas.h"
#include <iostream>
#include <mpi.h>

namespace pc = physical_constants;

nlte_gas::nlte_gas()
{
  e_gamma = 0;
  no_ground_recomb = 0;
}

//-----------------------------------------------------------
// set up all of the atoms making up the gas
//-----------------------------------------------------------
void nlte_gas::init
(std::string datadir, std::vector<int> e, std::vector<int> A, locate_array ng)
{
  // verbosity
  int my_rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
  verbose = (my_rank == 0);
  
  for (int i=0;i<e.size();i++) elem_Z.push_back(e[i]);
  for (int i=0;i<e.size();i++) elem_A.push_back(A[i]);
  mass_frac.resize(e.size());

  // copy the nugrid
  nu_grid.copy(ng);

  // read in the atom data
  atoms.resize(elem_Z.size());
  for (int i=0;i<atoms.size();i++) 
  {
    std::string fname = datadir + "/" + std::to_string(elem_Z[i]) + "/";
    fname = fname + std::to_string(elem_Z[i]) + ".comb";
    int success = atoms[i].Init(fname, ng); 
    if (!success) 
      if (verbose) std::cout << "# Can't open data file " << fname << "\n";
  }
}


//-----------------------------------------------------------
// read fuzz lines from a file
// for now, this is a ascii file; should really be hdf5
//-----------------------------------------------------------
int nlte_gas::read_fuzzfile(std::string fuzzfile)
{
  // make a map of the elements
  int map[100];
  for (int i=0;i<100;i++) map[i] = -1;
  for (int i=0;i<elem_Z.size();i++) map[atoms[i].Z] = i;

  // open up the fuzz file
  FILE *in = fopen(fuzzfile.c_str(),"r");
  if (in == NULL) return -1;

  int n_fuzz, n_store = 0;
  fscanf(in,"%d\n",&n_fuzz);

  double l,g,e;
  int i,el,ion;
  for (i=0;i<n_fuzz;i++)
  {
    fscanf(in,"%lf %lf %lf %d %d\n",&l,&g,&e,&el,&ion);
    double nu = pc::c/(l*pc::angs_to_cm);
    // find element in the gas, if it exists
    int ind = map[el];
    if (ind < 0) continue;
    // if ionization state and within wavelength bounds, go
    if ((ion < atoms[ind].n_ions)&&(nu >= nu_grid.minval())
	&&(nu <= nu_grid.maxval()))
    {
      fuzz_line f;
      f.nu  = nu;
      f.gf  = g;
      f.El  = e;
      f.bin = nu_grid.locate(nu);
      f.ion = ion;
      atoms[ind].fuzz_lines.push_back(f);
      n_store++;
      std::cout << l << "\t" << e << "\t" << ion << std::endl;
    }
  }

  return n_store;

}

//-----------------------------------------------------------
// Solve for the gas state (excitation/ionization)
//  either in LTE or NLTE
//-----------------------------------------------------------
void nlte_gas::solve_state(int lte)
{

  for (int i=0;i<atoms.size();i++)
  {
    atoms[i].n_dens  = dens*mass_frac[i]/(elem_A[i]*pc::m_p);
    atoms[i].e_gamma = e_gamma*mass_frac[i];
    atoms[i].no_ground_recomb = no_ground_recomb;
  }

  double max_ne = 10*dens/(A_mu*pc::m_p);
  double min_ne = 1e-6;
  double tol    = 1e-3;
  ne = ne_brent_method(min_ne,max_ne,tol,lte);
}

 


//-----------------------------------------------------------
// This is the function that represents the non-linear equation
// for N_e.  It is used in Brents method below to solve
// for the root, thus determining N_e.  This equation is
// basically just the one for charge conservation.
//-----------------------------------------------------------
double nlte_gas::charge_conservation(double ne, int lte)
{
  // start with charge conservation function f set to zero
  double f  = 0;
  // loop over all atoms
  for (int i=0;i<atoms.size();i++) 
  {
    // Solve the LTE or NLTE with this value of Ne
    if (lte) atoms[i].solve_lte (temp,ne,time);
    else     atoms[i].solve_nlte(temp,ne,time);

    // total electron donation from this atomic species
    f += dens*mass_frac[i]/(elem_A[i]*pc::m_p)*atoms[i].get_ion_frac();
  }

  //  printf("%e %e\n",f,ne);
  // total ionization density minus electron density should equal zero
  f = f - ne;
  // return to non-linear solver to iterate this to zero
  return f;
}


//-----------------------------------------------------------
// Brents method from Numerical Recipes to solve non-linear
// equation for electron density ne
//-----------------------------------------------------------
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
double nlte_gas::ne_brent_method(double x1,double x2,double tol,int lte)
{  
  int ITMAX = 100;
  double EPS = 3.0e-8;
  int iter;
  double a=x1,b=x2,c=x2,d,e,min1,min2;
  double fa=charge_conservation(a,lte);
  double fb=charge_conservation(b,lte);
  double fc,p,q,r,s,tol1,xm;
  
  //if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
  //printf("Root must be bracketed in zbrent\n");
  fc=fb;
  for (iter=1;iter<=ITMAX;iter++) {
    if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
      c=a;
      fc=fa;
      e=d=b-a;
    }
    if (fabs(fc) < fabs(fb)) {
      a=b;
      b=c;
      c=a;
      fa=fb;
      fb=fc;
      fc=fa;
    }
    tol1=2.0*EPS*fabs(b)+0.5*tol;
    xm=0.5*(c-b);
    if (fabs(xm) <= tol1 || fb == 0.0) return b;
    if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
      s=fb/fa;
      if (a == c) {
	p=2.0*xm*s;
	q=1.0-s;
      } else {
	q=fa/fc;
	r=fb/fc;
	p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
	q=(q-1.0)*(r-1.0)*(s-1.0);
      }
      if (p > 0.0) q = -q;
      p=fabs(p);
      min1=3.0*xm*q-fabs(tol1*q);
      min2=fabs(e*q);
      if (2.0*p < (min1 < min2 ? min1 : min2)) {
	e=d;
	d=p/q;
      } else {
	d=xm;
	e=d;
      }
    } else {
      d=xm;
      e=d;
    }
    a=b;
    fa=fb;
    if (fabs(d) > tol1)
      b += d;
    else
      b += SIGN(tol1,xm);
    fb=charge_conservation(b,lte);
  }
  printf("Maximum number of iterations exceeded in zbrent");
  return 0.0;
}



//-----------------------------------------------------------
// print out of the gas class
//-----------------------------------------------------------
void nlte_gas::print()
{
  std::cout << "# dens = " << dens << "\n";
  std::cout << "# temp = " << temp << "\n";
  std::cout << "# A_mu = " << A_mu << "\n";
  for (int i=0;i<elem_Z.size();i++)
    printf("%4d %12.4e %12.4e\n",elem_Z[i],mass_frac[i],atoms[i].get_ion_frac());
  for (int i=0;i<elem_Z.size();i++)
    atoms[i].print();
}
 
