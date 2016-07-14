#include <cstdlib>
#include <algorithm>
#include "physical_constants.h"
#include "nlte_gas.h"
#include <iostream>
#include <fstream>
#include <mpi.h>


namespace pc = physical_constants;

//----------------------------------------------------------------
// simple constructor
//----------------------------------------------------------------
nlte_gas::nlte_gas()
{
  e_gamma = 0;
  no_ground_recomb = 0;
}

//----------------------------------------------------------------
// initialize the gas by specifying the atoms that will
// compose it, along with datafile and freq. array
// inputs:
// std::string atomfile: name of atom data file (in hdf5)
// std::vector<int> e:  vector of atomic numbers
// std::vector<int> A:  vector of atomic weights (in atomic units)
// locate_array ng:  locate_array giving the freq. array
//---------------------------------------------------------------
void nlte_gas::initialize
(std::string atomfile, std::vector<int> e, std::vector<int> A, locate_array ng)
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

  // check if atomfile is there
  std::ifstream afile(atomfile);
  if (!afile) 
  {
    if (verbose) 
      std::cout << "Can't open atom datafile " << atomfile << "; exiting\n"; 
    exit(1); 
  }
  afile.close();

  // read in the atom data
  atoms.resize(elem_Z.size());
  int level_id = 0;
  for (int i=0;i<atoms.size();i++) 
  {
    int error = atoms[i].initialize(atomfile, elem_Z[i],ng,level_id); 
    if ((error)&&(verbose))
      std::cout << "# ERROR: incomplete data for atom Z=" << elem_Z[i] <<
	" in file " << atomfile << "\n";
  }
  

  // set up global list of levels (i.e., levels of all atoms)
  level_id = 0;
  for (int i=0;i<atoms.size();i++) 
  {
    for (int j=0;j<atoms[i].n_levels;j++)
    {
      globalLevelList_atom_.push_back(i);
      globalLevelList_index_.push_back(j);
    }
  }
  
  // set up global line list (lines from all atoms)
  for (int i=0;i<atoms.size();i++) 
  {
    for (int j=0;j<atoms[i].n_lines;j++)
    {
      nlteGlobalLine line;
      line.f_osc   = atoms[i].lines[j].f_lu;
      line.nu      = atoms[i].lines[j].nu;
      line.ion_weight = elem_A[i];
      line.iAtom   = i;
      line.iIndex  = j;

      int ll = atoms[i].lines[j].ll;
      int lu = atoms[i].lines[j].lu;
      //int global_ll = atoms[i].levels[ll].globalID;
      //int global_lu = atoms[i].levels[lu].globalID;
      line.iLowerLevel = ll; //global_ll;
      line.iUpperLevel = lu; //global_lu;
      
      // statstical weights
      int g1 = atoms[i].levels[ll].g;
      int g2 = atoms[i].levels[lu].g;
      line.g1_over_g2 = (1.0*g1)/(1.0*g2);
      
      globalLineList_.push_back(line);
    }
  }
    
  // now sort the global line list
  bool compareLineNu(const nlteGlobalLine &a, const nlteGlobalLine &b);
  std::sort(globalLineList_.begin(), globalLineList_.end(), compareLineNu);

  // sanity check
 // std::vector <nlteGlobalLine>::iterator iLine;
 //  std::cout << "nlines = " << globalLineList_.size() << "\n";
 //  int i = 0;
 //  for (iLine = globalLineList_.begin();iLine < globalLineList_.end();iLine++)
 //  {
 //    std::cout << i << " ";
 //    i++;
 //    std::cout << iLine->nu << "\t" << iLine->iIndex << " ";
 //    std::cout << iLine->iLowerLevel << "\t" << iLine->iUpperLevel << "\n";
 //  }
}


//----------------------------------------------------------------
// this is a simple comparison function used to sort the 
// global linelist in the initialize() function.
//----------------------------------------------------------------
bool compareLineNu(const nlteGlobalLine &a, const nlteGlobalLine &b)
{
  return a.nu < b.nu;
}

//-----------------------------------------------------------------
// get the data for the global line data
// input: vectors that will be resized to get lien data
//   std::vector<double> nu         frequncy of lines
//   std::vector<double> f_osc,     oscillator strength of lines
//   std::vector<int> lowerLevel    global index of lower level
//   std::vector<int> upperLevel)   global index of upper level
//-----------------------------------------------------------------
void nlte_gas::get_global_line_data(std::vector<double> nu,
				    std::vector<double> f_osc,
				    std::vector<double> g1_over_g2,
				    std::vector<int> iLowerLevel,
				    std::vector<int> iUpperLevel)
{
  int n_lines =  globalLineList_.size();
  nu.resize(n_lines);
  f_osc.resize(n_lines);
  iLowerLevel.resize(n_lines);
  iUpperLevel.resize(n_lines);

  for (int i=0;i<n_lines;i++)
  {
    nu[i]          = globalLineList_[i].nu;
    f_osc[i]       = globalLineList_[i].f_osc;
    g1_over_g2[i]  = globalLineList_[i].g1_over_g2;
    iLowerLevel[i] = globalLineList_[i].iLowerLevel;
    iUpperLevel[i] = globalLineList_[i].iUpperLevel;
  }
}


//-----------------------------------------------------------------
// Set mass fractions of each element in the gas
// this function will enforce that the mass fractions are
// normalize (i.e., add up to 1)
//
// input: 
// std::vector<double> x: vector of mass fractions of each element
//-----------------------------------------------------------------
void nlte_gas::set_mass_fractions(std::vector<double> x)
{ 
  double norm = 0.0;
  for (int i=0;i<mass_frac.size();i++) 
  {
    mass_frac[i] = x[i]; 
    norm += x[i];
  }
    
  // make sure it is normalized
  for (int i=0;i<mass_frac.size();i++) mass_frac[i] /= norm;
  
  // find mean weight of gas
  this->A_mu = 0;
  for (int i=0;i<mass_frac.size();i++) 
    A_mu += mass_frac[i]*elem_A[i]; 
}


//-----------------------------------------------------------
// read fuzz lines from a file
// input:
// std::string fuzzfile: name of hdf5 file with fuzz data
//-----------------------------------------------------------
int nlte_gas::read_fuzzfile(std::string fuzzfile)
{
  int n_tot = 0;
 
  // check if fuzzfile exists
  FILE *fin = fopen(fuzzfile.c_str(),"r");
  if (fin == NULL) return 0;
    
  for (int i=0;i<atoms.size();i++) n_tot += atoms[i].read_fuzzfile(fuzzfile);
  return n_tot;
}


//-----------------------------------------------------------
// return the ionization state, i.e., electron density
// ovr total ion number density
//-----------------------------------------------------------
double nlte_gas::get_ionization_state()
{
  double ni = dens/(A_mu*pc::m_p);
  return this->ne/ni;
}

//-----------------------------------------------------------
// Solve for the gas state (excitation/ionization)
// the level populations will be stored internally for
// further calculations
// input:
// int lte: 1 = do it in LTE, 0 = do it in NLTE
//-----------------------------------------------------------
int nlte_gas::solve_state(std::vector<real> J_nu)
{

  for (int i=0;i<atoms.size();i++)
  {
    atoms[i].n_dens  = dens*mass_frac[i]/(elem_A[i]*pc::m_p);
    atoms[i].e_gamma = e_gamma*mass_frac[i];
    atoms[i].no_ground_recomb = no_ground_recomb;
    // line widths
    double vd = sqrt(2*pc::k*temp/pc::m_p/elem_A[i]);
    if (line_velocity_width_ > 0) vd = line_velocity_width_;
    atoms[i].line_beta_dop_ = vd/pc::c;
  }

  double max_ne = 10*dens/(A_mu*pc::m_p);
  double min_ne = 1e-20*dens/(A_mu*pc::m_p);;
  double tol    = 1e-3;
  solve_error_  = 0;
  ne = ne_brent_method(min_ne,max_ne,tol,J_nu);

  return solve_error_;

}


//-----------------------------------------------------------
// This is the function that represents the non-linear equation
// for N_e.  It is used in Brents method below to solve
// for the root, thus determining N_e.  This equation is
// basically just the one for charge conservation.
//-----------------------------------------------------------
double nlte_gas::charge_conservation(double ne,std::vector<real> J_nu)
{
  // start with charge conservation function f set to zero
  double f  = 0;
  // loop over all atoms
  for (int i=0;i<atoms.size();i++) 
  {
    // Solve the LTE or NLTE with this value of Ne
    if (use_nlte_) atoms[i].solve_nlte(temp,ne,time,J_nu);
    else atoms[i].solve_lte (temp,ne,time);   

    // total electron donation from this atomic species
    f += dens*mass_frac[i]/(elem_A[i]*pc::m_p)*atoms[i].get_ion_frac();
  }

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
double nlte_gas::ne_brent_method(double x1,double x2,double tol,std::vector<real> J_nu)
{  
  int ITMAX = 100;
  double EPS = 3.0e-8;
  int iter;
  double a=x1,b=x2,c=x2,d,e,min1,min2;
  double fa=charge_conservation(a,J_nu);
  double fb=charge_conservation(b,J_nu);
  double fc,p,q,r,s,tol1,xm;
  
  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
    solve_error_ = 1;
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
    fb=charge_conservation(b,J_nu);
  }
  solve_error_ = 2;
  return 0.0;
}



//-----------------------------------------------------------
// return the fraction of atoms of index i that are in
// ionization state j.  
//-----------------------------------------------------------
double nlte_gas::get_ionization_fraction(int i, int j)
{
  if ((i < 0)||(i >= atoms.size()))    return -1;
  if ((j < 0)||(j >= atoms[i].n_ions)) return -1;
  return atoms[i].ionization_fraction(j);
}



//-----------------------------------------------------------
// return the fraction of atoms of index i that are in
// ionization state j.  
//-----------------------------------------------------------
double nlte_gas::get_level_fraction(int i, int j)
{
  if ((i < 0)||(i >= atoms.size()))    return -1;
  if ((j < 0)||(j >= atoms[i].n_levels)) return -1;
  return atoms[i].level_fraction(j);
}

//-----------------------------------------------------------
// return the fraction of atoms of index i that are in
// ionization state j.  
//-----------------------------------------------------------
double nlte_gas::get_level_departure(int i, int j)
{
  if ((i < 0)||(i >= atoms.size()))    return -1;
  if ((j < 0)||(j >= atoms[i].n_levels)) return -1;
  return atoms[i].level_depature(j);
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
 
