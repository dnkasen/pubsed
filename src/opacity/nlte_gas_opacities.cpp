#include "nlte_gas.h"
#include "physical_constants.h"
#include <iostream>

namespace pc = physical_constants;


//----------------------------------------------------------------
// calculate the total absorptive and scattering opacity
//----------------------------------------------------------------
void nlte_gas::computeOpacity(std::vector<double>& abs, 
			      std::vector<double>& scat, 
			      std::vector<double>& emis)
{
  int ns = nu_grid.size();
  std::vector<double> opac, evec;
  opac.resize(ns);
  evec.resize(ns);

  // zero out passed opacity arrays
  for (int i=0;i<ns;i++) {abs[i] = 0; scat[i] = 0;}
  
  //-----------------------------------------
  /// if grey opacity, just do simple thing
  //-----------------------------------------
  if (grey_opacity_ != 0)
  {
    double opac = dens*grey_opacity_;
    for (int i=0;i<ns;i++) {
      abs[i]  = opac*epsilon_;
      scat[i] = opac*(1-epsilon_);
    }
  } 

  //-----------------------------------------
  // else do all the other stuff
  //-----------------------------------------
  else
  {
    //---
    if (use_electron_scattering_opacity) 
    {
      double es_opac = electron_scattering_opacity();
      for (int i=0;i<ns;i++) {
        scat[i] += es_opac;
	       // debug -- a bit of small thermalization
	       abs[i] += es_opac*1e-10; }
    }

     //---
    if (use_free_free_opacity) 
    {
      free_free_opacity(opac);
      for (int i=0;i<ns;i++) {
         abs[i] += opac[i]; }
    }

     //---
    if (use_bound_free_opacity) 
    {
      bound_free_opacity(opac);
      for (int i=0;i<ns;i++) {
         abs[i] += opac[i]; }
    }

    //---
    if (use_bound_bound_opacity) 
    {
      bound_bound_opacity(opac);
      for (int i=0;i<ns;i++) {
         abs[i] += opac[i]; }
    }



    //---
    if (use_line_expansion_opacity) 
    {
      line_expansion_opacity(opac);
      for (int i=0;i<ns;i++) {
    	 abs[i]  += epsilon_*opac[i];
    	 scat[i] += (1-epsilon_)*opac[i];
        }
    }
  
    //---
    if (use_fuzz_expansion_opacity) 
    {
      fuzz_expansion_opacity(opac);
      for (int i=0;i<ns;i++) {
	     abs[i]  += epsilon_*opac[i];
	     scat[i] += (1-epsilon_)*opac[i];
      }
    }
  }

  
  //-- set emissivity, just thermal now
  for (int i=0;i<ns;i++)
  {
    double nu = nu_grid.center(i);
    double dnu = nu_grid.delta(i);

    double zeta = pc::h*nu/pc::k/temp;
    double bb = 2.0*nu*nu*nu*pc::h/pc::c/pc::c/(exp(zeta)-1);
    emis[i] = bb*abs[i]*dnu;
  }
}


//----------------------------------------------------------------
// simple electron scattering opacity
//----------------------------------------------------------------
double nlte_gas::electron_scattering_opacity()
{
  return pc::thomson_cs*ne;
}


//----------------------------------------------------------------
// free-free opacity (brehmstrahlung)
// note the gaunt factor is being set to 1 here
//----------------------------------------------------------------
void nlte_gas::free_free_opacity(std::vector<double>& opac)
{
  int npts   = nu_grid.size();
  int natoms = atoms.size();

  // calculate sum of n_ion*Z**2
  double fac = 0;
  for (int i=0;i<natoms;i++)
  {
    double Z_eff_sq = 0;
    for (int j=0;j<atoms[i].n_ions;j++) 
      Z_eff_sq += atoms[i].ionization_fraction(j)*j*j;

    double n_ion = mass_frac[i]*dens/(elem_A[i]*pc::m_p);
    fac += n_ion*Z_eff_sq;
  }
  // multiply by overall constants
  fac *= 3.7e8*pow(temp,-0.5)*ne;
  
  // multiply by frequency dependence
  for (int i=0;i<npts;i++)
  {
    double nu = nu_grid[i];
    opac[i] = fac/nu/nu/nu*(1 - exp(-pc::h*nu/pc::k/temp));
  }      
  
}


//----------------------------------------------------------------
// bound-free opacity (photoionization) from all elements
//----------------------------------------------------------------
void nlte_gas::bound_free_opacity(std::vector<double>& opac)
{
  // zero out opacity vector
  int ng = nu_grid.size();
  for (int j=0;j<ng;j++) opac[j] = 0;

  std::vector<double> atom_opac(ng);
  int na = atoms.size();

  // sum up the bound-free opacity from every atom
  for (int i=0;i<na;i++)
  {
    atoms[i].bound_free_opacity(atom_opac);
    for (int j=0;j<ng;j++) opac[j] += atom_opac[j];
  }
}



//----------------------------------------------------------------
// bound-bound opacity (lines)
// uses Doppler broadening for now
//----------------------------------------------------------------
void nlte_gas::bound_bound_opacity(std::vector<double>& opac)
{
  // zero out opacity vector
  int ng = nu_grid.size();
  for (int j=0;j<ng;j++) opac[j] = 0;
  
  std::vector<double> atom_opac(ng);
  int na = atoms.size();
  
  // sum up the bound-bound opacity from every atom
  for (int i=0;i<na;i++)
  {
    bound_bound_opacity(i,atom_opac);
    for (int j=0;j<ng;j++) opac[j] += atom_opac[j];
  }
}


//----------------------------------------------------------------
// bound-bound opacity (line) from single element iatom
//----------------------------------------------------------------
void nlte_gas::bound_bound_opacity(int iatom, std::vector<double>& opac)
{
  atoms[iatom].bound_bound_opacity(opac);
}





//----------------------------------------------------------------
// Calculate a binned expansion opacity based on nlte line data
// Passed:
//   opac -- double vector of the same size of the frequency
//   grid, which will be filled up with the opacities
// UNITS are cm^{-1} 
// So this is really an extinction coefficient
//----------------------------------------------------------------
void nlte_gas::line_expansion_opacity(std::vector<double>& opac)
{
  // zero out opacity array
  std::fill(opac.begin(),opac.end(),0);

  // loop over atoms
  for (int i=0;i<atoms.size();i++)
  {
    //compute line sobolev taus
    atoms[i].compute_sobolev_taus(time);

    // loop over all lines
    for (int j=0;j<atoms[i].n_lines;j++)
    {
      // add in this line to the sum
      double etau = atoms[i].lines[j].etau;
      opac[atoms[i].lines[j].bin] += (1 - etau);
    }
  }

  // renormalize opacity array
  for (int i=0;i<opac.size();i++) 
    opac[i] = opac[i]*nu_grid.center(i)/nu_grid.delta(i)/pc::c/time;
}

//----------------------------------------------------------------
// Calculate a binned expansion opacity based on fuzz line data
// Passed:
//   opac -- double vector of the same size of the frequency
//   grid, which will be filled up with the opacities
// UNITS are cm^{-1} 
// So this is really an extinction coefficient
//----------------------------------------------------------------
void nlte_gas::fuzz_expansion_opacity(std::vector<double>& opac)
{
  double exp_min = 1e-6;
  double exp_max = 100;

  // zero out opacity array
  for (int i=0;i<opac.size();i++) opac[i] = 0;

  // loop over atoms
  for (int i=0;i<atoms.size();i++)
  {
    // loop over lines in atom
    for (int j=0;j<atoms[i].fuzz_lines.n_lines;j++)
    {
      int   ion = atoms[i].fuzz_lines.ion[j];
      double gf = atoms[i].fuzz_lines.gf[j];
      double El = atoms[i].fuzz_lines.El[j];
      double nu = atoms[i].fuzz_lines.nu[j];

      // get sobolev tau
      double n_dens = mass_frac[i]*dens/(elem_A[i]*pc::m_p);
      double nion   = n_dens*atoms[i].ionization_fraction(ion);

      double nl = nion*exp(-1.0*El/pc::k_ev/temp)/atoms[i].partition(ion);
      double lam = pc::c/nu;
      double stim_cor = (1 - exp(-pc::h*nu/pc::k/temp));
      double tau = pc::sigma_tot*lam*nl*gf*stim_cor*time;
      
      // effeciently calculate exponential of tau
      double etau;
      if (tau < exp_min)      etau = 1 - tau;
      else if (tau > exp_max) etau = 0;
      else                    etau = exp(-tau);

      // add in this line to the sum
      int bin = atoms[i].fuzz_lines.bin[j];
      opac[bin] += (1 - etau);
    }
  }

  // renormalize opacity array
  for (int i=0;i<opac.size();i++) 
    opac[i] = opac[i]*nu_grid.center(i)/nu_grid.delta(i)/pc::c/time;
}

//---------------



//----------------------------------------------------------------
// Calculate a binned expansion opacity based on nlte line data
// You must multiply this by the line profile phi (units Hz^-1)
// in order to get the exctinction coefficient
// Passed:
//   opac -- double vector of the same size of the frequency
//   grid, which will be filled up with the opacities
// UNITS are cm^{-1} Hz 
// So this is really an extinction coefficient
//----------------------------------------------------------------
void nlte_gas::get_line_opacities(std::vector<double>& opac)
{
  int i = 0;

  std::vector<nlteGlobalLine>::iterator iLine = globalLineList_.begin();
  while (iLine != globalLineList_.end())
  {
    int iAtom = iLine->iAtom;
    int iLo   = iLine->iLowerLevel;
    int iHi   = iLine->iUpperLevel;
    double gr = iLine->g1_over_g2;
    
    double nd = atoms[iAtom].n_dens;
    double n1 = nd*atoms[iAtom].levels[iLo].n;
    double n2 = nd*atoms[iAtom].levels[iHi].n;
  
    double line_cs = pc::sigma_tot*iLine->f_osc;
    opac[i] = line_cs*n1*(1 - gr*n2/n1);
    if (n1 == 0) opac[i] = 0.0;

    i++;
    iLine++;
  }

}
