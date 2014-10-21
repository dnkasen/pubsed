#include "nlte_gas.h"
#include "physical_constants.h"
#include <iostream>

namespace pc = physical_constants;

double nlte_gas::electron_scattering_opacity()
{
  return pc::thomson_cs*ne;
}


void nlte_gas::line_expansion_opacity(std::vector<double>& opac)
{
  // zero out opacity array
  std::fill(opac.begin(),opac.end(),0);

  // loop over atoms
  for (int i=0;i<atoms.size();i++)
  {
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
    std::vector<fuzz_line>::iterator f;
    for (f = atoms[i].fuzz_lines.begin();f != atoms[i].fuzz_lines.end();f++)
    {
      // get sobolev tau
      double n_dens = mass_frac[i]*dens/(elem_A[i]*pc::m_p);
      double nion   = n_dens*atoms[i].ion_frac(f->ion);

      double nl = nion*exp(-1.0*f->El/pc::k_ev/temp)/atoms[i].partition(f->ion);
      double lam = pc::c/f->nu;
      double stim_cor = (1 - exp(-pc::h*f->nu/pc::k/temp));
      double tau = pc::sigma_tot*lam*nl*f->gf*stim_cor*time;

      // effeciently calculate exponential of tau
      double etau;
      if (tau < exp_min)      etau = 1 - tau;
      else if (tau > exp_max) etau = 0;
      else                    etau = exp(-tau);

      // add in this line to the sum
      opac[f->bin] += (1 - etau);
      //printf("%e %e %d %e\n",lam*1e8,f->El,f->bin,tau);
    }
  }

  // renormalize opacity array
  for (int i=0;i<opac.size();i++) 
    opac[i] = opac[i]*nu_grid.center(i)/nu_grid.delta(i)/pc::c/time;
}
