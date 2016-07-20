#include "nlte_atom.h"
#include "physical_constants.h"
#include <iostream>

namespace pc = physical_constants;



//---------------------------------------------------------
// calculate the bound-free extinction coefficient
// (units cm^{-1}) for all levels
//---------------------------------------------------------
void nlte_atom::bound_free_opacity(std::vector<double>& opac, std::vector<double>& emis)
{
  // zero out arrays
  for (int i=0;i<opac.size();i++) {opac[i] = 0; emis[i] = 0;}

  int ng = nu_grid.size();
  double lam_t   = sqrt(pc::h*pc::h/(2*pc::pi*pc::m_e*pc::k*gas_temp_));
  double lam_fac = 2*pc::h/pc::c/pc::c*lam_t*lam_t*lam_t/2.0;
  double kt = pc::k_ev*gas_temp_;

  for (int i=0;i<ng;i++)
  {
    opac[i] = 0;
    double nu    = nu_grid[i];
    double E     = pc::h*nu*pc::ergs_to_ev;
    double fac   = lam_fac*nu*nu*nu;

    for (int j=0;j<n_levels;j++)
    {
      // check if their is an ionization stage above
      int ic = levels[j].ic;
      if (ic == -1) continue;
      // check if above threshold
      double E_t = levels[j].E_ion;
      if (E < E_t) continue;
      // get extinction coefficient and emissivity
      double sigma = levels[j].s_photo.value_at_with_zero_edges(E);
      opac[i]  += n_dens*sigma*levels[j].n;
      double nc = n_dens*levels[ic].n;
      double gl_o_gc = (1.0*levels[j].g)/(1.0*levels[ic].g);

      // correction for stim emis  debug
      double se = (1 + 1.0/(exp(E/kt) - 1));
      emis[i]  += nc*sigma*gl_o_gc*fac*exp((E_t - E)/kt)*se;
    }
  }
}

//---------------------------------------------------------
// calculate the bound-free extinction coefficient
// (units cm^{-1}) for all levels
//---------------------------------------------------------
void nlte_atom::bound_bound_opacity(std::vector<double>& opac, std::vector<double>& emis)
{
  // zero out arrays
  for (int i=0;i<opac.size();i++) {opac[i] = 0; emis[i] = 0;}

  // loop over all lines
  for (int i=0;i<n_lines;i++)
  {
    int ll = lines[i].ll;
    int lu = lines[i].lu;

    double nlow  = levels[ll].n;
    double nup   = levels[lu].n;
    double glow  = levels[ll].g;
    double gup   = levels[lu].g;
    double nu_0  = lines[i].nu;

    double dnu = line_beta_dop_*nu_0;
    double gamma = lines[i].A_ul;
    double a_voigt = gamma/4/pc::pi/dnu;

    // extinction coefficient
    double alpha_0 = nlow*n_dens*gup/glow*lines[i].A_ul/(8*pc::pi)*pc::c*pc::c;
    // correction for stimulated emission
    alpha_0 = alpha_0*(1 - nup*glow/(nlow*gup));

    //if (alpha_0 < 0) std::cout << "LASER " << levels[ll].E << " " << levels[lu].E << "\n";
    //if (alpha_0 < 0) {std::cout << "LASER: " << nlow*gup << " " << nup*glow << "\n"; continue;}
    if (alpha_0 < 0) continue; 
    
    // region to add to -- hard code to 20 doppler widths
    double nu_1 = nu_0 - dnu*100;
    double nu_2 = nu_0 + dnu*100;
    int inu1 = nu_grid.locate(nu_1);
    int inu2 = nu_grid.locate(nu_2);

    // line emissivity: ergs/sec/cm^3/str
    // multiplied by phi below to get per Hz
    double line_j = lines[i].A_ul*nup*n_dens*pc::h/(4.0*pc::pi);
    for (int j = inu1;j<inu2;j++)
    {
      double nu = nu_grid[j];
      double x  = (nu_0 - nu)/dnu;
      double phi = voigt_profile_.getProfile(x,a_voigt)/dnu;
      opac[j] += alpha_0/nu/nu*phi;
      emis[j] += line_j*nu*phi;
    //  std::cout << nu_0 << " " << nu-nu_0 << " " << line_j*nu*phi/(alpha_0/nu/nu*phi) << " " 
    //    << 2*pc::h*nu*nu*nu/pc::c/pc::c/(exp(pc::h*nu/pc::k/gas_temp_) - 1) << "\n"; 
    }
    //printf("%e %e %e %d %d\n",lines[i].nu,nl,alpha_0,inu1,inu2,nu_grid[inu2]);
  } 
  
}




void nlte_atom::compute_sobolev_taus(double time)
{
  for (int i=0;i<n_lines;i++) compute_sobolev_tau(i,time);
}

double nlte_atom::compute_sobolev_tau(int i, double time)
{
  int ll = lines[i].ll;
  int lu = lines[i].lu;

  double nl = levels[ll].n;
  double nu = levels[lu].n;
  double gl = levels[ll].g;
  double gu = levels[lu].g;

  // check for empty levels
  if (nl < std::numeric_limits<double>::min())
  { 
    lines[i].tau  = 0;
    lines[i].etau = 1;
    lines[i].beta = 1;
    return 0;
  }

  double lam   = pc::c/lines[i].nu;
  double tau   = nl*n_dens*pc::sigma_tot*lines[i].f_lu*time*lam;
  // correction for stimulated emission
  tau = tau*(1 - nu*gl/(nl*gu));

  if (nu*gl > nl*gu) {
//    printf("laser regime, line %d, whoops\n",i);
    lines[i].tau  = 0;
    lines[i].etau = 1;
    lines[i].beta = 1;
    return 0; }

  double etau = exp(-tau);
  lines[i].etau = etau;
  lines[i].tau  = tau;
  lines[i].beta = (1-etau)/tau;
  return lines[i].tau;
}
