#include "AtomicSpecies.h"
#include "physical_constants.h"
#include <iostream>
#include <limits>

namespace pc = physical_constants;



//---------------------------------------------------------
// calculate the bound-free extinction coefficient
// (units cm^{-1}) for all levels
//---------------------------------------------------------
void AtomicSpecies::bound_free_opacity(std::vector<double>& opac, std::vector<double>& emis, double ne)
{
  // zero out arrays
  for (size_t i=0;i<opac.size();++i) {opac[i] = 0; emis[i] = 0;}

  int ng = nu_grid_.size();
  double kt_ev = pc::k_ev*gas_temp_;
  double lam_t   = sqrt(pc::h*pc::h/(2*pc::pi*pc::m_e* pc::k * gas_temp_));

  std::vector<double> nc_phifac(n_levels_);
  for (int j=0;j<n_levels_;++j)
  {
    int ic = levels_[j].ic;
    if (ic == -1) continue;
    double nc = n_dens_*levels_[ic].n;
    double gl_o_gc = (1.0*levels_[j].g)/(1.0*levels_[ic].g);
    double E_t = levels_[j].E_ion;
    nc_phifac[j] = nc*gl_o_gc/2. * lam_t * lam_t * lam_t;
  }

  for (int i=0;i<ng;++i)
  {
    opac[i] = 0;
    double nu    = nu_grid_.center(i);
    double E     = pc::h*nu*pc::ergs_to_ev;
    double emis_fac   = 2. * pc::h*nu*nu*nu / pc::c / pc::c;

    for (int j=0;j<n_levels_;++j)
    {
      // check if above threshold
      if (E < levels_[j].E_ion) continue;
      // check if there is an ionization stage above
      int ic = levels_[j].ic;
      if (ic == -1) continue;

      // get extinction coefficient and emissivity
      double zeta_net = (levels_[j].E_ion - E)/kt_ev;
      double ezeta_net = exp(zeta_net);
      double sigma = levels_[j].s_photo.value_at_with_zero_edges(E);
      double opac_fac = n_dens_ * levels_[j].n  - nc_phifac[j] * ne * ezeta_net;
      // kill maser
      if (opac_fac < 0)
      	opac_fac = 0.;
      opac[i]  += sigma * opac_fac;

      if (levels_[j].E == 0)
	{	
	if (no_ground_recomb_ == 0)
	  {
	    emis[i]  += emis_fac *sigma* nc_phifac[j] * ezeta_net; // ne gets multiplied at the end outside this funciton
	  }
	}
      else
	{
	  emis[i]  += emis_fac *sigma* nc_phifac[j] * ezeta_net; // ne gets multiplied at the end outside this funciton
	}
    }

  }
}

void AtomicSpecies::bound_free_opacity_for_cooling(std::vector<double>& emis, double ne, double T)
{
  // zero out arrays
  for (size_t i=0;i<emis.size();++i) {emis[i] = 0;}

  int ng = nu_grid_.size();
  double kt_ev = pc::k_ev*T;
  double lam_t   = sqrt(pc::h*pc::h/(2*pc::pi*pc::m_e* pc::k * T));

  std::vector<double> nc_phifac(n_levels_);
  for (int j=0;j<n_levels_;++j)
  {
    int ic = levels_[j].ic;
    if (ic == -1) continue;
    double nc = n_dens_*levels_[ic].n;
    double gl_o_gc = (1.0*levels_[j].g)/(1.0*levels_[ic].g);
    nc_phifac[j] = nc*gl_o_gc/2. * lam_t * lam_t * lam_t;
  }

  for (int i=0;i<ng;++i)
  {
    double nu    = nu_grid_.center(i);
    double E     = pc::h*nu*pc::ergs_to_ev;
    double emis_fac   = 2. * pc::h*nu*nu*nu / pc::c / pc::c;

    for (int j=0;j<n_levels_;++j)
    {
      // check if above threshold
      if (E < levels_[j].E_ion) continue;
      // check if there is an ionization stage above
      int ic = levels_[j].ic;
      if (ic == -1) continue;

      // get extinction coefficient and emissivity
      double zeta_net = (levels_[j].E_ion - E)/kt_ev;
      double ezeta_net = exp(zeta_net);
      double sigma = levels_[j].s_photo.value_at_with_zero_edges(E);

      if (levels_[j].E == 0)
	{
	  if (no_ground_recomb_ == 0)
	    {
	      emis[i]  += emis_fac *sigma* nc_phifac[j] * ezeta_net * (E - levels_[j].E_ion)/E; // ne gets multiplied at the end outside this funciton
	    }
	}
      else
	{
	  emis[i]  += emis_fac *sigma* nc_phifac[j] * ezeta_net * (E - levels_[j].E_ion)/E; // ne gets multiplied at the end outside this funciton
	}


    }

  }
}

void AtomicSpecies::bound_free_opacity_for_heating(std::vector<double>& heat_opac, double ne, double T)
{

    // zero out arrays
  for (size_t i=0;i<heat_opac.size();++i) {heat_opac[i] = 0.;}
  
  int ng = nu_grid_.size();
  double kt_ev = pc::k_ev*T;
  double lam_t   = sqrt(pc::h*pc::h/(2*pc::pi*pc::m_e* pc::k * T));

  std::vector<double> nc_phifac(n_levels_);
  for (int j=0;j<n_levels_;++j)
  {
    int ic = levels_[j].ic;
    if (ic == -1) continue;
    double nc = n_dens_*levels_[ic].n;
    double gl_o_gc = (1.0*levels_[j].g)/(1.0*levels_[ic].g);
    nc_phifac[j] = nc*gl_o_gc/2. * lam_t * lam_t * lam_t;
  }

  for (int i=0;i<ng;++i)
  {
    double nu    = nu_grid_.center(i);
    double E     = pc::h*nu*pc::ergs_to_ev;

    for (int j=0;j<n_levels_;++j)
    {
      // check if above threshold
      if (E < levels_[j].E_ion) continue;
      // check if there is an ionization stage above
      int ic = levels_[j].ic;
      if (ic == -1) continue;

      // get extinction coefficient and emissivity
      double zeta_net = (levels_[j].E_ion - E)/kt_ev;
      double ezeta_net = exp(zeta_net);
      double sigma = levels_[j].s_photo.value_at_with_zero_edges(E);
      double opac_fac = n_dens_ * levels_[j].n  - nc_phifac[j] * ne * ezeta_net;
      if (opac_fac < 0)
	opac_fac = 0.;
      
      heat_opac[i]  += sigma * (opac_fac) * (E - levels_[j].E_ion) * pc::ev_to_ergs; // including this energy difference for each frequency bin is the whole point of this function
    }

  }
}

double AtomicSpecies::collisional_net_cooling_rate(double ne, double T)
{
  
    double collisional_net_cooling = 0.;

  //  bound-bound collisional transitions
  for (int l=0;l<n_lines_;l++)
    {
      int lu  = lines_[l].lu;
      int ll  = lines_[l].ll;

      double ndown   = n_dens_ * levels_[ll].n;
      double nup     = n_dens_ * levels_[lu].n;

      double dE = (levels_[lu].E - levels_[ll].E)*pc::ev_to_ergs;

      double zeta = dE/pc::k/T; // note dE is in ergs
      double ezeta = exp(zeta);

      double effective_f_lu = 0.;
      if (lines_[l].f_lu < 1.e-3) effective_f_lu = 1.e-3; // just flooring the oscillator strengths so that forbidden lines can contribute. Should be improved
      else effective_f_lu = lines_[l].f_lu;

      double C_up = 3.9*pow(zeta,-1.)*pow(T,-1.5) / ezeta * ne * effective_f_lu;
      if (zeta > 700) C_up = 0.; // be careful about overflow
      double C_down = 3.9*pow(zeta,-1.)*pow(T,-1.5) * ne * effective_f_lu * levels_[ll].g/levels_[lu].g;
      
      collisional_net_cooling += dE * (ndown * C_up - nup * C_down); 

    }

  //bound-free collisional transitions:

    for (int i=0;i<n_levels_;++i)
      {
	int ic = levels_[i].ic;
	if (ic == -1) continue;

	    // ionization potential
	int istage  = levels_[i].ion;
	double chi  = (ions_[istage].chi - levels_[i].E) * pc::ev_to_ergs;
	double zeta = chi/pc::k/T; // note chi is now in ergs

	double nc = n_dens_ * levels_[ic].n;
	double ni = n_dens_ * levels_[i].n;

	// collisional ionization rate
	// needs to be multiplied by number of electrons in outer shell
	double C_ion = 2.7/zeta/zeta*pow(T,-1.5)*exp(-zeta)*ne;

	// collisional recombination rate
	int gi = levels_[i].g;
	int gc = levels_[ic].g;
	double C_rec = 5.59080e-16/zeta/zeta*pow(T,-3)*gi/gc*ne*ne;

	collisional_net_cooling += chi * ( ni * C_ion - nc * C_rec);
      }




  return collisional_net_cooling;

  
}

//---------------------------------------------------------
// calculate the bound-free extinction coefficient
// (units cm^{-1}) for all levels
//---------------------------------------------------------
void AtomicSpecies::bound_bound_opacity(std::vector<double>& opac, std::vector<double>& emis)
{
  // zero out arrays
  for (size_t i=0;i<opac.size();++i) {opac[i] = 0; emis[i] = 0;}

  // loop over all lines
  for (int i=0;i<n_lines_;++i)
  {
    int ll = lines_[i].ll;
    int lu = lines_[i].lu;

    double nlow  = levels_[ll].n;
    double nup   = levels_[lu].n;
    double glow  = levels_[ll].g;
    double gup   = levels_[lu].g;
    double nu_0  = lines_[i].nu;

    double dnu = line_beta_dop_*nu_0;
    double gamma = lines_[i].A_ul;
    double a_voigt = gamma/4/pc::pi/dnu;

    // extinction coefficient
    if (nlow == 0) continue;
    double alpha_0 = nlow*n_dens_*gup/glow*lines_[i].A_ul/(8*pc::pi)*pc::c*pc::c;
    // correction for stimulated emission
    alpha_0 = alpha_0*(1 - nup*glow/(nlow*gup));

    //if (alpha_0 < 0) std::cout << "LASER " << levels_[ll].E << " " << levels[lu].E << "\n";
    //if (alpha_0 < 0) {std::cout << "LASER: " << nlow*gup << " " << nup*glow << "\n"; continue;}
    if (alpha_0 <= 0) continue;

    //if (alpha_0/nu_0/nu_0/dnu*1e15 < 1e-10) continue;

    // don't bother calculating very small opacities
    if (alpha_0/(nu_0*nu_0*dnu) < minimum_extinction_) continue;

    // region to add to -- hard code to 20 doppler widths
    double nu_1 = nu_0 - dnu*5; //*30;
    double nu_2 = nu_0 + dnu*5; //*30; //debug
    int inu1 = nu_grid_.locate_within_bounds(nu_1);
    int inu2 = nu_grid_.locate_within_bounds(nu_2);


    // line emissivity: ergs/sec/cm^3/str
    // multiplied by phi below to get per Hz
    double line_j = lines_[i].A_ul*nup*n_dens_*pc::h/(4.0*pc::pi);
    for (int j = inu1;j<inu2;++j)
    {
      double nu = nu_grid_.center(j);
      double x  = (nu_0 - nu)/dnu;
      double phi = voigt_profile_.getProfile(x,a_voigt)/dnu;
      opac[j] += alpha_0/nu/nu*phi;
      emis[j] += line_j*nu*phi;

      //if (isnan(alpha_0/nu/nu*phi)) std::cout << alpha_0 << " " << nlow << " " << nup << "\n";
    //  std::cout << nu_0 << " " << nu-nu_0 << " " << line_j*nu*phi/(alpha_0/nu/nu*phi) << " "
    //    << 2*pc::h*nu*nu*nu/pc::c/pc::c/(exp(pc::h*nu/pc::k/gas_temp_) - 1) << "\n";
    }
    //printf("%e %e %e %d %d\n",lines_[i].nu,nl,alpha_0,inu1,inu2,nu_grid[inu2]);
  }
}




void AtomicSpecies::compute_sobolev_taus(double time)
{
  for (int i=0;i<n_lines_;++i) compute_sobolev_tau(i,time);
}

double AtomicSpecies::compute_sobolev_tau(int i, double time)
{
  int ll = lines_[i].ll;
  int lu = lines_[i].lu;

  double nl = levels_[ll].n;
  double nu = levels_[lu].n;
  double gl = levels_[ll].g;
  double gu = levels_[lu].g;

  // check for empty levels_
  if (nl < std::numeric_limits<double>::min())
  {
    lines_[i].tau  = 0;
    lines_[i].etau = 1;
    lines_[i].beta = 1;
    return 0;
  }

  double lam   = pc::c/lines_[i].nu;
  double tau   = nl*n_dens_*pc::sigma_tot*lines_[i].f_lu*time*lam;
  // correction for stimulated emission
  tau = tau*(1 - nu*gl/(nl*gu));

  if (nu*gl > nl*gu) {
//    printf("laser regime, line %d, whoops\n",i);
    lines_[i].tau  = 0;
    lines_[i].etau = 1;
    lines_[i].beta = 1;
    return 0; }

  double etau = exp(-tau);
  lines_[i].etau = etau;
  lines_[i].tau  = tau;
  lines_[i].beta = (1-etau)/tau;
  return lines_[i].tau;
}
