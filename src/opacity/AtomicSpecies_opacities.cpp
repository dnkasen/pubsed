#include "AtomicSpecies.h"
#include "physical_constants.h"
#include <iostream>
#include <limits>

namespace pc = physical_constants;

//---------------------------------------------------------
// calculate the bound-free extinction coefficient
// (units cm^{-1}) and emissivity from a sum over all
// atomic levels.  Uses current gas_temp_
// Results are put into opac and emis vectors
// passed the free electron density ne
//---------------------------------------------------------
void AtomicSpecies::bound_free_opacity
(std::vector<double>& opac, std::vector<double>& emis, double ne)
{
  bound_free_opacity_general(opac,emis,ne,gas_temp_,0);
}

//---------------------------------------------------------
// calculate the bound-free emissivity for the purposees
// of calculating cooling.
// Results are put into emis vector
// passed the free electron density ne and temperature T
//---------------------------------------------------------
void AtomicSpecies::bound_free_opacity_for_cooling
(std::vector<double>& emis, double ne, double T)
{
  // a place holder
  std::vector<double> opac;
  bound_free_opacity_general(opac,emis,ne,T,1);
}

//---------------------------------------------------------
// calculate the bound-free opacity for the purposees
// of calculating heating.
// Results are put into opac vector
// passed the free electron density ne and temperature T
//---------------------------------------------------------
void AtomicSpecies::bound_free_opacity_for_heating
(std::vector<double>& opac, double ne, double T)
{
  // a place holder
  std::vector<double> emis;
  bound_free_opacity_general(opac,emis,ne,T,2);
}


//---------------------------------------------------------
// calculate the bound-free extinction coefficient
// and emissivity from a sum over all atomic levels
// The coolheat flag provides different options
//   coolheat = 0  (calculate straight ahead opacity emissivity)
//   coolheat = 1  (calculate emissivity for cooling)
//   coolheat = 2  (calculate opacity for heating)
//---------------------------------------------------------
void AtomicSpecies::bound_free_opacity_general
(std::vector<double>& opac, std::vector<double>& emis, double ne, double T, int coolheat)
{
  // zero out arrays
  for (size_t i=0;i<opac.size();++i)
  {
    if ((coolheat == 0)||(coolheat == 2))
      opac[i] = 0;
    if ((coolheat == 0)||(coolheat == 1))
      emis[i] = 0;
  }

  int ng = nu_grid_.size();
  double kt_ev = pc::k_ev*T;
  double lam_t   = sqrt(pc::h*pc::h/(2*pc::pi*pc::m_e* pc::k * T));

  std::vector<double> nc_phifac(n_levels_);
  for (int j=0;j<n_levels_;++j)
  {
    int ic = adata_->get_lev_ic(j);
    if (ic == -1) continue;
    double nc = n_dens_*lev_n_[ic];
    int gl = adata_->get_lev_g(j);
    int gc = adata_->get_lev_g(ic);
    double gl_o_gc = (1.0*gl)/(1.0*gc);
    nc_phifac[j] = nc*gl_o_gc/2. * lam_t * lam_t * lam_t;
  }

  // loop over and set opac/emis for each frequency
  for (int i=0;i<ng;++i)
  {
    double nu    = nu_grid_.center(i);
    double E     = pc::h*nu*pc::ergs_to_ev;
    double emis_fac   = 2. * pc::h*nu*nu*nu / pc::c / pc::c;

    // summing the contribution of every level
    for (int j=0;j<n_levels_;++j)
    {
      // check if above threshold and ionization state above
      double Eion = adata_->get_lev_Eion(j);
      int ic = adata_->get_lev_ic(j);
      if (ic == -1) continue;
      if (E < Eion) continue;

      // get extinction coefficient and emissivity
      double zeta_net = (Eion - E)/kt_ev;
      double ezeta_net = exp(zeta_net);
      double sigma = adata_->get_lev_photo_cs(j,E);
      double opac_fac = n_dens_ * lev_n_[j]  - nc_phifac[j] * ne * ezeta_net;
      // kill maser
      if (opac_fac < 0)
      	opac_fac = 0.;

      // store opacity (with extra factor if wanted for heating)
      if (coolheat == 0)
        opac[i] += sigma * opac_fac;
      else if (coolheat == 2)
        opac[i] += sigma * (opac_fac) * (E - Eion) * pc::ev_to_ergs;

      // store emissivity (don't add in ground lev if flag set)
      if ((adata_->get_lev_E(j) == 0)&&(no_ground_recomb_))
        continue;

      if (coolheat == 0)
	      emis[i]  += emis_fac *sigma* nc_phifac[j] * ezeta_net;
      else if (coolheat == 1)
        emis[i]  += emis_fac *sigma* nc_phifac[j] * ezeta_net * (E - Eion)/E;
    }

  }
}

//---------------------------------------------------------
// calculate the net cooling rate for collisional
// processes (bound-bound and bound-free)
// given the passed temperature T and free electron density ne
//---------------------------------------------------------
double AtomicSpecies::collisional_net_cooling_rate(double ne, double T)
{
  double collisional_net_cooling = 0.;

  //  bound-bound collisional transitions
  for (int l=0;l<n_lines_;l++)
  {
    int lu  = adata_->get_line_u(l);
    int ll  = adata_->get_line_l(l);
    double El = adata_->get_lev_E(ll);
    double Eu = adata_->get_lev_E(lu);
    double gl = 1.0*adata_->get_lev_g(ll);
    double gu = 1.0*adata_->get_lev_g(lu);

    double ndown   = n_dens_ * lev_n_[ll];
    double nup     = n_dens_ * lev_n_[lu];

    double dE = (Eu - El)*pc::ev_to_ergs;
    double zeta = dE/pc::k/T;
    double ezeta = exp(zeta);

    // floor oscillator strengths so that forbidden lines can contribute.
    // should be improved by using real collisional rates for forbidden lines
    double f_lu = adata_->get_line_f(l);
    double effective_f_lu = f_lu;
    if (f_lu < 1.e-3) effective_f_lu = 1.e-3;

    double C_up = 3.9*pow(zeta,-1.)*pow(T,-1.5) / ezeta * ne * effective_f_lu;
    if (zeta > 700) C_up = 0.; // be careful about overflow
    double C_down = 3.9*pow(zeta,-1.)*pow(T,-1.5) * ne * effective_f_lu * gl/gu;

    collisional_net_cooling += dE * (ndown * C_up - nup * C_down);
  }

  //bound-free collisional transitions:
  for (int i=0;i<n_levels_;++i)
  {
	  int ic = adata_->get_lev_ic(i);
	  if (ic == -1) continue;

	  // ionization potential
	  double chi  = adata_->get_lev_Eion(i)* pc::ev_to_ergs;
	  double zeta = chi/pc::k/T; // note chi is now in ergs

	  double nc = n_dens_ * lev_n_[ic];
	  double ni = n_dens_ * lev_n_[i];

	  // collisional ionization rate
	  // needs to be multiplied by number of electrons in outer shell
	  double C_ion = 2.7/zeta/zeta*pow(T,-1.5)*exp(-zeta)*ne;

	  // collisional recombination rate
	  int gi = adata_->get_lev_g(i);
	  int gc = adata_->get_lev_g(ic);
	  double C_rec = 5.59080e-16/zeta/zeta*pow(T,-3)*gi/gc*ne*ne;

	  collisional_net_cooling += chi * ( ni * C_ion - nc * C_rec);
  }

  return collisional_net_cooling;
}

//---------------------------------------------------------
// calculate the bound-bound (i.e., line) extinction coefficient
// (units cm^{-1}) and emissivity for all lines
//---------------------------------------------------------
void AtomicSpecies::bound_bound_opacity(std::vector<double>& opac, std::vector<double>& emis)
{
  // zero out arrays
  for (size_t i=0;i<opac.size();++i) {opac[i] = 0; emis[i] = 0;}

  // loop over all lines
  for (int i=0;i<n_lines_;++i)
  {
    int lu  = adata_->get_line_u(i);
    int ll  = adata_->get_line_l(i);

    double nl    = lev_n_[ll];
    double nu    = lev_n_[lu];
    double gl    = 1.0*adata_->get_lev_g(ll);
    double gu    = 1.0*adata_->get_lev_g(lu);
    double nu_0  = adata_->get_line_nu(i);

    double dnu   = line_beta_dop_*nu_0;
    double A_ul  = adata_->get_line_A(i);
    double gamma = A_ul;
    double a_voigt = gamma/4/pc::pi/dnu;

    // extinction coefficient
    if (nl == 0) continue;
    double alpha_0 = nl*n_dens_*gu/gl*A_ul/(8*pc::pi)*pc::c*pc::c;
    // correction for stimulated emission
    alpha_0 = alpha_0*(1 - nu*gl/(nl*gu));

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
    double line_j = A_ul*nu*n_dens_*pc::h/(4.0*pc::pi);
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


//---------------------------------------------------------
// Calculate the extinction coefficient (units cm^{-1})
// for lines in the Sobolev expansion opacity formalism
// using the fuzz lines
//---------------------------------------------------------
void AtomicSpecies::fuzzline_expansion_opacity
(std::vector<double>& opac, double time)
{
  // zero out opacity array
  std::fill(opac.begin(),opac.end(),0);

  // loop over all lines
  int n_flines = adata_->get_n_fuzz_lines();
  for (int i=0;i<n_flines;++i)
  {
    // properties of this line
    double nu = adata_->get_fuzz_line_nu(i);
    double gf = adata_->get_fuzz_line_gf(i);
    double El = adata_->get_fuzz_line_El(i);
    int   ion = adata_->get_fuzz_line_ion(i);
    int   bin = adata_->get_fuzz_line_bin(i);

    // calculate Sobole optical depth
    double nion = n_dens_*ion_frac_[ion];
    double nl = nion*exp(-1.0*El/pc::k_ev/gas_temp_)/ion_part_[ion];
    double lam = pc::c/nu;
    double stim_cor = (1 - exp(-pc::h*nu/pc::k/gas_temp_));
    double tau = pc::sigma_tot*lam*nl*gf*stim_cor*time;

    // bin the lines
    double etau = exp(-tau);
    opac[bin] += (1 - etau);
  }

  // renormalize opacity array
  for (size_t i=0;i<opac.size();i++)
    opac[i] = opac[i]*nu_grid_.center(i)/nu_grid_.delta(i)/pc::c/time;
}


//---------------------------------------------------------
// Calculate the extinction coefficient (units cm^{-1})
// for lines in the Sobolev expansion opacity formalism
//---------------------------------------------------------
void AtomicSpecies::line_expansion_opacity
(std::vector<double>& opac, double time)
{
  // zero out opacity array
  std::fill(opac.begin(),opac.end(),0);

  // loop over all lines
  for (int i=0;i<n_lines_;++i)
  {
    int    ll  = adata_->get_line_l(i);
    int    lu  = adata_->get_line_u(i);
    double nl  = lev_n_[ll];
    double nu  = lev_n_[lu];
    double gl  = 1.0*adata_->get_lev_g(ll);
    double gu  = 1.0*adata_->get_lev_g(lu);

    // don't bother with unpopulated levels
    if (nl < std::numeric_limits<double>::min())
      continue;

    // Calculate Sobolev optical depth tau
    double nu_0  = adata_->get_line_nu(i);
    double lam   = pc::c/nu_0;
    double f_lu  = adata_->get_line_f(i);
    double tau   = nl*n_dens_*pc::sigma_tot*f_lu*time*lam;

    // correction for stimulated emission
    tau = tau*(1 - nu*gl/(nl*gu));
    if (nu*gl > nl*gu)
    {
       //   printf("laser regime, line %d, whoops\n",i);
       tau = 0;
    }

    // bin the lines
    double etau = exp(-tau);
    int    bin   = adata_->get_line_bin(i);
    if (!isnan(etau))
      opac[bin] += (1 - etau);
  }

  // renormalize opacity array
  for (size_t i=0;i<opac.size();i++)
     opac[i] = opac[i]*nu_grid_.center(i)/nu_grid_.delta(i)/pc::c/time;
}
