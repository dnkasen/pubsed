#include <iostream>
#include <fstream>
#include <stdio.h>
#include "GasState.h"
#include "locate_array.h"
#include "physical_constants.h"
#include "ParameterReader.h"
#include "sedona.h"

using namespace std;
namespace pc = physical_constants;

static locate_array nu_grid;
static GasState gas;
static ParameterReader params;
static std::vector<OpacityType> abs_opacity, scat_opacity, tot_opacity, emissivity;

// functions
/*
static double rosseland_mean(std::vector<OpacityType> opac);
static double planck_mean(std::vector<OpacityType> opac);
static void write_mesa_file(std::string);
static void write_frequency_file(std::string, int);
static void write_gas_state(std::string);
static void write_mean_opacities(std::string);
static int verbose;
*/

int main(int argc, char **argv)
{
}

/*
  verbose = 1;

  // open up the parameter reader
  std::string param_file = "param.lua";
  if( argc > 1 ) param_file = std::string( argv[ 1 ] );
  params.initialize(param_file,verbose);

  // initialize the frequency grid
  std::vector<double> nu_dims = params.getVector<double>("transport_nu_grid");
  if ((nu_dims.size() != 4)&&(nu_dims.size() != 3)) {
    cerr << "# improperly defined nu_grid; need {nu_1, nu_2, dnu, (log?)}; exiting" << endl;
    exit(1); }
  if (nu_dims.size() == 3)
    nu_grid.init(nu_dims[0],nu_dims[1],nu_dims[2]);
  if (nu_dims.size() == 4)
  {
    if (nu_dims[3] == 1) nu_grid.log_init(nu_dims[0],nu_dims[1],nu_dims[2]);
    else nu_grid.init(nu_dims[0],nu_dims[1],nu_dims[2]);
  }

  //-------------------------------------------------------
  // get the elements to use
  //-------------------------------------------------------
  std::vector<int> elems_Z = params.getVector<int>("elements_Z");
  std::vector<int> elems_A = params.getVector<int>("elements_A");

  //-------------------------------------------------------
  // Set up the gas class
  //-------------------------------------------------------
  std::string atomdata = params.getScalar<string>("data_atomic_file");
  gas.initialize(atomdata,elems_Z,elems_A,nu_grid);
  std::string fuzzfile = params.getScalar<string>("data_fuzzline_file");
  int nl = gas.read_fuzzfile(fuzzfile);
  std::cout << "# From fuzzfile \"" << fuzzfile << "\" " <<
     nl << " lines used\n";gas.read_fuzzfile(fuzzfile);
  std::vector<double> massfrac = params.getVector<double>("mass_fractions");
  gas.set_mass_fractions(massfrac);
  gas.time = params.getScalar<double>("time");

  // opacity parameters
  gas.use_nlte_ = 0;
  gas.use_electron_scattering_opacity
    = params.getScalar<int>("opacity_electron_scattering");
  gas.use_line_expansion_opacity
    = params.getScalar<int>("opacity_line_expansion");
  gas.use_fuzz_expansion_opacity
    = params.getScalar<int>("opacity_fuzz_expansion");
  gas.use_bound_free_opacity
    = params.getScalar<int>("opacity_bound_free");
  gas.use_bound_bound_opacity
    = params.getScalar<int>("opacity_bound_bound");
  gas.use_free_free_opacity
    = params.getScalar<int>("opacity_free_free");
  gas.grey_opacity_ = 0;
  gas.line_velocity_width_ = params.getScalar<double>("line_velocity_width");

  // ------------------------------------------
  // allocate vectors for opacities/emissivities
  // ------------------------------------------
  int ng = nu_grid.size();
  abs_opacity.resize(ng);
  scat_opacity.resize(ng);
  tot_opacity.resize(ng);
  emissivity.resize(ng);

  gas.print_properties();

  // ------------------------------------------
  // Do output requested
  // ------------------------------------------
  std::string mesafile = params.getScalar<string>("output_mesa_file");
  if (mesafile != "") write_mesa_file(mesafile);

  std::string frequencyfile = params.getScalar<string>("output_frequency_file");
  if (frequencyfile != "") write_frequency_file(frequencyfile,0);

  std::string lambdafile = params.getScalar<string>("output_wavelength_file");
  if (lambdafile != "") write_frequency_file(lambdafile,1);

  std::string statefile = params.getScalar<string>("output_gas_state");
  if (statefile != "") write_gas_state(statefile);

  std::string meanfile = params.getScalar<string>("output_mean_opacities");
  if (meanfile != "") write_mean_opacities(meanfile);

}


//*********************************************************
// Write opacities as function of frequency
//*********************************************************
void write_frequency_file(std::string outfile, int style)
{
  FILE *fout = fopen(outfile.c_str(),"w");
  if (fout == NULL)
  {
    if (verbose) std::cerr << "Can't open opacity output file " << outfile << std::endl;
    return;
  }

  vector<double> temperature_list = params.getArray<double>("temperature");
  vector<double> density_list     = params.getArray<double>("density");
  vector <double>::iterator dens, temp;
  for (dens=density_list.begin();dens != density_list.end(); ++ dens )
    for (temp=temperature_list.begin();temp != temperature_list.end(); ++ temp )
    {
      double this_dens = pow(10,*dens);
      gas.temp = *temp;
      gas.dens = this_dens;
      gas.solve_state();
      gas.computeOpacity(abs_opacity,scat_opacity,emissivity);

      double x;
      for (int i=0;i<nu_grid.size();++i)
      {
        int ind = i;
        if (style == 1) ind = nu_grid.size() - i - 1;
        tot_opacity[ind] = (abs_opacity[ind] + scat_opacity[ind])/this_dens;
        if (density_list.size() > 1) fprintf(fout,"%12.5e ",this_dens);
        if (style == 1) x = pc::c/nu_grid[ind]*1e8;
        else x = nu_grid[ind];
        if (temperature_list.size() > 1) fprintf(fout,"%12.5e ",*temp);
        fprintf(fout,"%15.8e %12.5e\n",x,tot_opacity[ind]);
      }
    }
}


//*********************************************************
// Write out the gas state
//*********************************************************
void write_gas_state(std::string statefile)
{
  FILE *fout = fopen(statefile.c_str(),"w");
  if (fout == NULL)
  {
    if (verbose) std::cerr << "Can't open output gas state file " << statefile << std::endl;
    return;
  }

  fprintf(fout,"---------------------------------------\n");
  fprintf(fout,"density (g/cc)          = %12.4e\n",gas.get_density());
  fprintf(fout,"temperature (K)         = %12.4e\n",gas.get_temperature());
  fprintf(fout,"mean atomic weight      = %12.4e\n",gas.get_mean_atomic_weight());
  fprintf(fout,"number density (g/cc)   = %12.4e\n",gas.get_density()/gas.get_mean_atomic_weight()/pc::m_p);
  fprintf(fout,"electron density (1/cc) = %12.4e\n",gas.get_electron_density());
  fprintf(fout,"---------------------------------------\n\n");

  for (int i=0;i<gas.get_number_atoms();i++)
  {
    fprintf(fout,"--------------------------------------------------\n");
    fprintf(fout,"element Z.A = %d.%d ; mass_fraction =  %.4e\n",
        gas.elem_Z[i],gas.elem_A[i],gas.mass_frac[i]);
    fprintf(fout,"-------------------------------------------------\n");
    fprintf(fout,"  ion   fraction    partition\n");
    for (int j=0;j<gas.get_number_ions(i);j++)
    {
      fprintf(fout,"%4d %12.4e %12.4e\n",j,gas.get_ionization_fraction(i,j),
        gas.get_partition_function(i,j));
    }
  }
}

//*********************************************************
// Write mean opacities
//*********************************************************
void write_mean_opacities(std::string outfile)
{
  FILE *fout = fopen(outfile.c_str(),"w");
  if (fout == NULL)
  {
    if (verbose) std::cerr << "Can't open mean opacity output file " << outfile << std::endl;
    return;
  }

  vector<double> temperature_list = params.getArray<double>("temperature");
  vector<double> density_list     = params.getArray<double>("density");
  int use_logR = params.getScalar<int>("use_logR");

  // write header
  fprintf(fout,"#  temp (K)     rho(g/cc)     ion_x      kappa_p       eps_p        kappa_r      eps_r\n");

  vector <double>::iterator dens, temp;
  for (dens=density_list.begin();dens != density_list.end(); ++ dens )
    for (temp=temperature_list.begin();temp != temperature_list.end(); ++ temp )
    {
      double this_dens = *dens;
      if (use_logR) this_dens = pow(10,*dens + 3*log10(*temp) - 18);
      else this_dens = pow(10,*dens);
      gas.temp = *temp;
      gas.dens = this_dens;
      std::vector<double> J_nu;
      gas.solve_state(J_nu);
      gas.computeOpacity(abs_opacity,scat_opacity,emissivity);

      for (int i=0;i<nu_grid.size();++i)
        tot_opacity[i] = abs_opacity[i] + scat_opacity[i];

      double kr    = rosseland_mean(tot_opacity);
      double eps_r = rosseland_mean(abs_opacity)/kr;
      double kp    = planck_mean(tot_opacity);
      double eps_p = planck_mean(abs_opacity)/kp;
      if (kr == 0) eps_r = 0;
      if (kp == 0) eps_p = 0;

      double ndens = gas.get_density()/gas.get_mean_atomic_weight()/pc::m_p;
      double x_e = gas.get_electron_density()/ndens;
      fprintf(fout,"%12.4e %12.4e %12.4e",gas.temp,gas.dens,x_e);
      fprintf(fout,"%12.4e %12.4e %12.4e %12.4e\n",kp,eps_p,kr,eps_r);
    }
  fclose(fout);
}




//*********************************************************
// Write an opacity table in mesa format
//*********************************************************
void write_mesa_file(std::string mesafile)
{
  FILE *mesaout = fopen(mesafile.c_str(),"w");
  fprintf(mesaout,"Opacities from snopac\n");

    // properties to use
  int use_logR = params.getScalar<int>("use_logR");
  int write_planck_mean = params.getScalar<int>("write_planck_mean");

  vector<double> temperature_list = params.getArray<double>("temperature");
  vector<double> density_list     = params.getArray<double>("density");

  //---------------------------------------------------------
  fprintf(mesaout,"form     version        X           Z          logRs ");
  fprintf(mesaout,"logR_min    logR_max       logTs    logT_min    logT_max\n");

  double X_hydrogen = 0.7;
  double Z_metals   = 0.02;
  fprintf(mesaout,"1   0   %f   %f ",X_hydrogen,Z_metals);
  fprintf(mesaout,"  %lu  %f  %f ",density_list.size(),density_list.front(),density_list.back());
  fprintf(mesaout,"  %lu  %f  %f\n",temperature_list.size(),log10(temperature_list.front()),log10(temperature_list.back()));
  fprintf(mesaout,"logT                       logR = logRho - 3*logT + 18\n");
  fprintf(mesaout,"1   0   %f   %f   %f ",gas.time/(60*60*24),X_hydrogen,Z_metals);

  vector <double>::iterator dens, temp;
  for (dens=density_list.begin();dens != density_list.end(); ++ dens )
    fprintf(mesaout,"%e ",(*dens));
  fprintf(mesaout,"\n\n");

  for (temp=temperature_list.begin();temp != temperature_list.end(); ++ temp )
  {
    fprintf(mesaout,"%f ",log10(*temp));
    for (dens=density_list.begin();dens != density_list.end(); ++ dens )
    {
      double this_dens = *dens;
      if (use_logR)
        this_dens = pow(10,*dens + 3*log10(*temp) - 18);
      else
        this_dens = pow(10,*dens);

      gas.temp = *temp;
      gas.dens = this_dens;
      std::vector<double> J_nu;
      gas.solve_state(J_nu);
      gas.computeOpacity(abs_opacity,scat_opacity,emissivity);

      for (int i=0;i<nu_grid.size();++i)
      {
        tot_opacity[i] = abs_opacity[i] + scat_opacity[i];
      }
      double kr,kp;
      kr = rosseland_mean(tot_opacity);
      if (write_planck_mean)
        kp = planck_mean(tot_opacity);

      if (write_planck_mean) fprintf(mesaout,"%f ",log10(kp));
      else fprintf(mesaout,"%f ",log10(kr));

    }
    fprintf(mesaout,"\n");
  }
  fclose(mesaout);
}

//*********************************************************
// Calculate a rosseland mean opacity
//*********************************************************
double rosseland_mean(std::vector<OpacityType> opac)
{
  // calculate rosseland opacity
  double norm = 0;
  double sum  = 0;
  for (int i = 1;i < nu_grid.size(); ++i)
  {
    double dnu = nu_grid[i] - nu_grid[i-1];
    double nu0 = 0.5*(nu_grid[i] + nu_grid[i-1]);
    double zeta = pc::h*nu0/pc::k/gas.temp;
    double ezeta = exp(zeta);
    double dBdT  = pow(nu0,4);
    dBdT *= ezeta/(ezeta - 1)/(ezeta-1);
    if (isnan(dBdT)) dBdT = 0;
    norm += dBdT*dnu;
    sum  += dBdT/(opac[i])*dnu;
  //  std::cout <<  nu0 << " " << opac[i] << " " << dBdT << "\n";
  }
  double kappa_R = norm/sum/gas.dens;
  return kappa_R;
}

//*********************************************************
// Calculate a planck mean opacity
//*********************************************************
double planck_mean(std::vector<OpacityType> opac)
{
  // calculate planck mean opacity
  double norm = 0;
  double sum  = 0;
  for (int i = 1;i < nu_grid.size(); ++i)
  {
    double dnu = nu_grid[i] - nu_grid[i-1];
    double nu0 = 0.5*(nu_grid[i] + nu_grid[i-1]);
    double zeta = pc::h*nu0/pc::k/gas.temp;
    double Bnu  = pow(nu0,3)/(exp(zeta) - 1);
    if (isnan(Bnu))Bnu = 0;
    norm += Bnu*dnu;
    sum  += Bnu*opac[i]*dnu;
  }
  double kappa_P = sum/norm/gas.dens;
  return kappa_P;
}
*/
