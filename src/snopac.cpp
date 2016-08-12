#include <iostream>
#include <fstream>
#include "nlte_gas.h"
#include "locate_array.h"
#include "physical_constants.h"
#include "ParameterReader.h"

#include <mpi.h>

using namespace std;
namespace pc = physical_constants;

locate_array nu_grid;
nlte_gas gas;

int main(int argc, char **argv)
{
  double rosseland_mean(std::vector<double> opac);

  // initialize MPI 
  MPI_Init( &argc, &argv );
  const int verbose = 0; 

  // open up the parameter reader
  std::string param_file = "param.lua";
  if( argc > 1 ) param_file = std::string( argv[ 1 ] );
  ParameterReader params(param_file,verbose);

  // initialize the frequency grid  
  std::vector<double> nu_dims = params.getVector<double>("transport_nu_grid");
  if ((nu_dims.size() != 4)&&(nu_dims.size() != 3)) {
    cout << "# improperly defined nu_grid; need {nu_1, nu_2, dnu, (log?)}; exiting\n";
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
  std::vector<double> elements = params.getVector<double>("elements");
  std::vector<int> elems_A, elems_Z;
  for (int i=0;i<elements.size();i++)  
  {
    std::string species = std::to_string(elements[i]);
    int ind1 = species.find(".");
    int ind2 = species.find("0");
    std::string el_Z = species.substr(0,ind1);
    std::string el_A = species.substr(ind1+1,ind2-ind1-1);
    elems_Z.push_back(std::stoi(el_Z));
    elems_A.push_back(std::stoi(el_A));
  }
  //-------------------------------------------------------

  // initialize nlte_gas class
  std::string atomdata = params.getScalar<string>("data_atomic_file");  
  gas.initialize(atomdata,elems_Z,elems_A,nu_grid);
  std::string fuzzfile = params.getScalar<string>("data_fuzzline_file");  
  gas.read_fuzzfile(fuzzfile);
  std::vector<double> massfrac = params.getVector<double>("mass_fractions");
  gas.set_mass_fractions(massfrac);

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

  // properties to use
  std::vector<double> dens_arr = params.getVector<double>("density");
  std::vector<double> temp_arr = params.getVector<double>("temperature");
  int use_logR = params.getScalar<int>("use_logR");

  double time = params.getScalar<double>("time");

  double dstart,dstop,ddel;
  if (dens_arr.size() == 1) 
  {
    dstart = dens_arr[0];
    dstop  = dens_arr[0];
    ddel   = 1;
  }
  else if (dens_arr.size() == 2)
  {
    std::cout << "wrong density array\n";
    exit(1);
  }
  else 
  {
    dstart = dens_arr[0];
    dstop  = dens_arr[1];
    ddel   = dens_arr[2];
  }

  double tstart,tstop,tdel;
  if (temp_arr.size() == 1) 
  {
    tstart = temp_arr[0];
    tstop  = temp_arr[0];
    tdel   = 1;
  }
  else if (temp_arr.size() == 2)
  {
    std::cout << "wrong temperature array\n";
    exit(1);
  }
  else
  {
    tstart = temp_arr[0];
    tstop  = temp_arr[1];
    tdel   = temp_arr[2];
  }


  // calculate the opacities/emissivities
  std::vector<double> abs_opacity, scat_opacity, tot_opacity, emissivity;
  int ng = nu_grid.size();
  abs_opacity.resize(ng);
  scat_opacity.resize(ng);
  tot_opacity.resize(ng);
  emissivity.resize(ng);

  for (double dens = dstart; dens <= dstop; )
  {
    for (double temp = tstart; temp <= tstop; )
    {
      double this_dens = dens;
      if (use_logR)
        this_dens = pow(10,dens + 3*log10(temp) - 18);

      gas.time = time;
      gas.temp = temp;
      gas.dens = this_dens;
      std::vector<double> J_nu;
      gas.solve_state(J_nu);
      gas.computeOpacity(abs_opacity,scat_opacity,emissivity);

      for (int i=0;i<ng;i++)
      {
        tot_opacity[i] = abs_opacity[i] + scat_opacity[i];
        //std::cout << dens << "\t" << temp << "\t";
        //std::cout << nu_grid[i] << "\t" << tot_opacity[i]/dens << "\n";
      }
      double kr = rosseland_mean(tot_opacity);
      std::cout << this_dens << "\t" << temp << "\t" << kr << "\n";

      if (temp_arr.size() == 4) temp = temp*(1 + tdel);
      else temp += tdel;

    }
    if (dens_arr.size() == 4) dens = dens*(1 + ddel);
    else dens += ddel;
  }

}

double rosseland_mean(std::vector<double> opac)
{
  // calculate rosseland opacity
  double norm = 0;
  double sum  = 0;
  for (int i = 1;i < nu_grid.size(); i++)
  {
    double dnu = nu_grid[i] - nu_grid[i-1];
    double nu0 = 0.5*(nu_grid[i] + nu_grid[i-1]);
    double zeta = pc::h*nu0/pc::k/gas.temp;
    double ezeta = exp(zeta);
    double dBdT  = 2*pc::h*pow(nu0,4)/(pc::c*pc::c*pc::k*gas.temp*gas.temp);
    dBdT *= ezeta/(ezeta - 1)/(ezeta-1);
    if (isnan(dBdT)) dBdT = 0;
    norm += dBdT*dnu;
    sum  += dBdT/(opac[i])*dnu;
  //  std::cout <<  nu0 << " " << opac[i] << " " << dBdT << "\n";
  }
  double kappa_R = norm/sum/gas.dens;
  return kappa_R;
}

  // std::vector<double> line_opac;
  // std::vector<double> fuzz_opac;
  // fuzz_opac.resize(n_pts);
  // line_opac.resize(n_pts);
  // gas.line_expansion_opacity(line_opac);
  // gas.fuzz_expansion_opacity(fuzz_opac);

  //ofstream outfile("Fe_exp_opac.txt");
  // for (int i=0;i<n_pts;i++)
  // {
  //   double lam = pc::c/(nu_grid[i]*pc::angs_to_cm);
  //   outfile << lam << "\t" << line_opac[i]/dens << " ";
  //   outfile << fuzz_opac[i]/dens << "\n";
  // }
  //outfile.close();
  



  
