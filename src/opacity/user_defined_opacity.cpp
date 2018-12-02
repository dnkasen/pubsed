#include "GasState.h"
#include "physical_constants.h"
#include <iostream>

namespace pc = physical_constants;

//----------------------------------------------------------------
// This is a function for
// As an example, uses
// kramer's opacity and e-scattering for pure hydrogen
//----------------------------------------------------------------
void GasState::get_user_defined_opacity
(std::vector<double>& opacity, std::vector<double>&abs_frac, std::vector<double>& emissivity)
{
	int n_freq_pts = nu_grid.size();

	// calculate if is the first time
	if (user_opacity_array_.size() == 0)
	{
		double numax = nu_grid.maxval();
		double numin = nu_grid.minval();
		double nu_space = (numax - numin)/1e4;
		double beta = line_velocity_width_/pc::c;
		double nu_line = numin;
		int bin = 0;

		for (int i=0;i<n_freq_pts;++i)
			opacity[i] = 0.0;


		while (nu_line < numax)
		{
			while ((nu_grid.right(bin) < nu_line)&&(bin < n_freq_pts))
				bin++;
			if (bin >= n_freq_pts) break;

			if (use_user_opacity_ == 2)
			{
				opacity[bin] += 1;
			}
			else
			{
				double dnu = nu_line*beta;
				int di = (int)(dnu/nu_grid.delta(bin)*6.0);
				int i1 = bin - di;
				int i2 = bin + di;
				if (i1 < 0) i1 = 0;
				if (i2 > n_freq_pts-1) i2 = n_freq_pts - 1;

				for (int i=i1;i<i2;++i)
				{
					double nu = nu_grid[i];
					double x = (nu - nu_line)/dnu;
					opacity[i] += exp(-1.0*x*x);
				}
			}
			nu_line += nu_space;
		}

		// save results
		user_opacity_array_.resize(n_freq_pts);
		for (int i=0;i<n_freq_pts;++i)
			user_opacity_array_[i] = opacity[i];
	}

	// set values

	double tau = 100*dens/1e-13;
	double etau = 1 - exp(-tau);
	double beta = line_velocity_width_/pc::c;
	double Ac  = tau/(beta*pc::c*time);

	for (int i=0;i<n_freq_pts;++i)
	{
		double fac = nu_grid.center(i)/nu_grid.delta(i)/pc::c/time;
		if (use_user_opacity_ == 2) opacity[i] = user_opacity_array_[i]*fac*etau;
		else opacity[i] = Ac*user_opacity_array_[i];


		double nu = nu_grid[i];
		double epsilon   = 1.0;
		// thermal emissivity = blackbody*alpha
		double ezeta = exp(1.0*pc::h*nu/pc::k/temp);
		double bb    = 2.0*nu*nu*nu*pc::h/pc::c/pc::c/(ezeta-1);
		double emis  = bb*opacity[i];

		// set the values
		abs_frac[i]   = epsilon;
		emissivity[i] = emis;
	}

}
