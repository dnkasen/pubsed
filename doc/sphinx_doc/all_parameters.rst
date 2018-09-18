===================================
All Runtime Parameters
===================================

|sedona| uses cgs units everywhere

|

.. list-table:: Radiating Core Parameters
	:widths: 15,10,40
	:header-rows: 1
	
	* - parameter
	  - values
	  - definition
	* - core_n_emit
	  - <integer>
	  - Number of particles to emit from core per time step (or iteration)
	* - core_radius
	  - <real>
	  - Radius (in cm) of emitting spherical core
	* - core_luminosity
	  - <real> or <function>
	  - Luminosity (in erg/s) emitted from core
	* - core_temperature 
	  - <real>
	  - Blackbody spectrum of core emission, if using blackbody emission
	* - core_photon_frequency
	  - <real>
	  - Frequency of photons emitted from core, if using monochromatic emission
	* - core_timescale
	  - <real>
	  - ?
	* - core_spectrum_file
	  - <string>
	  - filename of file to read to set spectrum of core emission
	* - core_fix_luminosity
	  - 0 = no | 1 = yes
	  - In steady state calculations, will rescale to fix output luminosity 


|




.. list-table:: Time Stepping Parameters
	:widths: 15,10,40
	:header-rows: 1

	* - parameter
	  - values
	  - definition
   	* - tstep_max_steps
   	  - <integer>
   	  - Maximum number of time steps to take before exiting
   	* - tstep_time_start
   	  - <real>
   	  - Start time (in seconds)
	* - tstep_time_stop
	  - <real>
	  - Stop time (in seconds)
	* - tstep_max_dt
	  - <real>
	  - Maximum value of a time step (in seconds)
	* - tstep_min_dt
	  - <real>
	  - Minimum value of a time step (in seconds)
	* - tstep_max_delta
	  - <real>
	  - Maximum fractional size of a timestep -- restricts dt to the specified value multiplied by the current time

|






.. list-table:: Radiation Transport Parameters
	:widths: 15,10,40
	:header-rows: 1

	* - parameter
	  - values
	  - definition
	* - transport_module
	  - "monte_carlo"
	  - What method to use for transport. Currently only monte carlo is implemented.
	* - transport_nu_grid
	  - <float vector>
	  - Define frequency grid used for transport and opacities
	* - transport_radiative_equilibrium
	  - 0 = no | 1 = yes
	  - Whether to solve for radiative equilibrium
	* - transport_steady_iterate
	  - <integer>
	  - Do a steady-state calculation with this number of iterations

|




.. list-table:: Opacity parameters
        :header-rows: 1
        :widths: 20,10,40

        * - parameter
          - values
          - definition
        * - opacity_grey_opacity        
          - <real>
          - value of grey opacity to use (in cm^2/g). Will overide all other opacity settings
        * - opacity_electron_scattering
          - 0 = no | 1 = yes
          - include electron scattering opacity
        * - opacity_bound_free 
          - 0 = no | 1 = yes
          - include bound-free (photoionization) opacity
        * - opacity_free_free 
          - 0 = no | 1 = yes
          - include free-free opacity
        * - opacity_bound_bound     
          - 0 = no | 1 = yes
          - include bound-bound (resolved line) opacity

|

