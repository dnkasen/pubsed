============================
Radiation Transport
============================


---------------------------
Controlling Transport
---------------------------

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

Note: The last parameter should really be under tstep, since it controls time evolution
not transport per se.


---------------------------
Controlling Output
---------------------------

Spectrum files and writeout files

-------------------------
Interaction Processes
-------------------------


^^^^^^^^^^^^^^^^^^^^^^^
Electron Scattering
^^^^^^^^^^^^^^^^^^^^^^^


^^^^^^^^^^^^^^^^^^^^^^^
Compton Scattering
^^^^^^^^^^^^^^^^^^^^^^^


^^^^^^^^^^^^^^^^^^^^^^^^^^^
Resonant Line Scattering
^^^^^^^^^^^^^^^^^^^^^^^^^^^



