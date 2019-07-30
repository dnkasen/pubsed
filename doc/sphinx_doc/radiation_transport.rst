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
        * - transport_boundary_in_reflect
          - 0 = no | 1 = yes
          -
        * - transport_boundary_out_reflect
          - 0 = no | 1 = yes
          -
        * - transport_store_Jnu
          - 0 = no | 1 = yes
          -
        * - transport_use_ddmc
          - 0 = no | 1 = yes
          - Whether to use discrete diffusion monte carlo
        * - transport_ddmc_tau_threshold
          - <float>
          - At what optical depth ddmc takes over
        * - transport_fleck_alpha
          - <float>
          - fleck alpha parameter (needs to be between 0.5 and 1 for ddmc)
        * - transport_solve_Tgas_with_updated_opacities
          - 0 = no | 1 = yes
          - whether to solve for Tgas after updating opacities
        * - transport_fix_Tgas_during_transport
          - 0 = no | 1 = yes
          - whether to fix Tgas
        * - transport_set_Tgas_to_Trad
          - 0 = no | 1 = yes
          - whether to set Tgas to Trad instead of solving for it




Note: The parameter transport_steady_iterate
.. ?
should really be under tstep, since it controls time evolution
not transport per se.


          -

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
