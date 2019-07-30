======================
Hydrodynamics
======================



------------------------------
Homologous Expansion
------------------------------



------------------------------
Lagrangian Hydrodynamics
------------------------------




.. list-table:: Hydrodynamics Parameters
        :widths: 15,10,40
        :header-rows: 1

        * - parameter
          - values
          - definition
        * - hydro_module
          - <string>
          - choose hydro module to evolve rho and T, options are homologous, none, 1D_lagrangian
        * - hydro_gamma_index
          - <float>
          - gamma used in the eos
        * - hydro_mean_particle_mass
          - <float>
          - mean particle mass mu, mass = mu * proton mass
        * - hydro_cfl
          - <float>
          - cfl parameter used in hydro simulation to control time steps
        * - hydro_v_piston
          - <float>
          -
        * - hydro_viscosity_parameter
          -
          -
        * - hydro_central_point_mass
          -
          -
        * - hydro_use_gravity
          -
          - Whether to include gravity for the hydro simulation
        * - hydro_use_transport
          - 0 = no | 1 = yes
          - Whether to use radiation transport (?)
        * - hydro_accrete_radius
          -
          -
        * - hydro_bomb_radius
          - <float>
          - Radius of the bomb in cm
        * - hydro_bomb_energy
          - <float>
          - Energy of the bomb in erg
        * - hydro_boundary_outflow
          - 0 = no | 1 = yes
          -
        * - hydro_boundary_rigid_outer_wall
          - 0 = no | 1 = yes
          -
