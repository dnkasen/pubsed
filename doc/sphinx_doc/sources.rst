===================
Radiation Sources
===================

-------------------
Radioactivity
-------------------


-------------------
Thermal Emission
-------------------

-------------------
Radiating Core
-------------------


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
        * - particles_max_total
          - <float>
          - maximum number of particles (photons) allowed on the grid at the same time
        * - particles_n_emit_radioactive
          - <integer>
          - number of particles emitted through radioactivity per timestep
        * - particles_n_emit_thermal
          - <integer>
          - number of thermal particles emitted per timestep
        * - particles_n_initialize
          - <integer>
          - number of particles used to initialize the simulation
        * - particles_n_emit_pointsources
          - <integer>
          -
        * - particles_pointsource_file
          - <string>
          -
        * - particles_last_iter_pump
          -
          -
        * - multiply_particles_n_emit_by_dt_over_dtmax
          - 0 = no | 1 = yes
          -
        * - force_rprocess_heating
          - 0 = no | 1 = yes
          -




-------------------------
Multiple Point Sources
-------------------------
