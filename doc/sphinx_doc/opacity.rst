====================
Opacity
====================


-----------------------
Atomic Data Files
-----------------------

Atomic data is store in hdf5 files

-----------------------
Opacity Settings
-----------------------


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

..
	opacity_line_expansion      = 0
	opacity_fuzz_expansion      = 0
	opacity_use_nlte            = 0
	opacity_atoms_in_nlte       = {}
	opacity_minimum_extinction  = 0
	opacity_maximum_opacity     = 1e40
	opacity_no_scattering       = 0
	dont_decay_composition      = 0


Though we use the word "opacity", the code actually calculates and stores an extinction coefficient.
The two are related by

.. math::

  \alpha = \kappa \rho

where...

-----------------------------------
Electron Scattering Opacity
-----------------------------------

This is calculated as

.. math ::

  \alpha_{\rm es} = \sigma_t n_e

where :math:`n_e` is the free electron density

Options for include Klein-Nishina corrections, Comptonization etc...

-----------------------------------
Bound-Free Opacity
-----------------------------------


-----------------------------------
Free-Free Opacity
-----------------------------------


-----------------------------------
Line Opacity
-----------------------------------

There are various options for treating lines. In general, only one of these approaches
should be used to treat lines, unless one can be certain that line opacity is not be
multiply counted.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Resolved Line Profiles
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This approach is selecting by setting the runtime parameter::

  opacity_bound_bound = 1

This approach treats the lines generally as Voigt profiles. The frequency good must be
fine enough that there are multiple grid points across to resolve the line profile.

The widths of lines can be artifically broadened using the runtime paramater::

  line_velocity_width = <real>

where <real> is a velocity (in cm/s) by which the lines should be Gaussian broadened.


^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Line Expansion Opacity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This approach is selected by setting the runtime parameter::
  
  opacity_line_expansion = 1

This approach bins lines into frequency bins, assuming a homologous flow.


This approach implies the Sobolev approximation. For each line, the code calculates the Sobolev optical depth

.. math::

  \tau_{\rm sob} = \frac{ \pi e^2}{m_e c} n_l f_{lu} \lambda_0 t_{\rm exp}

where...

The expansion opacity is then calculated by binning lines

.. math::

  \alpha_{\rm exp} = \frac{1}{c t_{\rm exp}} \frac{\lambda_0}{\Delta \lambda} \sum_i (1 - e^{-\tau_i})



^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Fuzz Expansion Opacity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This approach is selected by setting the runtime parameter::
  
  opacity_fuzz_expansion = 1

The physics here is identical to that of line expansion opacity, it is
just that the lines are read from an independent file.



^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Resonant Line Scattering
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is how we would treat scattering in strong individual lines like Lyman alpha.

-----------------------------------
LTE and NLTE settings
-----------------------------------

