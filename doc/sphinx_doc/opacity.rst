====================
Opacity
====================

There are several options for calculating and controlling the opacity
of the gas in the calculation.


Though we generally use the word "opacity", the code actually calculates and stores an extinction coefficient.
The two are related by

.. math::

  \alpha = \kappa \rho

where
:math:`\alpha` is the extinction coefficient (units cm:math:`^{-1}`),
:math:`\kappa` is the opacity (units cm:math:`^{2}` g:math:`^{-1}`),
and :math:`\rho` is the mass density


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
          - value of grey opacity to use (in cm^2/g). Will override all other opacity settings
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
        * - opacity_line_expansion
          - 0 = no | 1 = yes
          - include binned line expansion opacity
        * - opacity_fuzz_expansion
          - 0 = no | 1 = yes
          - include binned line expansion opacity, taken from a fuzz file

..
	opacity_line_expansion      = 0
	opacity_fuzz_expansion      = 0
	opacity_use_nlte            = 0
	opacity_atoms_in_nlte       = {}
	opacity_minimum_extinction  = 0
	opacity_maximum_opacity     = 1e40
	opacity_no_scattering       = 0
	dont_decay_composition      = 0




-----------------------------------
Grey Opacity
-----------------------------------

The grey opacity flags allow the user to specify a simple,
wavelength-independent opacity. This can be done in two ways.

To set a uniform grey opacity at all points in space, set the parameter::

  opacity_grey_opacity = 0.1

where, in this example, the code will set the
value :math:`\kappa_{g} = 0.1` cm:math:`^{2}` g:math:`^{-1}`
to all zones in the model. All other sources of opacity discussed below
(e.g., free-free, bound-free) will be ignored.

The user can set a spatially varying 

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
