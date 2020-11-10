====================
Temperature Solver
====================

There are several parameters that can be set in the .lua file for controlling how the gas temperature is updated as part of a transport step.

Several combinations of these parameters may conflict with each other. The code attempts to detect these conflicts near the start of execution, but it has not been tested in every case. Care should be taken whenever modifying these parameters.

....


.. code-block:: c

	transport_fix_Tgas_during_transport

This can prevent the gas temperature from being modified after particles are propagated in a transport step. Preventing such a temperature update can be useful when debugging problems where you want to use a known temperature, so that you can focus on other aspects of the calculation. This flag might also be needed for some versions of hydro calculations using [sedona] where you want to update the temperature in the hydro module based on absorbed radiation, but don't want the temperature to be updated based on radiative equilibrium or some other assumption during the transport step.

The default value for this paramaeter is 0, meaning the temperature *can* get updated during the transport step if you have otherwise instructed the code to do so. For example, you may have set transport_radiative_equilibrium = 1, which would ordinarily cause the temperature to be udpated shortly after particles are propagated during each transport step. Setting transport_fix_Tgas_during_transport parameter = 1 will prevent the temperature from being updated that way.

Setting the transport_fix_Tgas_during_transport parameter to 1 is not the same as setting the transport_radiative_equilibrium parameter to 0. That second flag controls not only whether to use radiative equlibrium to update the gas temperature during transport steps, but also whether to effectively scatter particles whenever they would be absorbed. Thus, if you want particles to be effectively scattered, but you don't want the temperature to be updated, you would set transport_radiative_equilibrium=1, while also setting transport_fix_Tgas_during_transport=1.

....


.. code-block:: c

	transport_set_Tgas_to_Trad

This determines whether to set the gas temperature to the radiation "temperature" as computed from the radiation energy density, as part of the transport step. This is an alternative to setting the temperature during the transport step based on radiative equilibrium. Using the radiation temperature leads to a simpler but usually less accurate gas temperature solution than by solving for the temperature by assuming radiative equilibrium.

The default value of this parmeter is 0, which means it is not used. It is activated by setting the parameter to 1.

....


.. code-block:: c

	transport_solve_Tgas_with_updated_opacities


This parameter controls whether temperature-dependent updates of the gas state are included in each iteration of the temperature solution. This can be helpful for ensuring stability in the temperature solution. Since this requires re-computing aspects of the gas state many times, this can signficantly slow down the calculation.

By default this parameter is set to 0, which means it is not used. When this parameter is is set to 1, its effect will depend on whether the calcultaion is LTE or NLTE. In LTE, activating this flag means that, during each iteration of the temperature solver, the absorption opacity will be re-computed with the most recent guess for Tgas. In NLTE, activating this flag means that the gas ionization state and bound electron populations will be re-computed updated during each iteration of the temperature solution - that is, the entire NLTE matrix equation, under the constraint of charge conservation, will be re-solved during each iteration of the temperature solver.
