.. _pop-params-guide:

================================================
POSDYON Population Synthesis Configuration Guide
================================================

This documentation provides a detailed overview of the configuration options available in the Posydon software package.


.. warning::
    The default values in the population_params.ini file included in POSYDON 
    have not been calibrated, validated, or are used by the POSYDON team. They are
    often an ad-hoc choice, and the user should carefully consider the values
    of each parameter for their science case.


Environment Variables
---------------------

The environment variable `PATH_TO_POSYDON` will be read from your shell session.

.. list-table::
   :widths: 50 50
   :header-rows: 1

   * - Parameter
     - Description/Value
   * - `PATH_TO_POSYDON`
     - `<PATH_TO_POSYDON>`


SimulationProperties
--------------------

After the environment variables, the next part of the ``population_params.ini`` file
contains how systems will move through the simulation steps [Flow Chart](flow_chart.rst) and the parameters for each step.

These parameters are read in and used to create a :class:`~posydon.binary_evol.simulationproperties.SimulationProperties`` object.

It allows for the addition of hooks before and after each step, and the evolution of a system (see :ref:`Custom Hooks <custom-hooks>`).
The ``SimulationProperties`` can be manually read and loaded in.
Note that the loading of the simulation steps is separate from creating the object.

.. code-block:: python

    from posydon.binary_evol.simulationproperties import SimulationProperties
    from posydon.popsyn.io import simprop_kwargs_from_ini

    # read from file
    sim_props = simprop_kwargs_from_ini('population_params.ini', vebose=False)
    
    # create SimulationProperties object
    sim = SimulationProperties(**sim_props)

    # load the steps
    sim.load_steps()


Flow Chart
~~~~~~~~~~

The flow chart is the core of POSYDON.
It controls the mapping between a POSYDON binary object and its step evolution, see the :ref:`Flow Chart Object <flow-chart>` page for more details.

.. list-table::
  :widths: 10 80 10
  :class: population-params-table
  :header-rows: 1
  

  * - Parameter
    - Description
    - Default Value

  * - ``import``
    - | The import path for the flow chart and the name of the flow chart function.
      | The flow chart function is used to determine the next step of the binary object.
    - ``['posydon.binary_evol.flow_chart', 'flow_chart']``
 
  * - ``absolute_import``
    - | An absolute import of a custom step. It follows the same structure as ``import``.
      | You can find an example in the tutorials.
    - ``None``

Step MESA (HMS-HMS, CO-HMS_RLO, CO-HeMS, CO-HeMS_RLO)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The MESA step is the most important step of POSYDON as it leverages the POSYDON MESA grids to evolve the binary object according to one of the supported MESA binary-star grids.

Below we give the default values for the HMS-HMS step, as an example.

.. list-table::
  :widths: 10 80 10
  :class: population-params-table
  :header-rows: 1
   
  * - Parameter
    - Description
    - Default Value
    
  * - ``import``
    - | The import path for the step and the name of the step class.
    - ``['posydon.binary_evol.DT.step_MESA', 'step_HMS_HMS']``

  * - ``absolute_import``
    - | An absolute import of a custom step. It follows the same structure as ``import``. 
      | ``['import.path', 'name_of_step']``
    - ``None``

  * - ``interpolation_path``
    - | The path to the interpolation file for the step. 
      | If None, the path is found by default with
      | ``PATH_TO_POSYDON_DATA`` and ``step_name``.
    - ``None``

  * - ``interpolation_filename`` 
    - | The name of the interpolation file for the step.
      | If None, the name is found by default with
      | ``PATH_TO_POSYDON_DATA``, ``step_name``, and ``metallicity``.
    - ``None``

  * - ``interpolation_method``
    - | The interpolation method used to interpolate the MESA grid.
     
      * ``'nearest_neighbour'`` : nearest neighbour interpolation
      * ``'linear3c_kNN'`` : linear interpolation with closest neighbours
      * ``'1NN_1NN'`` : 1NN interpolation with 1NN interpolation?
    - ``'linear3c_kNN'``

  * - ``save_initial_conditions``
    -  | Store the initial conditions when using nearest neighbour interpolation.
    - ``True``

  * - ``track_interpolation``
    - | Track the interpolation when using nearest neighbour interpolation.
    - ``False``

  * - ``stop_method``
    - | The method to stop the evolution of a binary
      
      * ``'stop_at_end'`` : stop at the end of the simulation
      * ``'stop_at_max_time'`` : stop at the maximum time
      * ``'stop_at_condition'`` : stop at a condition

    - ``'stop_at_max_time'``

  * - ``stop_star``
    - | Specifies the star to stop for the condition.  
      | Relevant when ``stop_method`` is ``'stop_at_condition'``.

      * ``'star_1'`` : stop condition applied to star 1
      * ``'star_2'`` : stop condition applied to star 2
    
    - ``'star_1'``

  * - ``stop_var_name``
    - | The variable name for the stop condition. 
      | Only applicable when ``stop_method`` is ``'stop_at_condition'``.
    - ``None``

  * - ``stop_value``
    - | The value at which to stop. 
      | Relevant when ``stop_method`` is ``'stop_at_condition'``.
    - ``None``

  * - ``stop_interpolate``
    - | Specifies whether to interpolate when stopping.
    - ``True``
  
  * - ``verbose``
    - | Enables verbose mode.
    - ``False``



Step Detached
~~~~~~~~~~~~~

The detached step uses analytical expressions and MESA single star grids to evolve the binary object when it is not interacting.
It evolves the binary object in isolation until Roche lobe overflow occurs.

.. list-table::
  :widths: 10 80 10
  :class: population-params-table
  :header-rows: 1
   
  * - Parameter
    - Description
    - Default Value

  * - ``import``
    - | The import path for the step and the name of the step class.
    - ``['posydon.binary_evol.DT.step_detached', 'detached_step']``

  * - ``absolute_import``
    - | An absolute import of a custom step. It follows the same structure as ``import``. 
      | ``['import.path', 'name_of_step']``
    - ``None``
  
  * - ``matching_method``
    - | The method to match the MESA single star grid to the binary object.
      
      * ``'minimize'`` : minimize the difference between the MESA single star grid and the binary object
      * ``'root'`` : find the root of the difference between the MESA single star grid and the binary object

    - ``'minimize'``

  * - ``do_wind_loss``
    - | Enable wind loss.
    - ``True``

  * - ``do_tides``
    - | Enable tides.
    - ``True``

  * - ``do_gravitational_radiation``
    - | Enable gravitational radiation.
    - ``True``
  
  * - ``do_magnetic_braking``
    - | Enable magnetic braking.
    - ``True``

  * - ``do_stellar_evolution_and_spin_from_winds``
    - | Enable stellar evolution and spin from winds.
    - ``True``

  * - ``RLO_orbit_at_orbit_with_same_am``
    - | Orbit at Roche lobe overflow with the same angular momentum.
      | Relevant in eccentric systems.
    - ``False``

  * - ``verbose``
    - | Enables verbose mode.
    - ``False``
    
Step Disrupted
~~~~~~~~~~~~~~

The dirtupted step evolves a system when a supernova has unbound the binary components.
This class inherits from the detached step, but only evolves the remaining star in isolation.
This means that this step uses the single stars loaded by the detached step.

.. list-table::
  :widths: 10 80 10
  :class: population-params-table
  :header-rows: 1

  * - Parameter
    - Description
    - Default Value

  * - ``import``
    - | The import path for the step and the name of the step class.
    - ``['posydon.binary_evol.DT.step_disrupted', 'DisruptedStep']``

Step Merged
~~~~~~~~~~~

In this step, the system has undergone a merger and is now a single star.
This class inherits from the detached step, but only evolves the remaining star in isolation.
This means that this step uses the single stars loaded by the detached step.

.. list-table::
  :widths: 10 80 10
  :class: population-params-table
  :header-rows: 1

  * - Parameter
    - Description
    - Default Value

  * - ``import``
    - | The import path for the step and the name of the step class.
    - ``['posydon.binary_evol.DT.step_merged', 'MergedStep']``

Step Initially Single
~~~~~~~~~~~~~~~~~~~~~

This step is used to evolve a single star system.
This class inherits from the detached step and evolves the star in isolation.
This means that this step uses the single stars loaded by the detached step.


.. list-table::
  :widths: 10 80 10
  :class: population-params-table
  :header-rows: 1

  * - Parameter
    - Description
    - Default Value

  * - ``import``
    - | The import path for the step and the name of the step class.
    - ``['posydon.binary_evol.DT.step_initially_single','InitiallySingleStep']``


Step Common Envelope
~~~~~~~~~~~~~~~~~~~~

The common envelope step is used to evolve a binary system when the primary or secondary star has initiated a common envelope phase.
It calculates the binding energy of the envelope and the energy available to eject it.
If the energy budget is greater than the binding energy, the envelope is ejected and the system is evolved to the next step.
If the energy budget is less than the binding energy, the system is merged.

.. list-table::
  :widths: 10 80 10
  :class: population-params-table
  :header-rows: 1

  * - Parameter
    - Description
    - Default Value

  * - ``import``
    - | The import path for the step and the name of the step class.
    - ``['posydon.binary_evol.CE.step_CEE', 'StepCEE']``

  * - ``absolute_import``
    - | An absolute import of a custom step. It follows the same structure as ``import``.
      | ``['import.path', 'name_of_step']``
    - ``None``

  * - ``prescription``
    - | The prescription used for the common envelope evolution.
    - ``'alpha-lambda'``

  * - ``common_envelope_efficiency``
    - | The efficiency of the common envelope ejection.
    - ``1.0``

  * - ``common_envelope_option_for_lambda``
    - | The option for calculating the lambda parameter.
     
      * ``'default_lambda'`` : default lambda value
      * ``'lambda_from_grid_final_values'`` : lambda from grid final values
      * ``'lambda_from_profile_gravitational'`` : lambda from profile gravitational
      * ``'lambda_from_profile_gravitational_plus_internal'`` : lambda from profile gravitational plus internal
      * ``'lambda_from_profile_gravitational_plus_internal_minus_recombination'`` : lambda from profile gravitational plus internal minus recombination
    - ``'lambda_from_grid_final_values'``

  * - ``common_envelope_lambda_default``
    - | The default lambda value, used only for the ``'default_lambda'`` option.
    - ``0.5``

  * - ``common_envelope_option_for_HG_star``
    - | The option for handling Hertzsprung gap stars in common envelope evolution.
     
      * ``'optimistic'`` : optimistic scenario
      * ``'pessimistic'`` : pessimistic scenario
    - ``'optimistic'``

  * - ``common_envelope_alpha_thermal``
    - | The alpha thermal parameter, used only for the ``'lambda_from_profile_gravitational_plus_internal'`` and ``'lambda_from_profile_gravitational_plus_internal_minus_recombination'`` options.
    - ``1.0``

  * - ``core_definition_H_fraction``
    - | The hydrogen fraction for defining the core.
     
      * ``0.01``
      * ``0.1``
      * ``0.3``
    - ``0.1``

  * - ``core_definition_He_fraction``
    - | The helium fraction for defining the core.
    - ``0.1``

  * - ``CEE_tolerance_err``
    - | The tolerance error for the common envelope evolution.
    - ``0.001``

  * - ``common_envelope_option_after_succ_CEE``
    - | The option for handling the system after a successful common envelope ejection.
     
      * ``'core_not_replaced_noMT'`` : core not replaced, no mass transfer
      * ``'core_replaced_noMT'`` : core replaced, no mass transfer
      * ``'core_not_replaced_stableMT'`` : core not replaced, stable mass transfer
      * ``'core_not_replaced_windloss'`` : core not replaced, wind loss
    - ``'core_not_replaced_noMT'``

  * - ``verbose``
    - | Enables verbose mode.
    - ``False``


Step Supernova
~~~~~~~~~~~~~~

This step performs the core-collapse supernova evolution of the star.
Multiple prescriptions are implemented and can be selected.
If a pretrained prescription is used (``use_interp_values=True``), interpolators are available to speed up the calculation and predict additional properties of the collapse, such as the spin.
The collection of trained prescriptions can be found in the ``MODELS.py`` file at ``posydon/grids/MODELS.py``.

.. warning::
  If you want to use a combination that does not have a trained model, please turn off ``use_interp_values``,
  otherwise your population runs will fail! Depending on ``use_core_masses``, 
  the core masses or stellar profiles are used to determine the supernova properties.





.. list-table::
  :widths: 10 80 10
  :class: population-params-table
  :header-rows: 1

  * - Parameter
    - Description
    - Default Value

  * - ``import``
    - | The import path for the step and the name of the step class.
    - ``['posydon.binary_evol.SN.step_SN', 'StepSN']``

  * - ``absolute_import``
    - | An absolute import of a custom step. It follows the same structure as ``import``.
      | ``['import.path', 'name_of_step']``
    - ``None``

  * - ``mechanism``
    - | The mechanism used for the supernova.
     
      * ``'direct'``
      * ``'Fryer+12-rapid'``
      * ``'Fryer+12-delayed'``
      * ``'Sukhbold+16-engine'``
      * ``'Patton&Sukhbold20-engine'``
    - ``'Patton&Sukhbold20-engine'``

  * - ``engine``
    - | The engine used for the supernova.
      | Relevant for ``'Sukhbold+16-engine'`` and ``'Patton&Sukhbold20-engine'`` mechanisms.
    - ``'N20'``

  * - ``PISN``
    - | The prescription used for pair-instability supernova.
     
      * ``None``
      * ``'Marchant+19'``
    - ``'Marchant+19'``

  * - ``ECSN``
    - | The prescription used for electron-capture supernova.
     
      * ``'Tauris+15'``
      * ``'Podsiadlowksi+04'``
    - ``'Podsiadlowksi+04'``

  * - ``conserve_hydrogen_envelope``
    - | Conserve the hydrogen envelope during the supernova.
    - ``True``

  * - ``max_neutrino_mass_loss``
    - | The maximum mass loss due to neutrinos.
    - ``0.5``

  * - ``max_NS_mass``
    - | The maximum mass of a neutron star.
    - ``2.5``

  * - ``use_interp_values``
    - | Use interpolation values for the supernova properties, which are calculated during ``step_MESA``.
    - ``True``

  * - ``use_profiles``
    - | Use profiles for the supernova.
    - ``True``

  * - ``use_core_masses``
    - | Use core masses for the supernova.
    - ``True``

  * - ``approx_at_he_depletion``
    - | Approximate at helium depletion.
    - ``False``

  * - ``kick``
    - | Apply a kick to the remnant.
    - ``True``

  * - ``kick_normalisation``
    - | The normalisation method for the kick.
     
      * ``'one_minus_fallback'``
      * ``'one_over_mass'``
      * ``'NS_one_minus_fallback_BH_one'``
      * ``'one'``
      * ``'zero'``
    - ``'one_over_mass'``

  * - ``sigma_kick_CCSN_NS``
    - | The standard deviation of the kick velocity for core-collapse supernova neutron stars.
    - ``265.0``

  * - ``sigma_kick_CCSN_BH``
    - | The standard deviation of the kick velocity for core-collapse supernova black holes.
    - ``265.0``

  * - ``sigma_kick_ECSN``
    - | The standard deviation of the kick velocity for electron-capture supernova.
    - ``20.0``

  * - ``verbose``
    - | Enables verbose mode.
    - ``False``

Step Double Compact Object
~~~~~~~~~~~~~~~~~~~~~~~~~~

In this step, the system has evolved to a double compact object system.
The merger time due to gravitational wave emission is calculated.

.. list-table::
  :widths: 10 80 10
  :class: population-params-table
  :header-rows: 1

  * - Parameter
    - Description
    - Default Value

  * - ``import``
    - | The import path for the step and the name of the step class.
    - ``['posydon.binary_evol.DT.double_CO', 'DoubleCO']``

  * - ``absolute_import``
    - | An absolute import of a custom step. It follows the same structure as ``import``.
      | ``['import.path', 'name_of_step']``
    - ``None``

  * - ``n_o_steps_interval``
    - | The number of steps interval for the double compact object evolution.
    - ``None``

Step End
~~~~~~~~

The final step of the simulation.
All succesfull systems should reach this step.

.. list-table::
  :widths: 10 80 10
  :class: population-params-table
  :header-rows: 1

  * - Parameter
    - Description
    - Default Value

  * - ``import``
    - | The import path for the step and the name of the step class.
    - ``['posydon.binary_evol.step_end', 'step_end']``

  * - ``absolute_import``
    - | An absolute import of a custom step. It follows the same structure as ``import``.
    - ``None``

Extra Hooks
~~~~~~~~~~~

It's possible to add custom hooks to the simulation steps.
A few example hooks are provides: ``TimingHooks`` and ``StepNamesHooks`` (See :ref:`Custom Hooks <custom-hooks>` for more details)
Each new hook has a unique import number, starting from 1.

.. list-table::
  :widths: 10 80 10
  :class: population-params-table
  :header-rows: 1

  * - Parameter
    - Description
    - Default Value

  * - ``import_1``
    - | The import path for the hook and the name of the hook class.
    - ``['posydon.binary_evol.simulationproperties', 'TimingHooks']``

  * - ``absolute_import_1``
    - | An absolute import of a custom hook. It follows the same structure as ``import``.
    - ``None``

  * - ``kwargs_1``
    - | Additional keyword arguments for the hook.
    - ``{}``

  * - ``import_2``
    - | The import path for the hook and the name of the hook class.
    - ``['posydon.binary_evol.simulationproperties', 'StepNamesHooks']``

  * - ``absolute_import_2``
    - | An absolute import of a custom hook. It follows the same structure as ``import``.
    - ``None``

  * - ``kwargs_2``
    - | Additional keyword arguments for the hook.
    - ``{}``



BinaryPopulation
----------------

A :class:`~posydon.popsyn.binarypopulation.BinaryPopulation` is created to evolve the population using a given :class:`~posydon.binary_evol.simulationproperties.SimulationProperties` object.

This class requires additional parameters, because it will require initial distributions 
for to sample :class:`~posydon.binary_evol.binarystar.BinaryStar` objects from,
such as the masses and orbital parameters.
Moreover, it contains the parameters for metallicity and the practicality of running populations.
This includes, the number of binaries, the metallicity, how often to save the population to file. 

When reading the binary population arguments from a ``population_params.ini`` file, the
 :class:`~posydon.binary_evol.simulationproperties.SimulationProperties` are read in automatically.

.. code-block:: python

    from posydon.popsyn.binarypopulation import BinaryPopulation
    from posydon.popsyn.io import binarypop_kwargs_from_ini

    # read from file
    pop_params = binarypop_kwargs_from_ini('population_params.ini', vebose=False)
    
    # create BinaryPopulation object
    pop = BinaryPopulation(**pop_params)



BinaryPopulation Options
~~~~~~~~~~~~~~~~~~~~~~~~

These parameters contain options on how the population is evolved, in practical terms, ie. the number of binaries.
It also contains which sampling distributions to use for the initial conditions of the binaries.

.. list-table::
  :widths: 10 80 10
  :class: population-params-table
  :header-rows: 1

  * - Parameter
    - Description
    - Default Value

  * - ``optimize_ram``
    - | Save population in batches.
    - ``True``

  * - ``ram_per_cpu``
    - | Set maximum RAM per CPU before batch saving (in GB).
    - ``None``

  * - ``dump_rate``
    - | Batch save after evolving N binaries.
    - | To facilitate I/O performance, this should not be set below 500 for populations of 100.000 binaries or more. 
    - ``2000``

  * - ``temp_directory``
    - | Folder for keeping batch files.
    - ``'batches'``

  * - ``tqdm``
    - | Enable progress bar.
    - ``False``

  * - ``breakdown_to_df``
    - | Convert BinaryStars into DataFrames after evolution.
    - ``True``

  * - ``use_MPI``
    - | If True, evolve with MPI (equivalent to: ``from mpi4py import MPI, comm = MPI.COMM_WORLD``).
    - ``True``

  * - ``metallicity``
    - | In units of solar metallicity. Supported values ``[2., 1., 0.45, 0.2, 0.1, 0.01, 0.001, 0.0001]``
    - ``[1.]``

  * - ``error_checking_verbose``
    - | If True, write all POSYDON errors to stderr at runtime
    - ``False``

  * - ``warnings_verbose``
    - | If True, write all POSYDON warnings to stderr at runtime
    - ``False``
    
  * - ``history_verbose``
    - | If True, record extra functional steps in the output DataFrames
    - | (These extra steps represent internal workings of POSYDON rather than physical phases of evolution)
    - ``False``

  * - ``entropy``
    - | Random Number Generation: uses system entropy.
    - ``None``

  * - ``number_of_binaries``
    - | Number of binaries to evolve.
    - ``10``

  * - ``star_formation``
    - | What star formation is used to sample the binaries.
      | Options:
    
      * ``'constant'``: sample with a constant rate over time
      * ``'burst'``: sample from a burst of star formation
      * ``'custom_linear'``: sample with a custom linear rate over time
      * ``'custom_log10'``: sample with a custom log10 rate over time
      * ``'custom_linear_histogram'``: sample from a custom linear histogram rate over time
      * ``'custom_log10_histogram'``: sample from a custom log10 histogram rate over time
    - ``'burst'``

  * - ``max_simulation_time``
    - | Maximum simulation time (in years).
    - ``13.8e9``

  * - ``binary_fraction``
    - | Fraction of binaries (0 < fraction <= 1).
    - ``1.0``

  * - ``primary_mass_scheme``
    - | Options:
    
      * ``Salpeter``: `Salpeter E. E., 1955, ApJ, 121, 161 <https://ui.adsabs.harvard.edu/abs/1955ApJ...121..161S/abstract>`_
      * ``'Kroupa1993'``: `Kroupa P., Tout C. A., Gilmore G., 1993, MNRAS, 262, 545 <https://ui.adsabs.harvard.edu/abs/1993MNRAS.262..545K/abstract>`_
      * ``'Kroupa2001'``: `Kroupa P., 2001, MNRAS, 322, 231 <https://ui.adsabs.harvard.edu/abs/2001MNRAS.322..231K/abstract>`_
    - ``'Kroupa2001'``

  * - ``primary_mass_min``
    - | Minimum primary mass (in solar masses)
      | limits: 0-300
    - ``7.0``

  * - ``primary_mass_max``
    - | Maximum primary mass (in solar masses).
      | Needs to be larger than ``secondary_mass_min``.
      | limits: 0-300
    - ``150.0``

  * - ``secondary_mass_scheme``
    - | Options:
    
      * ``'flat_mass_ratio'``: flat mass ratio distribution
      * ``'q=1'``: mass ratio of 1. Ignores ``secondary_mass_min/max`` and sets the secondary mass to the primary mass. 
    - ``'flat_mass_ratio'``

  * - ``secondary_mass_min``
    - | Minimum secondary mass (in solar masses).
      | Is required to be smaller than the minimum primary mass.
      | limits: 0-270
    - ``0.35``

  * - ``secondary_mass_max``
    - | Maximum secondary mass (in solar masses).
      | limits: 0-270
    - ``150.0``

  * - ``orbital_scheme``
    - | How to the orbital parameter is sampled.
      | Options:
    
      * ``'separation'``: use orbital separation
      * ``'period'``: use orbital period
    - ``'period'``

  * - ``orbital_period_scheme``
    - | Used only for ``orbital_scheme = 'period'``.
      | Options:
    
      * ``Sana+12_period_extended``: `Sana et al. 2012 <https://ui.adsabs.harvard.edu/abs/2012Sci...337..444S/abstract>`_
    - ``'Sana+12_period_extended'``

  * - ``orbital_period_min``
    - | Minimum orbital period (in days).
    - ``0.75``

  * - ``orbital_period_max``
    - | Maximum orbital period (in days).
    - ``6000.0``

  * - ``#orbital_separation_scheme``
    - | Used only for ``orbital_scheme = 'separation'``. 
      | Options:
      
      * ``'log_uniform'``: log-uniform distribution
      * ``'log_normal'``: log-normal distribution
    - ``'log_uniform'``

  * - ``#orbital_separation_min``
    - | Minimum orbital separation (in solar radii).
    - ``5.0``

  * - ``#orbital_separation_max``
    - | Maximum orbital separation (in solar radii).
    - ``1e5``

  * - ``#log_orbital_separation_mean``
    - | Used only for ``orbital_separation_scheme = 'log_normal'``.
    - ``None``

  * - ``#log_orbital_separation_sigma``
    - | Used only for ``orbital_separation_scheme = 'log_normal'``.
    - ``None``

  * - ``eccentricity_scheme``
    - | How the initial binary eccentricity is sampled.
      | Options:
    
      * ``'zero'`` : zero eccentricity
      * ``'thermal'``: thermal distribution
      * ``'uniform'``: uniform distribution
    
    - ``'zero'``


Saving Output
-------------

You can decide on your own output parameters for the population file.
The data is split in two different tables: the ``history`` table and the ``oneline`` table.

The ``history`` table contains values that change throughout the evolution of the system, 
while the ``oneline`` table contains values that are constant throughout the evolution of the system or only occur once.


BinaryStar Output
~~~~~~~~~~~~~~~~~

The :class:`~posydon.binary_evol.binarystar` class contains the binary systems and
the parameters here determine what output of that class will be outputted into the final population file.
`scalar_names` will only be stored in the oneline table.

.. list-table::
  :widths: 10 80 10
  :class: population-params-table
  :header-rows: 1

  * - Parameter
    - Description
    - Default Value

  * - ``extra_columns``
    - | Additional columns to include in the output.
      | Note: ``'step_times'`` requires ``TimingHooks`` from ``posydon.binary_evol.simulationproperties``.
    - .. code:: python
      
        {'step_names':'string', 'step_times':'float64'}

  * - ``only_select_columns``
    - | Columns to include in the output.
      | Example: ``['state', 'event', 'time', 'orbital_period', 'eccentricity', 'lg_mtransfer_rate']``
      | Options:

      * ``'state'``: The state of the binary system.
      * ``'event'``: The event occurring in the binary system.
      * ``'time'``: The time of the event in years.
      * ``'separation'``: The separation between the binary stars in solar radii.
      * ``'orbital_period'``: The orbital period of the binary system in days.
      * ``'eccentricity'``: The eccentricity of the binary orbit.
      * ``'V_sys'``: The systemic velocity of the binary system.
      * ``'rl_relative_overflow_1'``: The Roche lobe relative overflow for star 1.
      * ``'rl_relative_overflow_2'``: The Roche lobe relative overflow for star 2.
      * ``'lg_mtransfer_rate'``: The logarithm of the mass transfer rate.
      * ``'mass_transfer_case'``: The case of mass transfer.
      * ``'trap_radius'``: The trapping radius.
      * ``'acc_radius'``: The accretion radius.
      * ``'t_sync_rad_1'``: The synchronization timescale for radiative zones star 1.
      * ``'t_sync_conv_1'``: The synchronization timescale for convective zones star 1.
      * ``'t_sync_rad_2'``: The synchronization timescale for radiative zones star 2.
      * ``'t_sync_conv_2'``: The synchronization timescale for convective zones star 2.
      * ``'nearest_neighbour_distance'``: The distance to the nearest neighbour.
    - .. code:: python
      
        ['state',
        'event',
        'time',
        'orbital_period',
        'eccentricity',
        'lg_mtransfer_rate']

  * - ``scalar_names``
    - | Scalars of the binary to include in the output.
      | Options:
      
      * ``'interp_class_HMS_HMS'``: interpolation class for HMS-HMS
      * ``'interp_class_CO_HMS_RLO'``: interpolation class for CO-HMS-RLO
      * ``'interp_class_CO_HeMS'``: interpolation class for CO-HeMS
      * ``'interp_class_CO_HeMS_RLO'``: interpolation class for CO-HeMS-RLO
      * ``'mt_history_HMS_HMS'``: mass transfer history for HMS-HMS
      * ``'mt_history_CO_HMS_RLO'``: mass transfer history for CO-HMS-RLO
      * ``'mt_history_CO_HeMS'``: mass transfer history for CO-HeMS
      * ``'mt_history_CO_HeMS_RLO'``: mass transfer history for CO-HeMS-RLO

    - .. code:: python

        ['interp_class_HMS_HMS',
        'interp_class_CO_HMS_RLO',
        'interp_class_CO_HeMS',
        'interp_class_CO_HeMS_RLO',
        'mt_history_HMS_HMS',
        'mt_history_CO_HMS_RLO',
        'mt_history_CO_HeMS',
        'mt_history_CO_HeMS_RLO']
    
              
SingleStar 1 and 2 Output
~~~~~~~~~~~~~~~~~~~~~~~~~

This dictionary contains the parameters that will be saved in the output of the SingleStar objects in the system.
`only_select_columns` will be stored in the history table,  and the initial and final step will be stored in the oneline table with the prefix :code:`S1` or :code:`S2` depending on the star,
`scalar_names` will only be stored in the oneline table.

.. list-table::
  :widths: 10 80 10
  :class: population-params-table
  :header-rows: 1

  * - Parameter
    - Description
    - Default Value

  * - ``include_S1``
    - | Include SingleStar 1 in the output.
    - ``True``

  * - ``only_select_columns``
    - | Columns to include in the output for SingleStar 1.
      | Example: ``['state', 'mass', 'log_R', 'log_L', 'lg_mdot', 'he_core_mass', 'he_core_radius', 'co_core_mass', 'co_core_radius', 'center_h1', 'center_he4', 'surface_h1', 'surface_he4', 'surf_avg_omega_div_omega_crit', 'spin']``
      | Options:
      
      * ``'state'``: The state of the star.
      * ``'metallicity'``: The metallicity of the star.
      * ``'mass'``: The mass of the star.
      * ``'log_R'``: The logarithm of the radius of the star.
      * ``'log_L'``: The logarithm of the luminosity of the star.
      * ``'lg_mdot'``: The logarithm of the mass loss rate.
      * ``'lg_system_mdot'``: The logarithm of the system mass loss rate.
      * ``'lg_wind_mdot'``: The logarithm of the wind mass loss rate.
      * ``'he_core_mass'``: The mass of the helium core.
      * ``'he_core_radius'``: The radius of the helium core.
      * ``'c_core_mass'``: The mass of the carbon core.
      * ``'c_core_radius'``: The radius of the carbon core.
      * ``'o_core_mass'``: The mass of the oxygen core.
      * ``'o_core_radius'``: The radius of the oxygen core.
      * ``'co_core_mass'``: The mass of the carbon-oxygen core.
      * ``'co_core_radius'``: The radius of the carbon-oxygen core.
      * ``'center_h1'``: The hydrogen fraction at the center.
      * ``'center_he4'``: The helium fraction at the center.
      * ``'center_c12'``: The carbon-12 fraction at the center.
      * ``'center_n14'``: The nitrogen-14 fraction at the center.
      * ``'center_o16'``: The oxygen-16 fraction at the center.
      * ``'surface_h1'``: The hydrogen fraction at the surface.
      * ``'surface_he4'``: The helium fraction at the surface.
      * ``'surface_c12'``: The carbon-12 fraction at the surface.
      * ``'surface_n14'``: The nitrogen-14 fraction at the surface.
      * ``'surface_o16'``: The oxygen-16 fraction at the surface.
      * ``'log_LH'``: The logarithm of the hydrogen burning luminosity.
      * ``'log_LHe'``: The logarithm of the helium burning luminosity.
      * ``'log_LZ'``: The logarithm of the metallicity burning luminosity.
      * ``'log_Lnuc'``: The logarithm of the nuclear burning luminosity.
      * ``'c12_c12'``: The carbon-12 to carbon-12 reaction rate.
      * ``'center_gamma'``: The central gamma value.
      * ``'avg_c_in_c_core'``: The average carbon fraction in the carbon core.
      * ``'surf_avg_omega'``: The average surface angular velocity.
      * ``'surf_avg_omega_div_omega_crit'``: The ratio of the average surface angular velocity to the critical angular velocity.
      * ``'total_moment_of_inertia'``: The total moment of inertia.
      * ``'log_total_angular_momentum'``: The logarithm of the total angular momentum.
      * ``'spin'``: The spin of the star.
      * ``'conv_env_top_mass'``: The mass at the top of the convective envelope.
      * ``'conv_env_bot_mass'``: The mass at the bottom of the convective envelope.
      * ``'conv_env_top_radius'``: The radius at the top of the convective envelope.
      * ``'conv_env_bot_radius'``: The radius at the bottom of the convective envelope.
      * ``'conv_env_turnover_time_g'``: The turnover time of the convective envelope (gravity).
      * ``'conv_env_turnover_time_l_b'``: The turnover time of the convective envelope (local buoyancy).
      * ``'conv_env_turnover_time_l_t'``: The turnover time of the convective envelope (local turbulence).
      * ``'envelope_binding_energy'``: The binding energy of the envelope.
      * ``'mass_conv_reg_fortides'``: The mass of the convective region for tides.
      * ``'thickness_conv_reg_fortides'``: The thickness of the convective region for tides.
      * ``'radius_conv_reg_fortides'``: The radius of the convective region for tides.
      * ``'lambda_CE_1cent'``: The lambda parameter for common envelope evolution (1%).
      * ``'lambda_CE_10cent'``: The lambda parameter for common envelope evolution (10%).
      * ``'lambda_CE_30cent'``: The lambda parameter for common envelope evolution (30%).
      * ``'lambda_CE_pure_He_star_10cent'``: The lambda parameter for common envelope evolution for pure helium stars (10%).
      * ``'profile [not currently supported]'``: The profile (not currently supported).
    - .. code:: python

        ['state',
        'mass',
        'log_R',
        'log_L',
        'lg_mdot',
        'he_core_mass',
        'he_core_radius',
        'co_core_mass',
        'co_core_radius',
        'center_h1',
        'center_he4',
        'surface_h1',
        'surface_he4',
        'surf_avg_omega_div_omega_crit',
        'spin',]


  * - ``scalar_names``
    - | Scalars to include in the output for SingleStar 1.
      | Example: ``['natal_kick_array', 'SN_type', 'f_fb', 'spin_orbit_tilt_first_SN', 'spin_orbit_tilt_second_SN', 'm_disk_accreted', 'm_disk_radiated']``
      | Options:
      
      * ``'natal_kick_array'``: The array of natal kicks.
      * ``'SN_type'``: The type of supernova.
      * ``'f_fb'``: The fallback fraction.
      * ``'spin_orbit_tilt_first_SN'``: The spin-orbit tilt after the first supernova.
      * ``'spin_orbit_tilt_second_SN'``: The spin-orbit tilt after the second supernova.
      * ``'m_disk_accreted'``: The mass accreted onto the disk.
      * ``'m_disk_radiated'``: The mass radiated from the disk.
    - .. code:: python

        ['natal_kick_array',
        'SN_type',
        'f_fb',
        'spin_orbit_tilt_first_SN',
        'spin_orbit_tilt_second_SN',]



