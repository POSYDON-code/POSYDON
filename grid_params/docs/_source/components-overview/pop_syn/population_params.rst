.. _pop-params-guide:

================================================
POSDYON Population Synthesis Configuration Guide
================================================

This documentation provides a detailed overview of the configuration options available in the Posydon software package.

TODO: fill each table cell with a description of the parameter and the options

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


--------------------
SimulationProperties
--------------------

After the environment variables, the next part of the ``population_params.ini`` file
contains how systems will move through the simulation steps (Flow Chart) and the parameters for each step.

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
   :widths: 50 50
   :header-rows: 1

   * - Parameter
     - Description/Value
   * - `import`
     - ['posydon.binary_evol.flow_chart', 'flow_chart']
   * - absolute_import = None
     - 'package' (kwarg for importlib.import_module)

Step MESA (HMS-HMS, CO-HMS_RLO, CO-HeMS, CO-HeMS_RLO)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The MESA step is the most important step of POSYDON as it leverages the POSYDON MESA grids to evolve the binary object according to one of the supported MESA binary-star grids.

.. list-table::
   :widths: 50 50
   :header-rows: 1

   * - Parameter
     - Description/Value
   * - `import`
     - ['posydon.binary_evol.MESA.step_mesa', 'MS_MS_step']
   * - absolute_import = None
     - 'package' (kwarg for importlib.import_module)
   * - `interpolation_path`
     - None (found by default)
   * - `interpolation_filename`
     - None (found by default)
   * - `interpolation_method`
     - '1NN_1NN' ('nearest_neighbour', 'linear3c_kNN', '1NN_1NN' are options)
   * - `save_initial_conditions`
     - True (only for interpolation_method='nearest_neighbour')
   * - `track_interpolation`
     - False
   * - `stop_method`
     - 'stop_at_max_time' ('stop_at_end', 'stop_at_max_time', 'stop_at_condition' are options)
   * - `stop_star`
     - 'star_1' (only for stop_method='stop_at_condition', 'star_1' and 'star_2' are options)
   * - `stop_var_name`
     - None (only for stop_method='stop_at_condition', string)
   * - `stop_value`
     - None (only for stop_method='stop_at_condition', float)
   * - `stop_interpolate`
     - True
   * - `verbose`
     - False


Step Detached
~~~~~~~~~~~~~

The detached step uses analytical expressions and MESA single star grids to evolve the binary object when it is not interacting.


.. list-table::
   :widths: 50 50
   :header-rows: 1

   * - Parameter
     - Description/Value
   * - `import`
     - ['posydon.binary_evol.DT.step_detached', 'detached_step']
   * - `absolute_import`
     - None ('package' kwarg for importlib.import_module)
   * - `matching_method`
     - 'minimize' (options 'minimize' 'root')
   * - `do_wind_loss`
     - True
   * - `do_tides`
     - True
   * - `do_gravitational_radiation`
     - True
   * - `do_magnetic_braking`
     - True
   * - `do_stellar_evolution_and_spin_from_winds`
     - True
   * - `RLO_orbit_at_orbit_with_same_am`
     - False
   * - `verbose`
     - False
    
Step Disrupted
~~~~~~~~~~~~~~

The dirtupted step evolves a system when a supernova has unbound the binary components.
This class inherits from the detached step, but only evolves the remaining star in isolation.

.. list-table::
   :widths: 50 50
   :header-rows: 1

   * - Parameter
     - Description/Value
   * - `import`
     - ['posydon.binary_evol.DT.step_disrupted','DisruptedStep']

Step Merged
~~~~~~~~~~~

In this step, the system has undergone a merger and is now a single star.
This class inherits from the detached step, but only evolves the remaining star in isolation.


.. list-table::
   :widths: 50 50
   :header-rows: 1

   * - Parameter
     - Description/Value
   * - `import`
     - ['posydon.binary_evol.DT.step_merged','MergedStep']

Step Initially Single
~~~~~~~~~~~~~~~~~~~~~

This step is used to evolve a single star system.
This class inherits from the detached step and evolves the star in isolation.

.. list-table::
   :widths: 50 50
   :header-rows: 1

   * - Parameter
     - Description/Value
   * - `import`
     - ['posydon.binary_evol.DT.step_initially_single','InitiallySingleStep']

Step Common Envelope
~~~~~~~~~~~~~~~~~~~~

The common envelope step is used to evolve a binary system when the primary or secondary star has initiated a common envelope phase.
It calculates the binding energy of the envelope and the energy available to eject it.
If the energy budget is greater than the binding energy, the envelope is ejected and the system is evolved to the next step.
If the energy budget is less than the binding energy, the system is merged.

.. list-table::
   :widths: 50 50
   :header-rows: 1

   * - Parameter
     - Description/Value
   * - `import`
     - ['posydon.binary_evol.CE.step_CEE', 'StepCEE']
   * - `absolute_import`
     - None('package' kwarg for importlib.import_module)
   * - `prescription`
     - 'alpha-lambda'
   * - `common_envelope_efficiency`
     - 1.0 (float in [0, inf])
   * - `common_envelope_option_for_lambda`
     - 'lambda_from_grid_final_values' (options are: (1) 'default_lambda', (2) 'lambda_from_grid_final_values', (3) 'lambda_from_profile_gravitational', (4) 'lambda_from_profile_gravitational_plus_internal', (5) 'lambda_from_profile_gravitational_plus_internal_minus_recombination')
   * - `common_envelope_lambda_default`
     - 0.5 (float in [0, inf] used only for option (1))
   * - `common_envelope_option_for_HG_star`
     - 'optimistic' (options are 'optimistic', 'pessimistic')
   * - `common_envelope_alpha_thermal`
     - 1.0 (float in [0, inf] used only for option for (4), (5))
   * - `core_definition_H_fraction`
     - 0.1 (options are 0.01, 0.1, 0.3)
   * - `core_definition_He_fraction`
     - 0.1
   * - `CEE_tolerance_err`
     - 0.001 (float in [0, inf])
   * - `common_envelope_option_after_succ_CEE`
     - 'core_not_replaced_noMT' (options are 'core_not_replaced_noMT' 'core_replaced_noMT', 'core_not_replaced_stableMT' 'core_not_replaced_windloss')
   * - `verbose`
     - False

Step Supernova
~~~~~~~~~~~~~~

This steps performs the core-collapse supernova evolution of the star.
Multiple prescriptions are implemented and can be selected.
If a standard prescription is used, interpolators are available to speed up the calculation and predict additional properties of the collapse, such as the spin.


.. list-table::
   :widths: 50 50
   :header-rows: 1

   * - Parameter
     - Description/Value
   * - `import`
     - ['posydon.binary_evol.SN.step_SN', 'StepSN']
   * - `absolute_import`
     - None ('package' kwarg for importlib.import_module)
   * - `mechanism`
     - 'Patton&Sukhbold20-engine' (options are: 'direct', Fryer+12-rapid', 'Fryer+12-delayed', 'Sukhbold+16-engine', 'Patton&Sukhbold20-engine')
   * - `engine`
     - 'N20' (options are 'N20' for 'Sukhbold+16-engine', 'Patton&Sukhbold20-engine' or None for the others)
   * - `PISN`
     - 'Marchant+19' (options are None, "Marchant+19")
   * - `ECSN`
     - "Podsiadlowksi+04" (options are "Tauris+15", "Podsiadlowksi+04")
   * - `conserve_hydrogen_envelope`
     - True
   * - `max_neutrino_mass_loss`
     - 0.5 (float in [0,inf])
   * - `max_NS_mass`
     - 2.5 (float in [0,inf])
   * - `use_interp_values`
     - True (if True, use interpolation values for the SN properties, which are calculated during ``step_MESA``)
   * - `use_profiles`
     - True
   * - `use_core_masses`
     - True
   * - `approx_at_he_depletion`
     - False
   * - `kick`
     - True
   * - `kick_normalisation`
     - 'one_over_mass' (options are "one_minus_fallback", "one_over_mass", "NS_one_minus_fallback_BH_one", "one", "zero")
   * - `sigma_kick_CCSN_NS`
     - 265.0 (float in [0,inf])
   * - `sigma_kick_CCSN_BH`
     - 265.0 (float in [0,inf])
   * - `sigma_kick_ECSN`
     - 20.0 (float in [0,inf])
   * - `verbose`
     - False

Step Double Compact Object
~~~~~~~~~~~~~~~~~~~~~~~~~~

In this step, the system has evolved to a double compact object system.
The merger time due to gravitational wave emission is calculated.

.. list-table::
   :widths: 50 50
   :header-rows: 1

   * - Parameter
     - Description/Value
   * - `import`
     - ['posydon.binary_evol.DT.double_CO', 'DoubleCO']
   * - `absolute_import`
     - None ('package' kwarg for importlib.import_module)
   * - `n_o_steps_interval`
     - None

Step End
~~~~~~~~

The final step of the simulation.
All succesfull systems should reach this step.

.. list-table::
   :widths: 50 50
   :header-rows: 1

   * - Parameter
     - Description/Value
   * - `import`
     - ['posydon.binary_evol.step_end', 'step_end']
   * - `absolute_import`
     - None ('package' kwarg for importlib.import_module)

Extra Hooks
~~~~~~~~~~~

It's possible to add custom hooks to the simulation steps.
A few example hooks are provides: ``TimingHooks`` and ``StepNamesHooks`` (See :ref:`Custom Hooks <custom-hooks>` for more details)

.. list-table::
   :widths: 50 50
   :header-rows: 1

   * - Parameter
     - Description/Value
   * - `import`
     - ['posydon.binary_evol.simulationproperties', 'TimingHooks']
   * - `absolute_import_1`
     - None
   * - `kwargs_1`
     - {}
   * - `import`
     - ['posydon.binary_evol.simulationproperties', 'StepNamesHooks']
   * - `absolute_import_2`
     - None
   * - `kwargs_2`
     - {}



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
   :widths: 50 50
   :header-rows: 1

   * - Parameter
     - Description/Value
   * - `optimize_ram`
     - True (save population in batches)
   * - `ram_per_cpu`
     - None (set maximum ram per cpu before batch saving in GB)
   * - `dump_rate`
     - 1000 (batch save after evolving N binaries)
   * - `temp_directory`
     - 'batches' (folder for keeping batch files)
   * - `tqdm`
     - False (progress bar)
   * - `breakdown_to_df`
     - True (convert BinaryStars into DataFrames after evolution)
   * - `use_MPI`
     - True ( if True evolve with MPI, equivalent to the following: from mpi4py import MPI, comm = MPI.COMM_WORLD)
   * - `metallicity`
     - [2., 1., 0.45, 0.2, 0.1, 0.01, 0.001, 0.0001] (In units of solar metallicity)
   * - `entropy`
     - `None` (Random Number Generation: uses system entropy (recommended))
   * - `number_of_binaries`
     - 1000000 (int)
   * - `star_formation`
     - 'burst' (options are 'constant' 'burst' 'custom_linear' 'custom_log10' 'custom_linear_histogram' 'custom_log10_histogram')
   * - `max_simulation_time`
     - 13.8e9 (float in [0,inf])
   * - `binary_fraction`
     - 1 (float 0< fraction <=1)
   * - `primary_mass_scheme`
     - 'Kroupa2001' (options are 'Salpeter', 'Kroupa1993', 'Kroupa2001')
   * - `primary_mass_min`
     - 6.5 (float in [0,300])
   * - `primary_mass_max`
     - 250.0 (float in [0,300])
   * - `secondary_mass_scheme`
     - 'flat_mass_ratio' (options are 'flat_mass_ratio', 'q=1')
   * - `secondary_mass_min`
     - 0.35 (float in [0,300])
   * - `secondary_mass_max`
     - 250.0 (float in [0,300])
   * - `orbital_scheme``
     - 'period' (options are 'separation', 'period')
   * - `orbital_period_scheme`
     - 'Sana+12_period_extended' (used only for orbital_scheme = 'period')
   * - `orbital_period_min`
     - 0.75 (float i [0,inf])
   * - `orbital_period_max`
     - 6000.0 (float i [0,inf])
   * - `#orbital_separation_scheme`
     - 'log_uniform' (used only for orbital_scheme = 'separation', 'log_uniform', 'log_normal')
   * - `#orbital_separation_min`
     - 5.0 (float i [0,inf])
   * - `#orbital_separation_max`
     - 1e5 (float i [0,inf])
   * - `#log_orbital_separation_mean`
     - None (float i [0,inf] used only for orbital_separation_scheme ='log_normal')
   * - `#log_orbital_separation_sigma`
     - None (float i [0,inf] used only for orbital_separation_scheme ='log_normal')
   * - `eccentricity_scheme`
     - 'zero' (options are 'zero', 'thermal', 'uniform')


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

.. list-table::
   :widths: 50 50
   :header-rows: 1

   * - Parameter
     - Description/Value
   * - `extra_columns`
     - {'step_names':'string', 'step_times':'float64'} ('step_times' with from posydon.binary_evol.simulationproperties import TimingHooks)
   * - `only_select_columns`
     - ['state', 'event', 'time', 'orbital_period', 'eccentricity', 'lg_mtransfer_rate'] (all options: 'state', 'event', 'time', 'separation', 'orbital_period', 'eccentricity', 'V_sys', 'rl_relative_overflow_1', 'rl_relative_overflow_2', 'lg_mtransfer_rate', 'mass_transfer_case', 'trap_radius', 'acc_radius', 't_sync_rad_1', 't_sync_conv_1', 't_sync_rad_2', 't_sync_conv_2', 'nearest_neighbour_distance')
                      
              
SingleStar 1 and 2 Output
~~~~~~~~~~~~~~~~~~~~~~~~~

This dictionary contains the parameters that will be saved in the output of the SingleStar objects in the system.
`only_select_columns` will be stored in the history table,  and the initial and final step will be stored in the oneline table with the prefix :code:`S1` or :code:`S2` depending on the star,
`scalar_names` will only be stored in the oneline table.

.. list-table::
   :widths: 50 50
   :header-rows: 1

   * - Parameter
     - Description/Value
   * - `include_S1`
     - True
   * - `only_select_columns`
     - ['state', 'mass', 'log_R', 'log_L', 'lg_mdot', 'he_core_mass', 'he_core_radius', 'co_core_mass', 'co_core_radius', 'center_h1', 'center_he4', 'surface_h1', 'surface_he4', 'surf_avg_omega_div_omega_crit', 'spin',] (options are: 'state', 'metallicity', 'mass', 'log_R', 'log_L', 'lg_mdot', 'lg_system_mdot', 'lg_wind_mdot', 'he_core_mass', 'he_core_radius', 'c_core_mass', 'c_core_radius', 'o_core_mass', 'o_core_radius', 'co_core_mass', 'co_core_radius', 'center_h1', 'center_he4', 'center_c12', 'center_n14', 'center_o16', 'surface_h1', 'surface_he4', 'surface_c12', 'surface_n14', 'surface_o16', 'log_LH', 'log_LHe', 'log_LZ', 'log_Lnuc', 'c12_c12', 'center_gamma', 'avg_c_in_c_core', 'surf_avg_omega', 'surf_avg_omega_div_omega_crit', 'total_moment_of_inertia', 'log_total_angular_momentum', 'spin', 'conv_env_top_mass', 'conv_env_bot_mass', 'conv_env_top_radius', 'conv_env_bot_radius', 'conv_env_turnover_time_g', 'conv_env_turnover_time_l_b', 'conv_env_turnover_time_l_t', 'envelope_binding_energy', 'mass_conv_reg_fortides', 'thickness_conv_reg_fortides', 'radius_conv_reg_fortides', 'lambda_CE_1cent', 'lambda_CE_10cent', 'lambda_CE_30cent', 'lambda_CE_pure_He_star_10cent', 'profile [not currently supported]')           
   * - `scalar_names`
     - [ 'natal_kick_array', 'SN_type', 'f_fb', 'spin_orbit_tilt_first_SN','spin_orbit_tilt_second_SN', 'm_disk_accreted', 'm_disk_radiated']




