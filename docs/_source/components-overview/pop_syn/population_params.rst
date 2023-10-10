.. _pop-params:

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

SimulationProperties
--------------------

TODO: add description

Flow Chart
~~~~~~~~~~

The flow chart is the core of POSYDON. It controls the mapping between a POSYDON binary object and its step evolution, see the :ref:`Flow Chart Object <flow-chart>` page for more details.

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

TODO: add description

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

TODO: add description

.. list-table::
   :widths: 50 50
   :header-rows: 1

   * - Parameter
     - Description/Value
   * - `import`
     - ['posydon.binary_evol.DT.step_disrupted','DisruptedStep']

Step Merged
~~~~~~~~~~~

TODO: add description

.. list-table::
   :widths: 50 50
   :header-rows: 1

   * - Parameter
     - Description/Value
   * - `import`
     - ['posydon.binary_evol.DT.step_merged','MergedStep']

Step Initially Single
~~~~~~~~~~~~~~~~~~~~~

TODO: add description

.. list-table::
   :widths: 50 50
   :header-rows: 1

   * - Parameter
     - Description/Value
   * - `import`
     - ['posydon.binary_evol.DT.step_initially_single','InitiallySingleStep']

Step Common Envelope
~~~~~~~~~~~~~~~~~~~~

TODO: add description

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

TODO: add description

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
     - True
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

TODO: add description

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

TODO: add description

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

TODO: add description

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

TODO: add description

BinaryPopulation Options
~~~~~~~~~~~~~~~~~~~~~~~~

TODO: add description

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

TODO: add description

BinaryStar Output
~~~~~~~~~~~~~~~~~

TODO: add description

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

TODO: add description

.. list-table::
   :widths: 50 50
   :header-rows: 1

   * - Parameter
     - Description/Value
   * - `include_S1`
     - True
   * - `only_select_columns`
     - ['state', 'mass', 'log_R', 'log_L', 'lg_mdot', 'he_core_mass', 'he_core_radius', 'co_core_mass', 'co_core_radius', 'center_h1', 'center_he4', 'surface_h1', 'surface_he4', 'surf_avg_omega_div_omega_crit', 'spin',] (options are: 'state', 'metallicity', 'mass', 'log_R', 'log_L', 'lg_mdot', 'lg_system_mdot', 'lg_wind_mdot', 'he_core_mass', 'he_core_radius', 'c_core_mass', 'c_core_radius', 'o_core_mass', 'o_core_radius', 'co_core_mass', 'co_core_radius', 'center_h1', 'center_he4', 'center_c12', 'center_n14', 'center_o16', 'surface_h1', 'surface_he4', 'surface_c12', 'surface_n14', 'surface_o16', 'log_LH', 'log_LHe', 'log_LZ', 'log_Lnuc', 'c12_c12', 'center_gamma', 'avg_c_in_c_core', 'surf_avg_omega', 'surf_avg_omega_div_omega_crit', 'total_moment_of_inertia', 'log_total_angular_momentum', 'spin', 'conv_env_top_mass', 'conv_env_bot_mass', 'conv_env_top_radius', 'conv_env_bot_radius', 'conv_env_turnover_time_g', 'conv_env_turnover_time_l_b', 'conv_env_turnover_time_l_t', 'envelope_binding_energy', 'mass_conv_reg_fortides', 'thickness_conv_reg_fortides', 'radius_conv_reg_fortides', 'lambda_CE_1cent', 'lambda_CE_10cent', 'lambda_CE_30cent', 'lambda_CE_pure_He_star_10cent', 'profile')           
   * - `scalar_names`
     - [ 'natal_kick_array', 'SN_type', 'f_fb', 'spin_orbit_tilt', 'm_disk_accreted', 'm_disk_radiated']
