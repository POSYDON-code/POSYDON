.. _debugging_binaries:

#########################
Debugging Failed Binaries
#########################

During evolution in the ``BinaryPopulation``, all ``Exceptions`` are caught
to allow the population to keep evolving in the event a single binary
enocunters an error. Here we go over some methods to help debug errors.
We cause an error by using a broken evolutionary step:

.. code-block:: python
  :linenos:
  :emphasize-lines: 12-17,32

  from posydon.popsyn.binarypopulation import BinaryPopulation
  from posydon.binary_evol.simulationproperties import SimulationProperties
  from posydon.binary_evol.flow_chart import flow_chart
  from posydon.binary_evol.CE.step_CEE import StepCEE
  from posydon.binary_evol.SN.step_SN import StepSN
  from posydon.binary_evol.step_end import step_end
  from posydon.binary_evol.MESA.step_mesa import CO_HeMS_step, MS_MS_step, CO_HMS_RLO_step
  from posydon.binary_evol.DT.step_detached import detached_step
  from posydon.binary_evol.DT.double_CO import DoubleCO
  from posydon.binary_evol.simulationproperties import TimingHooks, StepNamesHooks

  class my_new_StepSN(StepSN):
      def __call__(self, binary):
          # use the StepSN call
          super().__call__(binary)
          # break here
          raise ValueError(f'{str(binary)} failed in my new step')

  mesa_step_kwargs = dict(
      interpolation_method='linear3c_kNN',
      save_initial_conditions=False,
      track_interpolation=False
  )

  sim_kwargs = dict(
      flow = (flow_chart, {}),
      step_HMS_HMS = (MS_MS_step, mesa_step_kwargs),
      step_CO_HeMS = (CO_HeMS_step, mesa_step_kwargs),
      step_CO_HMS_RLO = (CO_HMS_RLO_step, mesa_step_kwargs),
      step_detached = (detached_step, {}),
      step_CE = (StepCEE, dict(core_definition_H_fraction= 0.1,)),
      step_SN = (my_new_StepSN, {}),
      step_dco = (DoubleCO, {}),
      step_end = (step_end, {}),
      extra_hooks = [(TimingHooks, {}),(StepNamesHooks, {})]
  )

  kwargs = {'number_of_binaries' : 15,
            'primary_mass_min' : 7,
            'primary_mass_max' : 127,
            'secondary_mass_scheme' : 'flat_mass_ratio',
            'secondary_mass_min': 1,
            'secondary_mass_max': 127,
            'orbital_separation_min': 1,
            'orbital_separation_max': 3e3,
            'eccentricity_scheme':'zero',
            'extra_columns' : ['step_times','step_names'],
            'only_select_columns' : ['state', 'event', 'time', 'lg_mtransfer_rate',
                                     'orbital_period'],
            'include_S1' : True ,
            'S1_kwargs' : {'only_select_columns' : ['state', 'mass',
                                    'log_R', 'center_h1','center_he4',
                                    'he_core_mass', 'surface_he4', 'surface_h1',
                                    'lg_mdot'],
                          'scalar_names' : ['natal_kick_array', 'SN_type']},
            'include_S2' : True,
            'S2_kwargs' : {'only_select_columns' : ['state', 'mass',
                                    'log_R', 'center_h1','center_he4',
                                    'he_core_mass', 'surface_he4', 'surface_h1',
                                    'lg_mdot'],
                          'scalar_names' : ['natal_kick_array', 'SN_type']},
            'star_formation' : 'burst',
            'max_simulation_time' : 13.7e9,
           }

  pop = BinaryPopulation(
      population_properties=sim_prop,
      entropy = 12345678,
      **kwargs
  )

  pop.evolve(breakdown_to_df=False, from_hdf=False, tqdm=True, **kwargs)


After running the script above we can use the ``BinaryPopulation.find_failed``
instance method which returns a list of the failed binaries in the population.

.. code-block:: python

  failed_binaries = pop.find_failed()

  if failed_binaries:
      print( failed_binaries[0].traceback )

The ``BinaryStar`` instances that fail during evolution will have a
``traceback`` attribute showing the string representation (not the actual
stack
`traceback <https://docs.python.org/3/library/traceback.html>`_):

.. code-block:: text

  Traceback (most recent call last):
    File "/Users/kylerocha/Private/Github/POSYDON/posydon/popsyn/binarypopulation.py", line 153, in _safe_evolve
      binary.evolve()
    File "/Users/kylerocha/Private/Github/POSYDON/posydon/binary_evol/binarystar.py", line 183, in evolve
      self.run_step()
    File "/Users/kylerocha/Private/Github/POSYDON/posydon/binary_evol/binarystar.py", line 214, in run_step
      next_step(self)
    File "<ipython-input-1-3023b4d9482c>", line 32, in __call__
      raise ValueError(f'{str(binary)} failed in my new step')
  ValueError: BinaryStar(detached, None, p=6820.50, S1=(BH,M=26.79), S2=(H-rich_Shell_H_burning,M=35.61)) failed in my new step


If we want to go further we could also take the binary out of the population
and evolve it. If using a ``jupyter notebook``, you can enter the python
debugger by raising the error again and using the magic command ``%debug``
in another cell.


.. code-block:: ipython3

  failing_binary = failed_binaries[0]

  failing_binary.restore()

  failing_binary.evolve()

-----

.. code-block:: python

  %debug

See the debugger `commands
<https://docs.python.org/3/library/pdb.html#debugger-commands>`_ for help.
