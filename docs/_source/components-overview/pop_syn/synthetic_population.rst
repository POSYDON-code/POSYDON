.. _synthetic-population:

.. py:currentmodule:: posydon.popsyn.synthetic_population

Synthetic Populations
#####################


Population
==========

The :class:`~Population` object is an interface for the population data, which is stored in a HDF5 file.
You can find more information about the file structure here: :ref:`population-file-structure`.

The population file is created by a population synthesis run, and contains the history of each system in the population, and an oneline description of each system and its evolution.

The :class:`~Population` object is created by passing the path to the population file to the constructor.
It will not load in the data immediately, but only when requested.

.. code-block:: python

    from posydon.popsyn import Population

    pop = Population('path/to/population.h5')



There are three main components of the :class:`~Population` object:

- :class:`~Population.history`: the history of the population
- :class:`~Population.oneline`: the oneline data of the population
- :class:`~Population.formation_channels`: the formation channels of the population

The :code:`formation_channels` is not immediately available, but can be created using the
:func:`~Population.calculate_formation_channels` function of the :class:`~Population` object.


Additionally, :meth:`Population.mass_per_metallicity<Population.mass_per_metallicity>` contains some essential metadata calculated from the population synthesis run.
It is a pandas.DataFrame with the metallicity (in solar units) as the index and 3 columns:

1. :code:`count`, the number of systems at that metallicity in the file.
2. :code:`simulated_mass`, the total ZAMS mass of the population run at that metallicity.
3. :code:`underlying_mass`, the total mass of the population run at that metallicity, assuming a binary fraction of 1. This is calculated in :func:`~posydon.popsyn.normalized_pop_mass.initial_total_underlying_mass`.



history
--------

:class:`~Population.history` contains the full history of each system in the population.
The history is stored in a pandas DataFrame, where each row is a timestep in the evolution of a system.
The columns of the DataFrame are the parameters of the system at each timestep.

The history contains the data from the binary and the individual stars.
See :ref:`Single Star<single-star>` and :ref:`Binary Star<binary-star>` for more information about the columns in the history table.

The single star data is stored with the prefix :code:`S1` or :code:`S2` based on which star.

Accessing columns or binaries in the history table is easy.

.. code-block:: python

    binary_0 = pop.history[0]

This will return the complete history of the first binary in the population, including all its columns.

If you want to access a specific column, you can do so by using the column name.

.. code-block:: python

    mass_1 = pop.history['S1_mass']


A more powerful feature is the :func:`History.select` function, which allows you to select specific rows or columns from the history table. Here's an example:

.. note::
    The :func:`History.select` function only allows the use of :code:`where=`
    for specific columns and the index. The columns are limited to those containing strings.

.. code-block:: python

    # using where with the index
    mass_10 = pop.history.select(columns=['S1_mass'], where='index==10')

    # using where with a string column
    mass_ZAMS = pop.history.select(columns=['S1_mass'], where='event == "ZAMS"')
    
If you want to have a peak, you can use the :meth:`~Population.head` or :meth:`~Population.tail` functions.

.. code-block:: python
    pop.history.head(10)
    pop.history.tail(10)


Additional functions are made available for ease of use.

If you want to check the length of the history of a system, you can use :attr:`Population.history.lengths<History.lengths>` or :attr:`Population.history_lengths<Population.history_lengths>`.

.. code-block:: python

    print(pop.history.lengths)
    print(pop.history_lengths)

The total number of systems in the population can be found with :attr:`~Population.History.number_of_systems`.

.. code-block:: python

    print(pop.history.number_of_systems)

Similarly, if you would like to check the indices in the file, you can use :attr:`Population.indices<Population.indices>` or :attr:`Population.history.indices<History.indices>`.
The indices are useful in selecting systems from the population.

It's also possible to check the columns in the history table with :attr:`Population.columns<Population.columns>` or :attr:`Population.history.columns<History.columns>`.

.. code-block:: python

    print(pop.indices)
    print(pop.history.indices)

    print(pop.history.columns)
    print(pop.columns['history'])


oneline
--------

:meth:`~Population.oneline` contains a single line description of each system in the population.
This is useful for a quick inspection of the population.
The oneline data is stored in a pandas DataFrame, where each row is a system in the population.

Some properties over the evolution of the binary do not change, such as the natal kick properties or interpolation class.
Besides the initial and final properties of the system, this table also contains these data.

The initial-final properties are those in the history table, but with the postfix :code:`_i` and :code:`_f` depending on the initial or final value.
The additional values are the scalar values from the individual stars and the binary properties (See :ref:`Single Star<single-star>` and :ref:`Binary Star<binary-star>`).

Additionally, ``WARNING``, ``FAILED``, and metallicity columns are available in the oneline table.

.. csv-table:: Additional columns
  :header: "Properties", "Descriptions"
  :widths: 50, 150
  ``FAILED``, Indicates if the system failed during the population synthesis run.
  ``WARNING``, Indicates if there were any warnings for the system during the population synthesis run.
  ``metallicity``, The metallicity of the system.


Like the :code:`history` access, you can access the oneline data by using the index of the system or the columns.

.. code-block:: python

    binary_0 = pop.oneline[0]
    mass = pop.oneline['S1_mass_i']
    selection = pop.oneline.select(columns=['S1_mass_i'], where='index==10')

You can check the columns and indices of the oneline table with :attr:`Population.oneline.columns<Oneline.columns>` and :attr:`Population.columns['oneline']<Population.columns>`.

.. code-block:: python

    print(pop.oneline.columns)
    print(pop.columns['oneline'])

The number of systems in the population can be found with :attr:`Population.oneline.number_of_systems<Oneline.number_of_systems>`.
The length and indices of the oneline table can be found with :attr:`Population.oneline.lengths<Oneline.lengths>`, and :attr:`Population.oneline.indices<Oneline.indices>`, respectively.


.. code-block:: python

    print(pop.oneline.number_of_systems)
    print(pop.oneline.lengths)
    print(pop.oneline.indices)
  

formation_channels
------------------

:class:`~Population.formation_channels` contains the formation channels of each system in the population.
The formation channels are stored in a pandas DataFrame, where each row is a system in the population.

The formation channels are calculated by combining the `event` column in the history table into a single string using the :func:`~Population.calculate_formation_channels` function of the :class:`~Population` object.

Two columns are available in the formation channels table:

- `debug_channel` : A longer description of the formation channel, where additional events are included.

- `channel` : A cleaned-up version of the history events, where events are separated by a `-`. 


Additional Attributes
---------------------


- :attr:`~Population.ini_parameters`: The parameters for the initial sampling conditions of the population synthesis run.


- :attr:`~Population.mass_per_metallicity`: The mass per metallicity bin for the population synthesis run. 
  The `underlying_mass` is calculated with the assumption that binary fraction == 1.

- :attr:`~Population.history_lengths`: The length of the history of each system in the population. 
        This is created the first time the file is opened with the :class:`~Population` object.



Exporting part of the population
--------------------------------

The class function :func:`Population.export_selection<Population.export_selection>` allows you to export part of the population to a new HDF5 file.
It takes a list of indices of the systems you want to export, and the path to the new file.
This will copy the systems with the given indices to the new file, which includes their history, oneline data, and formation channels (if presen).

.. code-block:: python

    pop.export_selection([0, 1, 2], 'path/to/new_population.h5')


If the file already exists, the function will raise an error. If you want to overwrite or append to the file, you can use the :code:`overwrite` argument or :code:`append` argument.

.. code-block:: python

    pop.export_selection([0, 1, 2], 'path/to/new_population.h5', overwrite=True)
    pop.export_selection([0, 1, 2], 'path/to/new_population.h5', append=True)


TransientPopulation
===================

The :class:`~TransientPopulation` object is an interface for the transient data of the population.
It inherits from the :class:`~Population` object, but also allows access to the transient populations in the population file.

A transient population consists of instantaneous events in the population, such as supernovae, kilonovae, or gamma-ray bursts.
These have a single "moment" in time, and are not part of the evolution of the system.
The :class:`~TransientPopulation` class has been designed to handle these events, but future versions of the code may include more complex populations.

.. code-block:: python

    from posydon.popsyn.synthetic_population import TransientPopulation

    # where transient_name is the name of the transient population in the file
    trans_pop = TransientPopulation('path/to/population.h5', 'transient_name')


Creating a TransientPopulation
------------------------------

The transient population is created using the :func:`Population.create_transient_population` function of the :class:`~Population` object.
This function creates a separate table with each transient in the population file.
It loops over all the systems in the population in chunks and applies the given function to them.

The :func:`Population.create_transient_population` function takes a function as an argument: :code:`selection_function`.

The :code:`selection_function` takes 3 arguments: :code:`history_chunk`, :code:`oneline_chunk`, and :code:`formation_channels_chunk` (optional).
These chunks are cut based on a given chunksize, which is set to 1000000 by default, and are cut on system. 
This means that always a complete history of a system is passed to the function by :func:`Population.create_transient_population`.

:code:`selection_function` is a function you can adapt to your own needs, and
examples of building one are given in the `BBH <../../tutorials-examples/population-synthesis/bbh_analysis.html>`_ or `LGRB tutorial <../../tutorials-examples/population-synthesis/lgrb_pop_syn.html>`_.


.. note::

    The :code:`selection_function` should return a pandas DataFrame with at least a :code:`metallicity` and :code:`time` column of each event.
    Moreover, every row should only contain a single event.


We provide a few standard selection functions for the most common transient populations, such as binary black holes and long gamma-ray bursts.
These functions are available in the :mod:`posydon.popsyn.transient_select_funcs` module.


Accessing TransientPopulation
-----------------------------

After loading a transient population, you keep access to the history and oneline data of the population.
Now, you can access the transient data of the population using :attr:`TrannsientPopulation.population<TransientPopulation>`.


.. code-block:: python
  
    print(trans_pop.population)


Calculating Efficiencies
------------------------

With this population, you can calculate additional information, such as the efficiency over metallicity.

.. code-block:: python
  
    trans_pop.calculate_efficiency_over_metallicity(channels=True)

:code:`channels=True` includes the formation channels in the efficiency calculation.

Plotting
--------

The :class:`~TransientPopulation` contains a few plotting functions for ease.
    
.. code-block:: python
    # plots the efficiency over metallicity per channel
    trans_pop.plot_efficiency_over_metallicity(channel=True)


.. code-block:: python

    log_bins = np.logspace(6.5, 10.5, 51)
    trans_pop.plot_delay_time_distribution(metallicity=0.1, bins=log_bins)


The most useful function is :func:`plot_popsyn_over_grid_slice<TransientPopulation.plot_popsyn_over_grid_slice>`.
It allows you to overplot properties of your TransienPopulation onto the grids.

.. note::

  Make sure that you've set your :code:`PATH_TO_POSYDON` and :code:`PATH_TO_POSYDON_DATA` environment variables correctly.
  These are required to plot over the grid slices.


If you like to write to a folder, you can use :code:`plot_dir='path/to/dir'` and use :code:`save_fig=True`.

.. code-block:: python
  
    # plot the HMS-HMS grid at 1e-4 with S1_spin and q=0.7
    plot_popsyn_over_grid_slice('HMS-HMS', 1e-4, slices=[0.7], prop='S1_spin', prop_range=[0,0.3], save_fig=False, channel='ZAMS_oRLO1_CC1_oRLO2_CC2_END')


Rates
=====

The :class:`~Rates` object inherits from the :class:`~TransientPopulation` object 
and is used to access the cosmic rate data of the transient population.

It also allows the user to calculate the intrinsic rate density of the events in the population, and apply observational effects to the population.

.. code-block:: python

    from posydon.popsyn.synthetic_population import Rates

    rates = Rates('path/to/population.h5', 'transient_name', 'SFH_identifier')


Creating a Rates object
-----------------------

Cosmic weights are added to the population file using the :func:`~TransientPopulation.calculate_cosmic_weights` function.
This function calculates the cosmic weights of the events in the population based on the birth redshifts and the population weight.
The function takes an ``SFH_identifier``, which is where the cosmic weights are stored in the population file.
The ``MODEL_in`` argument is used to specify the model parameters for the rate calculation.

The table below shows the Default values and the supported values.

.. csv-table:: MODEL_in
  :header: "Parameter", "Value", "Description"
  :widths: 30, 30, 150

  "delta_t", 100, "The time interval to split the birth times into"
  "SFR", "IllustrisTNG", 'The star formation history identifier [IllustrisTNG/Madau+Fragos17/Madau+Dickinson14/Neijssel+19]""
  "sigma_SFR", None, "The uncertainty in the SFR (float) or the identifier of the SFR uncertainty [Bavera+20/Neijssel+19](str)"
  "Z_max", 1.0, "The maximum metallicity to consider [0.0-1.0]"
  "select_one_met", False, "Select one metallicity. Requires only one metallicity in the population"
  "dlogZ", None, "The metallicity bin width when selecting a single bin (float), the bin edges (tuple(float, float)), or if ``None`` with select_one_met=True, then all metallicities is used."
  "Zsun", Zsun, "The solar metallicity"


Accessing rates data
-----------------------------

The cosmic rate data is stored in 3 different tables in the population file:

1. :code:`birth` : A table containing the birth redshifts and lookback times used in the rate calculation.
2. :code:`z_events` : The redshifts of the events in the population and the birth redshifts of the events.
3. :code:`weights` : The weights of each event based on their birth redshifts and their population weight.


You can calculate the intrinsic rate density of the events in the population using :func:`Rates.calculate_intrinsic_rate_density`.
This populates the :code:`intrinsic_rate_density` table in the population file.

.. code-block:: python

    # calculate the intrinsic rate density per formation channel
    rates.calculate_intrinsic_rate_density(mt_channels=True)

    rates.plot_intrinsic_rate_density()

The :class:`~Rates` object also contains information about the metallicity and redshift bins and edges.

.. code-block:: python

    print(rates.centers_metallicity_bins)
    print(rates.edges_metallicity_bins)
    print(rates.centers_redshift_bins)
    print(rates.edges_redshift_bins)

Plotting Rates
-----------------------------

Besides plotting the intrinsic rate, you can plot the distribution of properties of the population.
You can use any property in the TransientPopulation table.

.. code-block:: python

    rates.plot_hist_properties('S1_mass', intrinsice=True, label='S1', show=True)


Applying observational effects
------------------------------

Although the intrinsic rate density is a useful quantity, it is not directly observable, especially for binary black holes.

As such, we also include the possibility to apply observational effects to the population.
This is done using the :func:`Rates.calculate_observable_population` function.
It reweights the event weights based on the detection efficiency of the event.
The function takes a function as an argument: :code:`observable_func` and a ``observable_identifier``.

The ``observable_func`` you give the function should take 3 arguments:
1. transient_chunk : The transient data of the population.
2. z_events_chunk : The redshifts of the events in the population.
3. weights_chunk : The weights of each event based on their birth redshifts and their population weight.

The ``observable_func`` should take these arguments and use them to determine the detection efficiency of the event.
We have included an example in the :func:`posydon.popsyn.transient_select_funcs.DCO_detactability`.

However, since that function requires a detection argument, it requires a wrapper to work with our function here. 

.. code-block:: python

    from posydon.popsyn.transient_select_funcs import DCO_detactability
    def DCO_wrapper(transient_chunk, z_events_chunk, weights_chunk):
      sensitivity = 'design_H1L1V1'
      return DCO_detactability(sensitivity, transient_chunk, z_events_chunk, weights_chunk, verbose=False)

    # We also give it a name, which is used as an identifier in the file
    rates.calculate_observable_population(DCO_wrapper, 'design_H1L1V1')


Accessing Observable Population data
------------------------------------

The observable population is accesed through the :code:`observable_population` attribute of the :class:`~Rates` object.
You require to know the observable_identifier to access the data, which can be accessed with :attr:`Rates.observable_population_names<Rates.observable_population_names>`

.. code-block:: python

    print(rates.observable_population_names)
    print(rates.observable_population('design_H1L1V1'))

Plotting the observable population
-----------------------------------

The observable population can be plotted in the same way as the intrinsic rate density.
However, you require to define which observable population you want to plot.

.. code-block:: python

    rates.plot_hist_properties('S1_mass', intrinsic=False, observable='design_H1L1V1', label='S1', show=True)

If you like to overplot multiple properties, you can set ``show=False`` and manually provide an axis.

.. code-block:: python

    import matplotlib.pyplot as plt
    import numpy as np
    bins = np.linspace(0,100,101)
    fig, ax = plt.subplots(1,1)

    rates.plot_hist_properties('S1_mass', intrinsice=True, observable='design_H1L1V1', bins=bins, ax = ax, label='S1', show=False)
    rates.plot_hist_properties('S2_mass', intrinsice=True, observable='design_H1L1V1', bins=bins, ax = ax, label='S1', show=False)

    

.. _population-file-structure:

The Structure of Generated Population Files 
###########################################

The main output of a population synthesis run is a HDF5 population file.

Each element in the file is stored as a pandas DataFrame.
While some elements are always present, because they're calculated as part of the population synthesis run,
other elements are optional and can be added by the user.

The tables describe the location of the data inside the population file.
This is only necessary if you want to access the data directly from the file.
If you use the :meth:`~Population` object, you can access the data directly from the object.

.. list-table:: Standard Components of a Population file
    :widths: 50 150
    :header-rows: 1

    * - Path
      - Description
    * - `history`
      - The history of each system (single star or binary) in the population, where each system has a unique index.
    * - `oneline`
      - The oneline data of each system in the population. A description of the system in a single line, which is useful for quick inspection of the population.
    * - `ini_parameters`
      - The parameters for the initial sampling conditions of the population synthesis run.
    * - `mass_per_metallicity`
      - The mass per metallicity bin for the population synthesis run. 
        The `underlying_mass` is calculated with the assumption that binary fraction == 1.

As you work with your population, you can add additional components to the population file.
Based on the components and the user given identifiers, the data is stored in the following locations in the population file.

.. list-table:: Additional components
    :widths: 50 150
    :header-rows: 1

    * - Path
      - Description
    * - `history_lengths`
      - The length of the history of each system in the population. This is created the first time the file is opened with the :class:`~Population` object.
    * - `formation_channels`
      -  The formation channels of each system in the population. This combines the `event` column in the history table into a single string. :func:`~Population.calculate_formation_channels` is used to create this component.
    * - `transiens/{transient_name}`
      - The transient data of each system in the population. The transient data is stored in a separate table for each transient. This is created by :func:`~Population.create_transient_population`.
    * - `transiens/{transient_name}/efficiencies`
      - The transient efficiencies over metallicity. This is calculated with :func:`~TransientPopulation.get_efficiency_over_metallicity`.
    * - `transiens/{transient_name}/rates/{SFH_identifier}/MODEL`
      - The MODEL parameters for the specific transient rate calculations done with :func:`~TransientPopulation.calculate_cosmic_weights`.
    * - `transiens/{transient_name}/rates/{SFH_identifier}/birth`
      - A table containing the birth redshifts and lookback times used in the rate calculation.
    * - `transiens/{transient_name}/rates/{SFH_identifier}/z_events`
      - The redshifts of the events in the population and the birth redshifts of the events.
    * - `transiens/{transient_name}/rates/{SFH_identifier}/weights`
      - The weights of each event based on their birth redshifts and their population weight.
    * - `transiens/{transient_name}/rates/{SFH_identifier}/intrinsic_rate_density`
      - The intrinsic rate density of the events in the population, calculated with :func:`~Rates.calculate_intrinsic_rate_density`.
    


