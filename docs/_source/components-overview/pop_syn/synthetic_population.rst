.. _synthetic-population:



The Population Object 
===============================

The :class:`~posydon.popsyn.synthetic_population.Population` object is an interface for the population data, which is stored in a HDF5 file.
You can find more information about the file structure here: :ref:`population-file-structure`.

The population file is created by a population synthesis run, and contains the history of each system in the population, and an oneline description of each system and its evolution.

The :class:`~posydon.popsyn.synthetic_population.Population` object is created by passing the path to the population file to the constructor.
It will not load in the data immediately, but only when requested.

.. code-block:: python

    from posydon.popsyn import Population

    pop = Population('path/to/population.h5')



There are three main components of the :class:`~posydon.popsyn.synthetic_population.Population` object:

- :class:`~posydon.popsyn.synthetic_population.Population.history`: the history of the population
- :class:`~posydon.popsyn.synthetic_population.Population.oneline`: the oneline data of the population
- :class:`~posydon.popsyn.synthetic_population.Population.formation_channels`: the formation channels of the population

The :code:`formation_channels` is not immediately available, but can be created using the
:func:`~posydon.popsyn.synthetic_population.Population.calculate_formation_channels` function of the :class:`~posydon.popsyn.synthetic_population.Population` object.


history
--------

:class:`~posydon.popsyn.synthetic_population.Population.history`contains the full history of each system in the population.
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


A more powerful feature is the :func:`posydon.popsyn.synthetic_population.History.select` function, which allows you to select specific rows or columns from the history table. Here's an example:

.. note::
    The :func:`posydon.popsyn.synthetic_population.History.select` function only allows the use of :code:`where=`
    for specific columns and the index. The columns are limited to those containing strings.

.. code-block:: python

    # using where with the index
    mass_10 = pop.history.select(columns=['S1_mass'], where='index==10')

    # using where with a string column
    mass_ZAMS = pop.history.select(columns=['S1_mass'], where='event == "ZAMS"')
    
If you want to have a peak, you can use the :meth:`~posydon.popsyn.synthetic_population.Population.head` or :meth:`~posydon.popsyn.synthetic_population.Population.tail` functions.

.. code-block:: python
    pop.history.head(10)
    pop.history.tail(10)


Additional functions are made available for easy of use.

If you want to check the length of the history of a system, you can use :attr:`Population.history.lengths<posydon.popsyn.synthetic_population.History.lengths>` or :attr:`Population.history_lengths<posydon.popsyn.synthetic_population.Population.history_lengths>`.

.. code-block:: python

    print(pop.history.lengths)
    print(pop.history_lengths)

The total number of systems in the population can be found with :attr:`~posydon.popsyn.synthetic_population.Population.History.number_of_systems`.

.. code-block:: python

    print(pop.history.number_of_systems)

Similarly, if you would like to check the indices in the file, you can use :attr:`Population.indices<posydon.popsyn.synthetic_population.Population.indices>` or :attr:`Population.history.indices<posydon.popsyn.synthetic_population.History.indices>`.
The indices are useful in selecting systems from the population.

It's also possible to check the columns in the history table with :attr:`Population.columns<posydon.popsyn.synthetic_population.Population.columns>` or :attr:`Population.history.columns<posydon.popsyn.synthetic_population.History.columns>`.

.. code-block:: python

    print(pop.indices)
    print(pop.history.indices)

    print(pop.history.columns)
    print(pop.columns['history'])


oneline
--------

:meth:`~posydon.popsyn.synthetic_population.Population.oneline` contains a single line description of each system in the population.
This is useful for a quick inspection of the population.
The oneline data is stored in a pandas DataFrame, where each row is a system in the population.

Some properties over the evolution of the binary do not change, such as the natal kick properties or interpolation class.
Besides the initial and final properties of the system, this table also contains these data.

The initial-final properties are those in the history table, but with the postfix :code:`_i` and :code:`_f` depending on the initial or final value.
The additional values are the scalar values from the individual stars and the binary properties (See :ref:`Single Star<single-star>` and :ref:`Binary Star<binary-star>`).

Additionally, WARNING, FAILED, and metallicity columns are available in the oneline table.

.. csv-table:: Additional columns
  :header: "Properties", "Descriptions"
  :widths: 50, 150
  `FAILED`, Indicates if the system failed during the population synthesis run.
  `WARNING`, Indicates if there were any warnings for the system during the population synthesis run.
  `metallicity`, The metallicity of the system.


Like the :code:`history` access, you can access the oneline data by using the index of the system or the columns.

.. code-block:: python

    binary_0 = pop.oneline[0]
    mass = pop.oneline['S1_mass_i']
    selection = pop.oneline.select(columns=['S1_mass_i'], where='index==10')

You can check the columns and indices of the oneline table with :attr:`Population.oneline.columns<posydon.popsyn.synthetic_population.Oneline.columns>` and :attr:`Population.columns['oneline']<posydon.popsyn.synthetic_population.Population.columns`.

.. code-block:: python

    print(pop.oneline.columns)
    print(pop.columns['oneline'])

The number of systems in the population can be found with :attr:`Population.oneline.number_of_systems<posydon.popsyn.synthetic_population.Oneline.number_of_systems>`.
The length and indices of the oneline table can be found with :attr:`Population.oneline.lengths<posydon.popsyn.synthetic_population.Oneline.lengths>`, and :attr:`Population.oneline.indices<posydon.popsyn.synthetic_population.Oneline.indices>`, respectively.


.. code-block:: python

    print(pop.oneline.number_of_systems)
    print(pop.oneline.lengths)
    print(pop.oneline.indices)
  



formation_channels
------------------

:class:`~posydon.popsyn.synthetic_population.Population.formation_channels` contains the formation channels of each system in the population.
The formation channels are stored in a pandas DataFrame, where each row is a system in the population.

The formation channels are calculated by combining the `event` column in the history table into a single string using the :func:`~posydon.popsyn.synthetic_population.Population.calculate_formation_channels` function of the :class:`~posydon.popsyn.synthetic_population.Population` object.

Two columns are available in the formation channels table:

- `debug_channel` : A longer description of the formation channel, where additional events are included.

- `channel` : A cleaned-up version of the history events, where events are separated by a `-`. 


# Exporting part of the population


# Transient population creation functions


# Star Foramtion Rate functions


# Observability functions








.. _population-file-structure:

The Structure of Population Files 
=================================

The main output of a population synthesis run is a HDF5 population file.

Each element in the file is stored as a pandas DataFrame.
While some elements are always present, because they're calculated as part of the population synthesis run,
other elements are optional and can be added by the user.

The tables describe the location of the data inside the population file.
This is only necessary if you want to access the data directly from the file.
If you use the :meth:`~posydon.popsyn.synthetic_population.Population` object, you can access the data directly from the object.

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
      - The length of the history of each system in the population. This is created the first time the file is opened with the :class:`~posydon.popsyn.synthetic_population.Population` object.
    * - `formation_channels`
      -  The formation channels of each system in the population. This combines the `event` column in the history table into a single string. :func:`~posydon.popsyn.synthetic_population.Population.calculate_formation_channels` is used to create this component.
    * - `transiens/{transient_name}`
      - The transient data of each system in the population. The transient data is stored in a separate table for each transient. This is created by :func:`~posydon.popsyn.synthetic_population.Population.create_transient_population`.
    * - `transiens/{transient_name}/efficiencies`
      - The transient efficiencies over metallicity. This is calculated with :func:`~posydon.popsyn.synthetic_population.TransientPopulation.get_efficiency_over_metallicity`.
    * - `transiens/{transient_name}/rates/{SFH_identifier}/MODEL`
      - The MODEL parameters for the specific transient rate calculations done with :func:`~posydon.popsyn.synthetic_population.TransientPopulation.calculate_cosmic_weights`.
    * - `transiens/{transient_name}/rates/{SFH_identifier}/birth`
      - A table containing the birth redshifts and lookback times used in the rate calculation.
    * - `transiens/{transient_name}/rates/{SFH_identifier}/z_events`
      - The redshifts of the events in the population and the birth redshifts of the events.
    * - `transiens/{transient_name}/rates/{SFH_identifier}/weights`
      - The weights of each event based on their birth redshifts and their population weight.
    * - `transiens/{transient_name}/rates/{SFH_identifier}/intrinsic_rate_density`
      - The intrinsic rate density of the events in the population, calculated with :func:`~posydon.popsyn.synthetic_population.Rates.calculate_intrinsic_rate_density`.
    


