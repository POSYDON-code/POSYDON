.. _synthetic-population:


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
    :widths: 20 80
    :header-rows: 1

    * - Component
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
    :widths: 20 80
    :header-rows: 1

    * - Component
      - Description
    * - `history_lengths`
      - The length of the history of each system in the population. This is created the first time the file is opened with the :class:`~posydon.popsyn.synthetic_population.Population` object.
    * - `formation_channels`
      -  The formation channels of each system in the population. This combines the `event` column in the history table into a single string.
    * - `transiens/{transient_name}`


The Population Object 
===============================

The `Population` object is an interface for the population data, such that 
you can work with parts of the data and perform your own analysis. 


There are three main components of the `Population` object:
- `history` : the history of the population
- `oneline` : the oneline data of the population
- `formation_channels` : the formation channels of the population

The `formation_channels` is not immediately available, but can be created using the
:func:`~posydon.popsyn.synthetic_population.Population.calculate_formation_channels` function of the :class:`~posydon.popsyn.synthetic_population.Population` object.







ini_params

mass_per_metallicity

metallicities
solar_metallicities

indices





The Structure of Population Files
=================================





