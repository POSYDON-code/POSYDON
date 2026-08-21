.. _population_synthesis:

Population Synthesis
====================
The population synthesis framework is built on top of the `_stellar-binary-simulation`_ framework, and
evolves a population of binaries using the `BinaryPopulation`_ class.

It uses the same configuration file as the `_stellar-binary-simulation`_ framework, but it also allows to specify the initial conditions of the population, such as the initial mass function, the initial mass ratio distribution, and the initial orbital period distribution. See `SimulationProperties`_ for more information on the configuration parameters.



The Binary Population Object
----------------------------

The `BinaryPopulation` object contains a list of `BinaryStar` objects and the `SimulationProperties` object which contains the information about the population synthesis model.

.. toctree::
    :maxdepth: 1

    pop_syn/binary_population


----

The Synthetic Population Object
-------------------------------

The `SyntheticPopulation` object contains a collection of `BinaryPopulation` objects run, e.g. at different metallicities or with different model assumptions. It also provide the POSYDON API interface to analyse, process, and visualize the results of a POSYDON population synthesis model.

.. toctree::
    :maxdepth: 1

    pop_syn/synthetic_population


Reweighting Populations
-----------------------
POSYDON provides a reweighting framework that allows users to reweight the results of a population synthesis model to match a different set of initial conditions. This is useful for exploring the effects of different initial conditions on the results of a population synthesis model without having to run a new simulation.

.. toctree::
    :maxdepth: 1

    pop_syn/reweighting



The Star Formation History
---------------------------

The star formation history is a key component in population synthesis, since it
determined the amount of stars that are formed at each moment in time.

.. toctree::
    :maxdepth: 1

    pop_syn/star_formation_history
