.. _stellar-binary-simulation:

################################
Stellar & Binary-star Simulation
################################


Why the object hierarchy matters
==================================

The whole of POSYDON's population-synthesis machinery is built on a small set of
Python classes that wrap the physics of stars and binaries. If you want to
understand the code -- rather than just run it -- it pays to first soak in this
idea:

* a star is a ``SingleStar`` (its structure and evolutionary history);
* a binary is a ``BinaryStar`` (the pair plus the two ``SingleStar`` objects);
* a collection of binaries evolved together is a ``BinaryPopulation``;
* one or more populations (e.g. at different metallicities) are combined into a
  ``SyntheticPopulation``;
* the knobs that control *how* each binary evolves live in the
  ``SimulationProperties`` / flow-chart and configuration objects.

Each object delegates its behavior to the physics modules (detached, common
envelope, supernova) through the binary's evolving "state" and "event", and that
is what the rest of this section documents object by object.



POSYDON: Building Populations with Object-Oriented Design
=========================================================

POSYDON provides a fully customizable framework to population synthesis modeling. By leveraging the power and flexibility of object-oriented design, it offers a modular and scalable architecture that enables researchers to customize, extend, and integrate various astrophysical processes with ease. This section delves into the intricacies of the POSYDON hierarchy, outlining the key classes, their interrelationships, and their roles in crafting sophisticated population synthesis models.

The Single Star Object
----------------------

The `SingleStar` object contains the stellar properties and its evolutionary history. The single star object is passed to the :class:`~.binary_evol.BinaryStar` class to create a `BinaryStar` object.

.. toctree::
    :maxdepth: 1

    pop_syn/single_star

----

The Binary Star Object
----------------------

The `BinaryStar` object contains the binary properties and its evolutionary history as well as two `SingleStar`. The binary star object is passed to the :class:`~.populations.BinaryPopulation` class to create a population synthesis model.

.. toctree::
    :maxdepth: 1
    :titlesonly:

    pop_syn/binary_star


----

The Binary Population Object
----------------------------

The `BinaryPopulation` object contains a list of `BinaryStar` objects and the `SimulationProperties` object which contains the information about the population synthesis model.

.. toctree::
    :maxdepth: 1

    pop_syn/binary_population


----


The Simulation Properties object
--------------------------------


POSYDON Population Synthesis Configuration Guide
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


The `SimulationProperties` objects contains the POSYDON population synthesis parameter configurations determining the evolution of each binaries in the `BinaryPopulation` class. The following guide will walk you through the configuration parameters of a POSYDON binary population synthesis model.

.. toctree::
    :maxdepth: 1

    pop_syn/population_params


The Flow Chart Object
~~~~~~~~~~~~~~~~~~~~~

The flow chart object is the main object that is used to configure the population synthesis model. It is a dictionary that contains all the parameters that are used to configure the model. The flow chart object is passed to the :class:`~.populations.BinaryPopulation` class to create a population synthesis model.

.. toctree::
    :maxdepth: 1

    pop_syn/flow_chart


Evolutionary Hooks
~~~~~~~~~~~~~~~~~~

The evolutionary hooks allows the use to execute code between POSYDON steps. The following guide will walk you through the evolutionary hooks that are available in POSYDON.

.. toctree::
    :maxdepth: 1

    pop_syn/custom_hooks


User Functions and Modules
~~~~~~~~~~~~~~~~~~~~~~~~~~

To allow users to easier add functions for their own needs without changing the main code, we have a dedicated directory for :samp:`user_modules`.

.. toctree::
    :maxdepth: 1

    pop_syn/user_modules

----


The Synthetic Population Object
-------------------------------

The `SyntheticPopulation` object contains a collection of `BinaryPopulation` objects run, e.g. at different metallicities or with different model assumptions. It also provide the POSYDON API interface to analyse, process, and visualize the results of a POSYDON population synthesis model.

.. toctree::
    :maxdepth: 1

    pop_syn/synthetic_population


The Star Formation History
---------------------------

The star formation history is a key component in population synthesis, since it
determined the amount of stars that are formed at each moment in time.

.. toctree::
    :maxdepth: 1

    pop_syn/star_formation_history
