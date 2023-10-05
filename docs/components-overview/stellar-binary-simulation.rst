.. _stellar-binary-simulation:

################################
Stellar & Binary-star Simulation
################################

POSYDON: Building Populations with Object-Oriented Design
=========================================================

POSYDON provide a fully customizable framework to population synthesis modeling. By leveraging the power and flexibility of object-oriented design, it offers a modular and scalable architecture that enables researchers to customize, extend, and integrate various astrophysical processes with ease. This section delves into the intricacies of the POSYDON hierarchy, outlining the key classes, their interrelationships, and their roles in crafting sophisticated population synthesis models. 

The Single Star Object
----------------------

The `SingleStar` object contains the stellar properties and its evolutionary history. The single star object is passed to the :class:`~.binary_evol.BinaryStar` class to create a `BinaryStar` object.

    - :ref:`The Single Star Object <single-star>`

The Binary Star Object
----------------------

The `BinaryStar` object contains the binary properties and its evolutionary history as well as two `SingleStar`. The binary star object is passed to the :class:`~.populations.BinaryPopulation` class to create a population synthesis model.

    - :ref:`The Binary Star Object <binary-star>`


The Binary Population Object
----------------------------

The `BinaryPopulation` object contains a list of `BinaryStar` objects and the `SimulationProperties` objcect which contains the information about the population synthesis model.

    - :ref:`The Binary Population Object <binary-population>`

The Simulation Properties object 
--------------------------------


POSDYON Population Synthesis Configuration Guide
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The `SimulationProperties` objects contains the POSYDON population synthesis parameter configurations determining the evolution of each binaryes in the `BinaryPopulation` class. The following guide will walk you through the configuration parameters of a POSYDON binary population synthesis model.

    - :ref:`POSDYON Population Synthesis Parameters <pop-params>`


The Flow Chart Object
~~~~~~~~~~~~~~~~~~~~~

The flow chart object is the main object that is used to configure the population synthesis model. It is a dictionary that contains all the parameters that are used to configure the model. The flow chart object is passed to the :class:`~.populations.BinaryPopulation` class to create a population synthesis model.

    - :ref:`The Flow Chart Object <flow-chart>`

Evolutionary Hooks
~~~~~~~~~~~~~~~~~~

The evolutionary hooks allows the use to execute code between POSYDON steps. The following guide will walk you through the evolutionary hooks that are available in POSYDON.

    - :ref:`Evolutionary Hooks <custom-hooks>`


The Synthetic Population Object
-------------------------------

The `SyntheticPopulation` object contains a collection of `BinaryPopulation` objects run, e.g. at different metallicities or with different model assumptions. It also provide the POSYDON API interface to analyse, process, and visualize the results of a POSYDON population synthesis model.

    - :ref:`The Synthetic Population Object <synthetic-population>`


Debugging POSYDON population synthesis models
=============================================

POSYDON provides a set of tools to debug population synthesis models. The following guide will walk you through the debugging tools that are available in POSYDON.

    - :ref:`Debugging POSYDON population synthesis models <pop-syn-debugging>`

