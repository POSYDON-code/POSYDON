.. _first-population:

Quick Start: Your First Population
==================================

Installation and Data Download
------------------------------

First, make sure POSYDON has been installed on a conda environment. If you have not yet completed this step, please take a look at our installation guide `here <installation-guide>`_.

Additionally, POSYDON requires data to be downloaded from Zenodo. For our simple population, we will run 100 binaries at solar metallicity, so we only need the 1 Z☉ dataset, which can be downloaded using the ``get-posydon-data`` command in a terminal as follows: 

.. note:: 
    This dataset is 10 GB, so please make sure you have enough disk space available. If you have downloaded the full DR2 dataset using the command ``get-posydon-data DR2`` or simply ``get-posydon-data``, then you should already have the 1 Z☉ dataset and you can skip this step.

.. code-block:: bash

    get-posydon-data DR2_1Zsun

Getting everything set up
-------------------------

First, create a new working directory for your population:

.. code-block:: bash

    mkdir my_test_population

Next, set the environmental variables so they reference the directory of your POSYDON installation and POSYDON data download:

.. code-block:: bash

    %env PATH_TO_POSYDON=/YOUR/POSYDON/PATH/
    %env PATH_TO_POSYDON_DATA=/YOUR/POSYDON_DATA/PATH/

Next, create our python script to read in the .ini file and generate a population:

.. code-block:: bash

    touch my_population.py

Within this file, load up the necessary libraries and create an instance of the `PopulationRunner` class to handle the population:

.. code-block:: python

    from posydon.popsyn.synthetic_population import PopulationRunner

    poprun = PopulationRunner('$PATH_TO_POSYDON/posydon/popsyn/population_params_default.ini', verbose=True)
    poprun.evolve()

Evolve the population:

.. code-block:: bash

    python my_population.py

After the population has finished running, you can take a look at the output. More information on how to do this is available :ref:`here </components-overview/pop_syn/synthetic_population.rst>`.
