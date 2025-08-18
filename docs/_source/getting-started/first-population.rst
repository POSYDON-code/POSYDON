.. _first-population:

Running your First Population
=============================

Installation and Data Download
------------------------------

First, make sure POSYDON has been installed on a conda environment. If you have not yet completed this step, please take a look at our installation guide `here <installation-guide>`_.

Additionally, POSYDON requires data to be downloaded from Zenodo. For our simple population, we will run 100 binaries at solar metallicity, so we only need the file with our grids at 1 Z☉. This file is 10 GB, so please make sure you have enough disk space available. If you have downloaded all of DR2, using the command ``get-posydon-data DR2``, then you should be good to go. If not, then for now, you only need to download the 1 Z☉ grids using the command:

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

After the population has finished running, you take a look at the output. More information on how to do this is available :ref:`here <synthetic_population>`_
