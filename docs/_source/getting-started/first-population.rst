Running your First Population
=============================

Installation and Data Download
------------------------------

First, make sure POSYDON has been installed on a conda environment. If you have not yet completed this step, please take a look at our installation guide `here <installation-guide>`_.

Additionally, POSYDON requires data to be downloaded from Zenodo. For our simple population, we will run 100 binaries at solar metallicity, so we only need the file with our grids at Zsun. This file is 10 GB, so please make sure you have enough disk space available. If you have downloaded all of DR2, using the command ``get-posydon-data DR2``, then you should be good to go. If not, then for now, you only need to download the Zsun grids using the command:

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

Next, we'll create our python script to read in the .ini file and generate a population:

.. code-block:: bash

    touch my_population.py

Within this file, we'll start by loading up the necessary libraries and creating an instance of our class to handle populations:

.. code-block:: python

    from posydon.popsyn.synthetic_population import PopulationRunner

    poprun = PopulationRunner('$PATH_TO_POSYDON/posydon/popsyn/population_params_default.ini', verbose=True)
    poprun.evolve()

Now, we can evolve our population:

.. code-block:: bash

    python my_population.py

After our population has run, we can take a look at the output.
