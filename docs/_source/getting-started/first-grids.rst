.. _first-grids:

Looking at our POSYDON grids
============================

Installation and Data Download
------------------------------

First, make sure POSYDON has been installed on a conda environment. If you have not yet completed this step, please take a look at our installation guide `here <installation-guide>`_.

Additionally, POSYDON requires data to be downloaded from Zenodo. For our simple population, we will run 100 binaries at solar metallicity, so we only need the file with our grids at 1 Z☉. This file is 10 GB, so please make sure you have enough disk space available. If you have downloaded all of DR2, using the command ``get-posydon-data DR2``, then you should be good to go. If not, then for now, you only need to download the 1 Z☉ grids using the command:

.. code-block:: bash

    get-posydon-data DR2_1Zsun

Loading a grid
--------------

You can access our grids of MESA simulations using the ``PSyGrid`` module. To load a grid, first import the module, instantiate a new ``PSyGrid`` object, and pass the argument ``filepath`` to open a grid. In this example we load the HMS-HMS grid at solar metallicity:

.. code-block:: python

    from posydon.grids.psygrid import PSyGrid
    PATH_TO_POSYDON_DATA = '/YOUR/POSYDON_DATA/PATH'
    grid_path = f'{PATH_TO_POSYDON_DATA}/POSYDON_data/HMS-HMS/1e+00_Zsun.h5'
    mygrid = PSyGrid(filepath=grid_path)

Grid contents
-------------
You can check the contents of the ``PSyGrid`` object with a print command: 

.. code-block:: python

    print(mygrid)

This will provide a summary of grid metadata including how many runs are in the grid and what information is stored about them. 

You can get a list of the :ref:`grid configuration properties <tab_grid_properties>` available for your :samp:`PSyGrid` object with:

.. code-block:: python

    mygrid.config.keys()

You can access the value of any grid configuration property "PROPERTY" with 
:samp:`mygrid.config[{PROPERTY}]`.

Next, you can look at the initial and final values of the runs. All the values
are available at :samp:`mygrid.initial_values` and :samp:`mygrid.final_values`,
respectively. To get a tuple of all the available values use

.. code-block:: python

    mygrid.initial_values.dtype.names
    mygrid.final_values.dtype.names

You can access the initial value of any physical grid property "PHYS" with 
:samp:`mygrid.initial_values[{PHYS}]`. It will return a numpy array with the 
values of this property for all the runs. 
Note that these physical properties of the binaries in the grid are different 
from the grid configuration properties listed above. 
Then, you can find the initial mass of star 1 in the third MESA run with

.. code-block:: python

    mygrid.initial_values['star_1_mass'][2]

.. note::
    Remember that the first run has the index :samp:`0` and the last one
    :samp:`len(mygrid)-1`.

You can retrieve individual runs by index. :samp:`mygrid[{IDX}]` is a
:samp:`PSyRunView` object, which contains the data of the run of index 
:samp:`IDX`. The :samp:`PSyRunView` object contains seven components:

.. table:: :samp:`PSyRunView` object components

    ================  ===========
    Component         Description
    ================  ===========
    'initial_values'  all initial values of the run
    'final_values'    all final values of the run including termination flags
    'binary_history'  the binary history
    'history1'        the history of star 1
    'history2'        the history of star 2
    'final_profile1'  the final profile of star 1
    'final_profile2'  the final profile of star 2
    ================  ===========

Again, you can check for the contents of the individual runs with
:samp:`dtype.names`, e.g.

.. code-block:: python

    myrun = mygrid[0]
    myrun['binary_history'].dtype.names

The example above finds the initial mass of star 1 in the third MESA run by 
indexing the list :samp:`mygrid.initial_values`. 
You can get the same value from the list of initial values associated with a 
single MESA run: 

.. code-block :: python

    mygrid[2]['initial_values']['star_1_mass']

Plot a `PSyGrid` object
-----------------------

There are three main plotting functionalities available
to display the content of a :samp:`PSyGrid` object:

- :samp:`plot`: This creates a one-dimensional plot from the :samp:`PSyGrid`.
  An example can be found in the :ref:`tutorials <plot_1d>`. The code details
  are available in the
  :py:func:`PSyGrid.plot <posydon.grids.psygrid.PSyGrid.plot>` code and the
  :py:class:`visualization <posydon.visualization.plot1D>` library.
- :samp:`plot2D`: This creates a two-dimensional representation from the
  :samp:`PSyGrid`. Again, an example can be found in the
  :ref:`tutorials <plot_2d>`. The code details are available in the
  :py:func:`PSyGrid.plot <posydon.grids.psygrid.PSyGrid.plot2D>` code and the
  :py:class:`visualization <posydon.visualization.plot2D>` library.
- :samp:`HR`: This is similar to :samp:`plot` but specialized for producing
  Hertzsprung–Russell diagrams.

More in-depth documentation about the ``PSyGrid`` module and its full functionality is available :ref:`here <psygrid>`_. 
