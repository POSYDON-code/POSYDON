.. _first-grids:

Quick Start: Examining the MESA Grids
============================

Installation and Data Download
------------------------------

First, make sure POSYDON has been installed in a ``conda`` environment. If you have not yet completed this step, please take a look at our installation guide `here <installation-guide>`_.

Additionally, POSYDON requires data to be downloaded from Zenodo. For our simple population, we will run 100 binaries at solar metallicity, so we only need the 1 Z☉ dataset, which can be downloaded using the ``get-posydon-data`` command in a terminal as follows: 

.. note:: 
    This dataset is 10 GB, so please make sure you have enough disk space available. If you have downloaded the full DR2 dataset using the command ``get-posydon-data DR2`` or simply ``get-posydon-data``, then you should already have the 1 Z☉ dataset and you can skip this step.

.. code-block:: bash

    get-posydon-data DR2_1Zsun

Loading a MESA grid
--------------

You can access our grids of MESA simulations using the ``PSyGrid`` module (for an overview of this module, see :ref:`here </components-overview/post_processing/psygrid.rst>`). We will run through some simple commands to orient ourselves with the ``PSyGrid`` object below. 

To load a MESA grid, first import the module, instantiate a new ``PSyGrid`` object, and pass the argument ``filepath`` to open a grid. In this example we load the HMS-HMS grid at solar metallicity:

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
:samp:`mygrid.config[{PROPERTY}]`. Viewing this information can be useful if you are curious about how the grid was run and which options were included when postprocessing the grid (this is an advanced topic that you can read more about :ref:`here </components-overview/post_processing/pipeline.rst>`).

To maintain a relatively low data storage footprint, the PSyGrid object does not store all of the raw MESA stellar evolution data, but you can use it to look at the initial (first time step) and final (last time step) values and downsampled histories of the MESA simulations. All the initial and final simulation values
are available at :samp:`mygrid.initial_values` and :samp:`mygrid.final_values`,
respectively. To get a tuple of all the available column names use

.. code-block:: python

    mygrid.initial_values.dtype.names
    mygrid.final_values.dtype.names

You can access the initial value of any physical grid property ``var_name`` (which is a ``str`` corresponding to one of the available data columns above) with 
:samp:`mygrid.initial_values[{var_name}]`. Doing this will return a ``numpy`` array with the 
initial values of this property for all the MESA runs in this grid. 
Note that these physical properties of the binaries in the grid are different 
from the grid configuration properties listed above and are the actual data outputs of the MESA simulations. 

As an example, you can find the initial mass of the primary star (``star_1_mass``) in one of the grid's MESA binary star simulations with

.. code-block:: python

    mygrid.initial_values['star_1_mass'][2]

.. note::
    In the example, we used an arbitrary index `2` to select a binary, but the indices run from :samp:`0` to 
    :samp:`len(mygrid)-1`.

You can also retrieve individual runs (binaries) by index using :samp:`mygrid[{i}]`. This will return a
:samp:`PSyRunView` object which contains data pertaining to that particular run of index 
:samp:`i`. The :samp:`PSyRunView` object contains seven components:

.. table:: :samp:`PSyRunView` object components

    ================  ===========
    Component         Description
    ================  ===========
    'initial_values'  all initial (i.e., at the first time step of the MESA simulation) values of the run
    'final_values'    all final (i.e., at the final time step of the MESA simulation) values of the run. (This also includes termination flags, which you can read more about :ref:`here </tutorials-examples/generating-datasets/plot_2d.ipynb>`.)
    'binary_history'  the downsampled binary star evolution data
    'history1'        the downsampled evolution data of star 1
    'history2'        the downsampled evolution data of star 2
    'final_profile1'  the profile (properties vs. mass or radius) of star 1 at the final timestep
    'final_profile2'  the profile of star 2 at the final timestep
    ================  ===========

Again, you can check for the valid column names of the individual runs with
:samp:`dtype.names`, e.g.

.. code-block:: python

    myrun = mygrid[0]
    myrun['binary_history'].dtype.names

Before, we retrieved the initial mass of star 1 in the binary indexed at ` 2` by
slicing the list :samp:`mygrid.initial_values` and the selecting the value at index `2`. 
You can get the same value from the list of initial values associated with a 
single MESA run using the `PSyRunView` object via: 

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

More in-depth documentation about the ``PSyGrid`` module and its full functionality is available :ref:`here <psygrid>`. 
