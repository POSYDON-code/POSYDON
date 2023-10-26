.. _psygrid:

##############
PSyGrid object
##############

The `PSyGrid` object contains the data coming from detailed MESA simulations.

Import a `PSyGrid` object by using:

.. code-block:: python

  from posydon.grids.psygrid import PSyGrid

Initializing a `PSyGrid` object
-------------------------------
One can simply get a new `PSyGrid` object via

.. code-block:: python

    mygrid = PSyGrid()

This calls the initializer, which can take two optional arguments

- :samp:`filepath`: it is the path to the associated h5 file; if the file exists the `PSyGrid` object in there will be loaded
- :samp:`verbose`: it is a boolean indicating whether detailed output should be given when calling functions of the `PSyGrid` object

Creating a `PSyGrid` object
---------------------------
The `create` function of a `PSyGrid` object will read MESA output data into the
`PSyGrid` object. Thus it has a required argument and several optional
arguments. The required one is :samp:`MESA_grid_path`, which is the path to the
directory, where the MESA runs are in.

.. table:: Optional arguments of the PSyGrid creation
    :widths: 20,10,90

    =============  =========  ===========
    Argument       Default    Description
    =============  =========  ===========
    psygrid_path   None       the path to the associated h5 file (it needs to be given if none was specified during initialization)
    overwrite      False      if `True` overwrite a potentially already existing file
    slim           False      if `True` only inital and final values as well as the meta data will be stored
    warn           "end"      if "normal" warnings are printed, when they arise
                              
                              if "end" all warnings are printed at the end
                              
                              if "suppress" no warnings are printed
    fmt            "posydon"  grid format; only "posydon" is currently supported
    **grid_kwargs             further grid properties can be specified in a dictionary
    =============  =========  ===========

.. table:: Grid properties

    ==============================  ============  ===========
    Property                        Default       Description
    ==============================  ============  ===========
    'description'                   ""            a description text
    'max_number_of_runs'            None          the maximum number of runs
    'format'                        "hdf5"        file format; only "hdf5" is currently supported
    'compression'                   "gzip9"       the compression (of the hdf5 file)
    'history_DS_error'              None          the maximum error allowed when downsampling the history
    'history_DS_exclude'            default list  the history columns to exclude from downsampling (default list: ["model_number", "age", "star_age"])
    'profile_DS_error'              None          the maximum error allowed when downsampling the final profile
    'profile_DS_interval'           None          the maximum change in an downsampled interval relative to the change from initial to final
    'profile_DS_exclude'            default list  the profile columns to exclude from downsampling (default list: ["mass" "star_mass",])
    'star1_history_saved_columns'   "minimum"     specifies which history columns of star 1 should be read
                                                  
                                                  if "all" read all the columns in the MESA output
                                                  
                                                  if "minimum" use the default
                                                  
                                                  if a tuple of column names read only those columns
                                                  
                                                  if a list of column names read the default and those columns
    'star2_history_saved_columns'   "minimum"     specifies which history columns of star 2 should be read having the same options as 'star1_history_saved_columns'
    'binary_history_saved_columns'  "minimum"     specifies which binary history columns should be read having the same options as 'star1_history_saved_columns'
    'star1_profile_saved_columns'   "minimum"     specifies which profile columns of star 1 should be read having the same options as 'star1_history_saved_columns'
    'star2_profile_saved_columns'   "minimum"     specifies which profile columns of star 2 should be read having the same options as 'star1_history_saved_columns'
    'initial_value_columns'         None          history columns to store initial values from (currently not in use, instead all specified history columns are used and additionally the abundances X, Y, and Z)
    'final_value_columns'           None          history columns to store final values from (currently not in use, instead all specified history columns are used and additionally termination flags and for binaries the interpolation class)
    'start_at_RLO'                  False         specifies whether to crop the history to start at RLO
    'stop_before_carbon_depletion'  False         specifies whether to crop the history of massive stars (>100 Msun) to stop at 10% central carbon and after helium is depleted
    'binary'                        True          specifies whether a gird evolved binaries; put `False` for single stars
    'eep'                           None          path to directory with EEP files (for single stars only)
    'initial_RLO_fix'               False         specifies whether the boundary of initial RLO should be determined to flag all systems below as initial RLO independent of the MESA output
    'He_core_fix'                   True          specifies to ensure that the He core is always larger or equal to the carbon-oxygen core
    'accept_missing_profile'        False         specifies whether try to include all data from MESA runs without final profiles
    ==============================  ============  ===========

You can read the MESA data into an existing `PSyGrid` object, which may
overwrites data:

.. code-block:: python

    mygrid.create(MESA_grid_path=".")

or combine the initialization with the creation:

.. code-block:: python

    mygrid = PSyGrid().create(MESA_grid_path=".")

Loading a `PSyGrid` object
--------------------------
You can load an existing h5 file (e.g. "myPSyGrid.h5") into a `PSyGrid` object
by

.. code-block:: python

    mygrid.load(filepath="myPSyGrid.h5")

It is more convinient to load the file directly when initializing the `PSyGrid`
object

.. code-block:: python

    mygrid = PSyGrid(filepath="myPSyGrid.h5")

Contents of a `PSyGrid` object
------------------------------

Print a `PSyGrid` object
~~~~~~~~~~~~~~~~~~~~~~~~
TODO: __str__()
TODO: __len__()

Accessing data in a `PSyGrid` object
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TODO: __getitem__(), PSyRunView()

Plot a `PSyGrid` object
~~~~~~~~~~~~~~~~~~~~~~~
TODO: plot
TODO: plot2D
TODO: HR

Work on/with a `PSyGrid` object
-------------------------------

Loop on a `PSyGrid` object
~~~~~~~~~~~~~~~~~~~~~~~~~~
TODO: __iter__()

Expand a `PSyGrid` object
~~~~~~~~~~~~~~~~~~~~~~~~~
TODO: add_column()

Join two or more `PSyGrid` objects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TODO: join_grids()

Extract the initial and final values as a pandas data frame
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TODO: get_pandas_initial_final()

Get reruns from a `PSyGrid` object
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TODO: rerun()

Close associated hdf5 file
~~~~~~~~~~~~~~~~~~~~~~~~~~
TODO: close()
TODO: __del__()

The code summary of the `PSyGrid` object can be found at the
:ref:`dedicated reference page <code_psygrid>`.
