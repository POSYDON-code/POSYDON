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

- :samp:`filepath`: it is the path to the associated h5 file; if the file
  exists the `PSyGrid` object in there will be loaded
- :samp:`verbose`: it is a boolean indicating whether detailed output should be
  given when calling functions of the `PSyGrid` object


Creating a `PSyGrid` object
---------------------------

The `create` function of a `PSyGrid` object will read MESA output data into the
`PSyGrid` object. Thus it has a required argument and several optional
arguments. The required one is :samp:`MESA_grid_path`, which is the path to the
directory, where the MESA runs are in.

.. table:: Optional arguments of the PSyGrid creation
    :widths: 18,10,72

    =============  =========  ===========
    Argument       Default    Description
    =============  =========  ===========
    psygrid_path   None       the path to the associated h5 file (it needs to be given if none was specified during initialization)
    overwrite      False      if `True` overwrite a potentially already existing file
    slim           False      if `True` only initial and final values as well as the meta data will be stored
    warn           "end"      if "normal" warnings are printed, when they arise
                              
                              if "end" all warnings are printed at the end
                              
                              if "suppress" no warnings are printed
    fmt            "posydon"  grid format; only "posydon" is currently supported
    **grid_kwargs             further grid properties can be specified in a dictionary
    =============  =========  ===========

.. _tab_grid_properties:

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
    'profile_DS_exclude'            default list  the profile columns to exclude from downsampling (default list: ["mass", "star_mass"])
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

It is more convenient to load the file directly when initializing the `PSyGrid`
object

.. code-block:: python

    mygrid = PSyGrid(filepath="myPSyGrid.h5")


Contents of a `PSyGrid` object
------------------------------

Print a `PSyGrid` object
~~~~~~~~~~~~~~~~~~~~~~~~

To check the content of the `PSyGrid` object you can simply print it:

.. code-block:: python

    print(mygrid)

This will provide a summary, which tell you

- to which hdf5 file it is connected
- how many runs are in there and have
 
  - a binary history
  - a history of star 1
  - a history of star 2
  - a final profile of star 1
  - a final profile of star 2
   
- the fields in the each of the histories/profiles of the last run
- the fields of the initial and final values
- information on the configuration
- a shorthand list of the MESA directories (the locations of the data the runs
  where extracted from)

To access single runs, it is important to know how many are there to avoid to
call for a non existing one. Hence, you can simply get the number of runs via:

.. code-block:: python

    len(mygrid)

.. note::
    Alternatively you can request the length from the internal number stored in
    :samp:`mygrid.n_runs`.


Accessing data in a `PSyGrid` object
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The first data you may want to check are the
:ref:`grid properties <tab_grid_properties>`. You can get a list of the
properties available for your `PSyGrid` object simply with

.. code-block:: python

    mygrid.config.keys()

By providing one of the properties to :samp:`mygrid.config[{PROPERTY}]` you can
access its value.

Next, you can look at the initial and final values of the runs. All the values
are available at :samp:`mygrid.initial_values` and :samp:`mygrid.final_values`,
respectively. To get a tuple of all the available values use

.. code-block:: python

    mygrid.initial_values.dtype.names
    mygrid.final_values.dtype.names

Each value you then get for example via :samp:`mygrid.initial_values[{VALUE}]`.
It will return a numpy array with the this value for all the runs. So you get
the initial mass of star 1 in the third run with

.. code-block:: python

    mygrid.initial_values['star_1_mass'][2]

.. note::
    Remember, that the first run has the index :samp:`0` and the last one
    :samp:`len(mygrid)-1`.

The each initial and final value will have the same number and order of run
entries. This holds for the number of list entries of MESA directories, too.

.. code-block:: python

    mygrid.MESA_dirs

You can get the individual runs via its index. :samp:`mygrid[{IDX}]` is a
`PSyRunView` object, which contains the data of the run specified with `IDX`.
The `PSyRunView` object contains seven components:

.. table:: `PSyRunView` object components

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

Again you can check for the contents of the components with
:samp:`dtype.names`, e.g.

.. code-block:: python

    myrun = mygrid[0]
    myrun['binary_history'].dtype.names

Now you know a second way to get the initial mass of star 1 in the third run
with a one-liner

.. code-block :: python

    mygrid[2]['initial_values']['star_1_mass']

You would may think of a third way being

.. code-block :: python

    mygrid[2]['binary_history']['star_1_mass'][0]

But this will not give the initial value, while it is close to it. The reason
for this not being the same is in MESA having a slightly different value in the
first line of the history files compared to the given initial value. The final
values and the derived initial values instead are the same as the last or first
values in the corresponding history.

.. note::
    For efficiency reasons not all the `PSyGrid` object is loaded into RAM.
    Instead parts are reads from the associated hdf5 file if needed. This has
    the consequence, that it is discouraged to refer to the same values more
    than once in a code. If you need the same value more often, you should
    store it in a local variable.


Plot a `PSyGrid` object
~~~~~~~~~~~~~~~~~~~~~~~

Beside getting the values itself there are plotting functionalities available
to display the content of a `PSyGrid` object. There are three main plotting
functionalities:

- :samp:`plot`: This creates a one dimensional plot from the `PSyGrid`. An
  example can be found in the :ref:`tutorials <plot_1d>`. The code details are
  available in the :ref:`PSyGrid code <code_psygrid>` and the
  :ref:`visualisation libary <vis_plot1D>`.
- :samp:`plot2D`: This creates a two dimensional representation from the
  `PSyGrid`. Again, an example can be found in the :ref:`tutorials <plot_2d>`.
  The code details are available in the :ref:`PSyGrid code <code_psygrid>`
  and the :ref:`visualisation libary <vis_plot2D>`.
- :samp:`HR`: This is similar to :samp:`plot` but specialized on producing
  Hertzsprung–Russell diagrams.


Work on/with a `PSyGrid` object
-------------------------------

Loop on a `PSyGrid` object
~~~~~~~~~~~~~~~~~~~~~~~~~~

Similarly to accessing a single value in the `PSyGrid` object we can loop over
a `PSyGrid` object, which will loop over the individual runs in the `PSyGrid`
object. Hence the following two codes will produce the same output. The first
one loops through the numpy array if the initial companion mass

.. code-block:: python

    for mass in mygrid.initial_values['star_2_mass']:
        print(mass)

while the second one loops through the runs and prints the initial companion
mass

.. code-block:: python

    for run in mygrid:
        print(run['initial_values']['star_2_mass'])



Expand a `PSyGrid` object
~~~~~~~~~~~~~~~~~~~~~~~~~

Because of the complexity of the `PSyGrid` object we encourage the user to only
use our dedicated functions to add content to the object. There is a function
to add an extra column to the :samp:`final_values`. Here is an example how to
add a new column which contains the final orbital period in units of years
instead of days:

.. code-block:: python

    new_column_data = mygrid.final_values['period_days']/365.25
    mygrid.add_column('period_years', new_column_data, where='final_values', overwrite=False)

The four arguments are a string with the name of the new field, the data to be
stored in there, to which component of the `PSyGrid` object it should get
added, and whether a field with the same name should be overwritten, if it
already exists.

.. warning::
    The new data has to have as many entries as the `PSyGrid` object has runs.

.. note::
    Currently, the parameter :samp:`where` only supports the value 'final_values'.


Join two or more `PSyGrid` objects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are different reasons, why you create several `PSyGrid` objects which you
would like to combine to a single one later, e.g. adding reruns. There is a
functionality to do this for you. To avoid too many conflicts of possible
modifications of already loaded `PSyGrid` objects, this function is not part of
the `PSyGrid`-object class. Instead it take a list of paths to the hdf5 files
containing `PSyGrid` objects to be combined to a new one. Those are the two
required arguments of the :samp:`join_grids` function. Additionally, you can
specify the arguments :samp:`compression`, :samp:`description`, and
:samp:`verbose`. The :samp:`join_grids` function will check, whether the grids
are compatible and join them if possible.

.. note::
    If there are common systems in two or more grids, this routine will only
    put the last run with same initial conditions in the newly combined
    `PSyGrid` object.

We recommend to use the :ref:`post-processing pipeline <pipeline>` to create
and join grids.

..
    Extract the initial and final values as a pandas data frame
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TODO: get_pandas_initial_final()


Get reruns from a `PSyGrid` object
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Usually, not all runs of a grid will be successfully run in MESA. Hence, one
may wants to rerun some of them with changed parameters. There is a function,
to export runs from a `PSyGrid` object. There are two general way to specify,
which systems should be exported to rerun:

1. Write your own logic and create a numpy array with the indexes of the systems, you would like to run again.
2. Specify, which termination flag(s) the systems should have to be rerun.

.. table:: Arguments of the :samp:`rerun` function

    =================  =======  ===========
    Argument           Default  Description
    =================  =======  ===========
    path_to_file       './'     where to create the file(s) for the rerun
    runs_to_rerun      None     a numpy array containing the indexes of the runs in the `PSyGrid` object (if given, leave :samp:`termination_flags=None`)
    termination_flags  None     a single termination flag code or a list of them (if given, leave :samp:`runs_to_rerun=None`)
    new_mesa_flag      None     dictionary with the names and the values of MESA parameters to be changed for the inlists of the new runs
    =================  =======  ===========

The :ref:`post-processing pipeline <pipeline_stepR>` already provides some pre
defined rerun options.


Close associated hdf5 file
~~~~~~~~~~~~~~~~~~~~~~~~~~

Finally, you can close the hdf5 file, which is recommended to ensure that all
your changes on the `PSyGrid` object is safely written into the file.

.. code-block:: python

    mygrid.close()

This is done as well, in the case you call the destructor of the `PSyGrid`
object.

.. code-block:: python

    del mygrid


The code summary of the `PSyGrid` object can be found at the
:ref:`dedicated reference page <code_psygrid>`.
