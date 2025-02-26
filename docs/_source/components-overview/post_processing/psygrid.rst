.. _psygrid:

################
`PSyGrid` object
################

The :samp:`PSyGrid` object contains the data coming from detailed MESA
simulations.

Import the :samp:`PSyGrid` module with:

.. code-block:: python

  from posydon.grids.psygrid import PSyGrid


Initializing a `PSyGrid` object
-------------------------------

You can instantiate a new :samp:`PSyGrid` object with

.. code-block:: python

    mygrid = PSyGrid()

This calls the initializer, which can take two optional arguments

- :samp:`filepath`: this is the path to the associated h5 file; if the file
  exists, the initializer will load the contained `PSyGrid` object
- :samp:`verbose`: this is a boolean indicating whether detailed output should be
  given when calling functions of the :samp:`PSyGrid` object


Creating a `PSyGrid` object
---------------------------

The :samp:`create` function of a :samp:`PSyGrid` object will read MESA output
data into the :samp:`PSyGrid` object. It has a required argument, :samp:`MESA_grid_path`, 
which is the path to the directory containing the MESA runs. There are several optional 
arguments as well:

.. table:: Optional arguments of the :samp:`PSyGrid` creation
    :widths: 18,10,72

    ===============  =========  ===========
    Argument         Default    Description
    ===============  =========  ===========
    psygrid_path     None       the path to the associated h5 file (required if none was specified during initialization)
    overwrite        False      if :samp:`True`, overwrite any file that exists
    slim             False      if :samp:`True`, only store initial/final values and metadata
    warn             "end"      if :samp:`normal`, warnings are printed, when they arise
                              
                                if :samp:`end`, all warnings are printed at the end
                              
                                if :samp:`suppress`, no warnings are printed
    fmt              "posydon"  grid format; only :samp:`posydon` is currently supported
    \*\*grid_kwargs             further grid properties can be specified in a dictionary
    ===============  =========  ===========

The :samp:`PSyGrid` object has the following grid properties: 

.. _tab_grid_properties:

.. table:: Grid properties

    ==============================  ============  ===========
    Property                        Default       Description
    ==============================  ============  ===========
    'description'                   ""            description text
    'max_number_of_runs'            None          the maximum number of runs
    'format'                        "hdf5"        file format; only :samp:`hdf5` is currently supported
    'compression'                   "gzip9"       the compression (of the hdf5 file)
    'history_DS_error'              None          the maximum error allowed when downsampling the history
    'history_DS_exclude'            default list  the history columns to exclude from downsampling (default list = ["model_number", "age", "star_age"])
    'profile_DS_error'              None          the maximum error allowed when downsampling the final profile
    'profile_DS_interval'           None          the maximum change in an downsampled interval relative to the change from initial to final
    'profile_DS_exclude'            default list  the profile columns to exclude from downsampling (default list = ["mass", "star_mass"])
    'star1_history_saved_columns'   "minimum"     specifies which history columns of star 1 should be read
                                                  
                                                  if :samp:`all`, read all the columns in the MESA output
                                                  
                                                  if :samp:`minimum`, use the default
                                                  
                                                  if a tuple of column names, read only those columns
                                                  
                                                  if a list of column names, read the default and those columns
    'star2_history_saved_columns'   "minimum"     specifies which history columns of star 2 should be read (same options as star1_history_saved_columns)
    'binary_history_saved_columns'  "minimum"     specifies which binary history columns should be read (same options as star1_history_saved_columns)
    'star1_profile_saved_columns'   "minimum"     specifies which profile columns of star 1 should be read (same options as star1_history_saved_columns)
    'star2_profile_saved_columns'   "minimum"     specifies which profile columns of star 2 should be read (same options as star1_history_saved_columns)
    'initial_value_columns'         None          history columns from which to store initial values (currently not in use, instead all specified history columns are used as well as the abundances X, Y, and Z)
    'final_value_columns'           None          history columns from which to store final values (currently not in use, instead all specified history columns are used as well as termination flags and for binaries the interpolation class)
    'start_at_RLO'                  False         specifies whether to crop the history to start at RLO
    'stop_before_carbon_depletion'  False         specifies whether to crop the history of massive stars (>100 Msun) to stop at 10% central carbon and after helium is depleted
    'binary'                        True          specifies whether a grid evolved binaries; put :samp:`False` for single stars
    'eep'                           None          path to directory with EEP files (for single stars only)
    'initial_RLO_fix'               False         specifies whether the boundary of initial RLO should be determined to flag all systems below as initial RLO independent of the MESA output
    'He_core_fix'                   True          specifies to ensure that the helium core is always larger or equal to the carbon-oxygen core
    'accept_missing_profile'        False         specifies whether try to include all data from MESA runs without final profiles
    ==============================  ============  ===========

You can read the MESA data into an existing :samp:`PSyGrid` object, which may
overwrite data:

.. code-block:: python

    mygrid.create(MESA_grid_path=".")

Alternatively, you can combine the initialization with creation of the grid based on MESA data:

.. code-block:: python

    mygrid = PSyGrid().create(MESA_grid_path=".")


Loading a `PSyGrid` object
--------------------------

You can load an existing h5 file (e.g. "myPSyGrid.h5") into a :samp:`PSyGrid`
object:

.. code-block:: python

    mygrid.load(filepath="myPSyGrid.h5")

It may be more convenient to load the file directly when initializing the
:samp:`PSyGrid` object

.. code-block:: python

    mygrid = PSyGrid(filepath="myPSyGrid.h5")


Contents of a `PSyGrid` object
------------------------------

Print a `PSyGrid` object
~~~~~~~~~~~~~~~~~~~~~~~~

You can check the contents of the :samp:`PSyGrid` object with a print command:

.. code-block:: python

    print(mygrid)

This will provide a summary, which tell you:

- to which hdf5 file it is connected
- how many runs are in the grid and how many have
 
  - a binary history
  - a history of star 1
  - a history of star 2
  - a final profile of star 1
  - a final profile of star 2
   
- which histories/profile fields are included in the last run
- which initial and final value are stored in the grid
- information about the grid configuration
- a shorthand list of the MESA directories (the locations of the data the runs
  where extracted from)

To access single runs, it is important to know how many are there to avoid calling a nonexisting run. You can find the number of runs with:

.. code-block:: python

    len(mygrid)

.. note::
    Alternatively, you can request the length internally with
    :samp:`mygrid.n_runs`.


Accessing data in a `PSyGrid` object
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The first data you may want to check are the
:ref:`grid properties <tab_grid_properties>`. You can get a list of the
properties available for your :samp:`PSyGrid` object with

.. code-block:: python

    mygrid.config.keys()

You can access the value of any grid configuration property "PROP" with :samp:`mygrid.config[{PROP}]`.

Next, you can look at the initial and final values of the runs. All the values
are available at :samp:`mygrid.initial_values` and :samp:`mygrid.final_values`,
respectively. To get a tuple of all the available values use

.. code-block:: python

    mygrid.initial_values.dtype.names
    mygrid.final_values.dtype.names

You can access the initial value of any individual grid property "PROP" with :samp:`mygrid.initial_values[{PROP}]`.
It will return a numpy array with the values of this property for all the runs. 
Then, you can find the initial mass of star 1 in the third MESA run with

.. code-block:: python

    mygrid.initial_values['star_1_mass'][2]

.. note::
    Remember that the first run has the index :samp:`0` and the last one
    :samp:`len(mygrid)-1`.

Each grid property will have the same number and order of MESA run entries in the initial and final values.
This holds for the list of MESA directories from which the runs are extracted, too.

.. code-block:: python

    mygrid.MESA_dirs

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
You can get the same value from the list of initial values associated with a single MESA run: 

.. code-block :: python

    mygrid[2]['initial_values']['star_1_mass']

You can get something close to the initial value with:

.. code-block :: python

    mygrid[2]['binary_history']['star_1_mass'][0]

But this will not give the initial value, while it is close to it. 
This is because MESA has a slightly different value in the
first line of the history files compared to the given initial value. The final
values and the derived initial values instead are the same as the last or first
values in the corresponding history.

.. note::
    For efficiency reasons not all the :samp:`PSyGrid` object is loaded into
    RAM. Instead parts are reads from the associated hdf5 file if needed. 
    For this reason, it is discouraged to refer to the same values
    more than once in a code. If you need the same value more often, you should
    store it in a local variable.


Plot a `PSyGrid` object
~~~~~~~~~~~~~~~~~~~~~~~

Beside getting the values itself there are plotting functionalities available
to display the content of a :samp:`PSyGrid` object. There are three main plotting
functionalities:

- :samp:`plot`: This creates a one dimensional plot from the :samp:`PSyGrid`.
  An example can be found in the :ref:`tutorials <plot_1d>`. The code details
  are available in the
  :py:func:`PSyGrid.plot <posydon.grids.psygrid.PSyGrid.plot>` code and the
  :py:class:`visualization <posydon.visualization.plot1D>` library.
- :samp:`plot2D`: This creates a two dimensional representation from the
  :samp:`PSyGrid`. Again, an example can be found in the
  :ref:`tutorials <plot_2d>`. The code details are available in the
  :py:func:`PSyGrid.plot <posydon.grids.psygrid.PSyGrid.plot2D>` code and the
  :py:class:`visualization <posydon.visualization.plot2D>` library.
- :samp:`HR`: This is similar to :samp:`plot` but specialized for producing
  Hertzsprung–Russell diagrams.


Work on/with a `PSyGrid` object
-------------------------------

Loop over a `PSyGrid` object
~~~~~~~~~~~~~~~~~~~~~~~~~~

Similarly to accessing a single value in the :samp:`PSyGrid` object, we can loop
over a :samp:`PSyGrid` object, which will loop over the individual runs in the
:samp:`PSyGrid` object. Hence the following two code snippets will produce the same
output. The first one loops through the numpy array of the initial companion masses:

.. code-block:: python

    for mass in mygrid.initial_values['star_2_mass']:
        print(mass)

while the second one loops through the runs and prints the initial companion
mass of each one: 

.. code-block:: python

    for run in mygrid:
        print(run['initial_values']['star_2_mass'])



Expand a `PSyGrid` object
~~~~~~~~~~~~~~~~~~~~~~~~~

Because of the complexity of the :samp:`PSyGrid` object, we encourage users
to only use our dedicated functions to add content to the object. There is a
function to add an extra column to the :samp:`final_values`. Here is an example
of how to add a new column that contains the final orbital period, in units of
years instead of days:

.. code-block:: python

    new_column_data = mygrid.final_values['period_days']/365.25
    mygrid.add_column('period_years', new_column_data, where='final_values', overwrite=False)

The four arguments are a string with the name of the new field, the data to be
stored in the column, the component of the :samp:`PSyGrid` object to which the column will be
added, and a boolean indicating whether the column should overwrite any existing column with the same name.

.. warning::
    The new data has to have as many entries as the :samp:`PSyGrid` object has
    runs.

.. note::
    Currently, the parameter :samp:`where` only supports the value
    'final_values'.


Join two or more `PSyGrid` objects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are different reasons why you might have several :samp:`PSyGrid` objects
that you would like to combine into a single grid later, e.g. adding reruns.
POSYDON has a function called :samp:`join_grids` to do this for you. 
To avoid too many conflicts with possible modifications of already loaded 
:samp:`PSyGrid` objects, this function is not part of the :samp:`PSyGrid` 
object class. Instead, it takes as an argument the list of
paths to the hdf5 files containing :samp:`PSyGrid` objects to be combined to a
new one, and then a path to the hdf5 file of the new grid to be generated. 
The :samp:`join_grids` function will check whether the grids are compatible 
and join them if possible. Additionally, you can optionally specify the 
arguments :samp:`compression`, :samp:`description`, and :samp:`verbose`.

.. note::
    If there are common systems in two or more grids, this routine will only
    put the last run with same initial conditions in the newly combined
    :samp:`PSyGrid` object.

We recommend that you use the :ref:`post-processing pipeline <pipeline>` to create
and join grids.

..
    Extract the initial and final values as a pandas data frame
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TODO: get_pandas_initial_final()


Get reruns from a `PSyGrid` object
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Usually, not all runs of a grid will be successfully run in MESA. You
may want to rerun some of them with changed parameters. The function :samp:`rerun`
exports runs from a :samp:`PSyGrid` object to be run again. There are two options:

1. Write your own logic and create a numpy array with the indices of the systems that you would like to run again.
2. Specify which termination flag(s) necessitate a rerun of the system.

..
    .. table:: Arguments of the :samp:`rerun` function

    =================  =======  ===========
    Argument           Default  Description
    =================  =======  ===========
    path_to_file       './'     where to create the file(s) for the rerun
    runs_to_rerun      None     a numpy array containing the indices of the runs in the :samp:`PSyGrid` object (if given, leave :samp:`termination_flags=None`)
    termination_flags  None     a single termination flag code or a list of them (if given, leave :samp:`runs_to_rerun=None`)
    new_mesa_flag      None     dictionary with the names and the values of MESA parameters to be changed for the inlists of the new runs
    =================  =======  ===========

.. table:: Arguments of the :samp:`rerun` function

    =================  =======  ===========
    Argument           Default  Description
    =================  =======  ===========
    path_to_file       './'     where to create the file(s) for the rerun
    runs_to_rerun      None     a list containing the indices of the runs in the :samp:`PSyGrid` object
    termination_flags  None     a single termination flag code, or a list of them
    new_mesa_flag      None     dictionary with the names and the values of MESA parameters to be changed for the inlists of the new runs
    flags_to_check     None     a termination flag key or a list of them (if :samp:`None`, check only 'termination_flag_1')
    =================  =======  ===========

.. note::
    If both :samp:`runs_to_rerun` and :samp:`termination_flags` are given, all
    systems matching at least one of the two conditions will be selected for rerun.

The :ref:`post-processing pipeline <pipeline_stepR>` provides some pre-defined rerun options.


Close associated hdf5 file
~~~~~~~~~~~~~~~~~~~~~~~~~~

Finally, you can close the hdf5 file, which is recommended to ensure that all
your changes on the :samp:`PSyGrid` object are safely written into the file.

.. code-block:: python

    mygrid.close()

This is also done if you call the destructor of the
:samp:`PSyGrid` object.

.. code-block:: python

    del mygrid


The code summary of the :samp:`PSyGrid` object can be found at the
:py:class:`~posydon.grids.psygrid` reference page.
