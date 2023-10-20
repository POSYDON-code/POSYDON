.. _pipeline_steps:

##############
Pipeline steps
##############

The pipeline is devided into several steps, which build up on eachother. Each
step will take a csv file as input. The name of this file is used to tell the
pipeline, which step should be performed.

The run-pipeline script takes four arguments:

.. code-block:: bash

    run-pipeline PATH_TO_GRIDS PATH_TO_CSV_FILE DATA_ID VERBOSE

1. [path] The path to the girds main directory (currently not used)
2. [path] The path to the csv file
3. [int] An index indicating the data entry to read from the csv file
4. [int] Where one wants verbose output (1) or not (0)

.. note::
    The current directory will be used as working directory, hence navigate to
    your work directory first.

Step1: creating a `psygrid` object
----------------------------------

Frist, we need to create the `psygird` object. To do so, the pipeline needs to
now the directory which contains the MESA runs, the compression, and whether to
crop the history for some certain runs. Hence, the :samp:`step_1.csv` file
should have those columns:

.. code-block::

    path_to_grid,compression,stop_before_carbon_depletion

And the lines below contain the data for each unique combination of the three
parameters to be processed. Here the :samp:`DATA_ID` simply refers to the line
below the header starting by 0. Thus, the second line in the file has the index
0, the thrid one has index 1 and so on.

The currently supported compression types are:

.. table:: Compression types

    ========  ===========
    Type      Description
    ========  ===========
    ORIGINAL  It keeps all columns and history entries given by MESA
    LITE      It discards some columns and reduces the history and final profiles to an maximum error of 0.1 and limit profiles to contain in maximum 200 data points
    \*_RLO    The `_RLO` can be put on any of the previous types to crop the history before the onset of Roche-lobe overflow for grids containing a compact object
    ========  ===========

Step2: combining `psygrid` objects
----------------------------------

Usually, the girds are slitt into batches or reruns are done. In those cases,
there will be several `psygrid` objects created for one gird. This step will
join them into one. The :samp:`step_2.csv` file should have a matrix structure.
The columns contain the girds which should be combined to the one specified in
the header (first) row. The :samp:`DATA_ID` corresponds here to the column
number (starting with 0). Here an example:

.. code-block::

    NEW_H5_FILE1,NEW_H5_FILE2
    OLD_H5_FILE11,OLD_H5_FILE21
    OLD_H5_FILE12,OLD_H5_FILE22
    ,OLD_H5_FILE23

.. warning::
    The data will be put on top of eachother. E.g. if there is the same initial
    system in :samp:`OLD_H5_FILE11` and :samp:`OLD_H5_FILE12`, the one in
    :samp:`OLD_H5_FILE11` will be discarded and only the one in
    :samp:`OLD_H5_FILE12` will end up in :samp:`NEW_H5_FILE1`.

Step3: calculating extra values from detailed data
--------------------------------------------------

In this step we calculate extra quanties from the final profiles, which are
important for the core collapse of the star to a compact object.

TBC

Step4: training of the interpolators
------------------------------------

TODO

Step9: exporting the data set
-----------------------------

TODO

StepR: exporting a rerun
------------------------

TODO
