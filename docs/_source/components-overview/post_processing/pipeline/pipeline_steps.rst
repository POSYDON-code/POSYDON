.. _pipeline_steps:

##############
Pipeline steps
##############
The post-processing pipeline is divided into several steps which build 
on each other. Each step will take a :samp:`csv` file as input. The name of this 
file determines which pipeline step should be performed.

The script to run the pipeline takes four arguments:

.. code-block:: bash

    posydon-run-pipeline PATH_TO_GRIDS PATH_TO_CSV_FILE DATA_ID VERBOSE

1. [path] The path to the grids main directory (currently not used)
2. [path] The path to the :samp:`csv` file
3. [int] An index indicating the data entry to read from the :samp:`csv` file
4. [int] Whether one wants verbose output (1) or not (0)

.. note::
    The current directory will be used as working directory, so navigate to
    the correct location first.

.. _pipeline_step1:

Step 1: Creating a `PSyGrid` object
----------------------------------

First, we need to create the :samp:`PSygrid` object. To do so, the pipeline
needs to know the directory which contains the MESA runs, the compression, the
grid type, and whether to crop the history for some certain runs. Hence, the
:samp:`step_1.csv` file should have the following columns:

.. code-block::

    path_to_grid,compression,grid_type,stop_before_carbon_depletion

And the lines below contain the data for each unique combination of the three
parameters to be processed. Here the :samp:`DATA_ID` simply refers to the line
below the header starting by 0. Thus, the second line in the file has the index
0, the third one has index 1, and so on.

The currently supported compression types are:

.. table:: Compression types

    ========  ===========
    Type      Description
    ========  ===========
    ORIGINAL  Keeps all columns and history entries given by MESA
    LITE      Discards some columns and reduces the history and final profiles to an maximum error of 0.1, and limits profiles to contain a maximum of 200 data points
    \*_RLO    The :samp:`_RLO` can be put on any of the previous types to crop the history before the onset of Roche-lobe overflow for grids containing a compact object
    ========  ===========

.. _pipeline_step2:

Step 2: Combining `PSyGrid` objects
----------------------------------

Often, grids are split into batches, or reruns are done. In those cases,
there will be individual :samp:`PSyGrid` objects created for each batch or rerun. This step
will join them into a single combined :samp:`PSyGrid` object representing the complete grid. The :samp:`step_2.csv` file should have a matrix
structure. The columns contain the grids which should be combined to the one
specified in the header (first) row. The :samp:`DATA_ID` corresponds here to
the column number (starting with 0). Here is an example:

.. code-block::

    NEW_H5_FILE1,NEW_H5_FILE2
    OLD_H5_FILE11,OLD_H5_FILE21
    OLD_H5_FILE12,OLD_H5_FILE22
    ,OLD_H5_FILE23

.. warning::
    The data will be layered on top of each other. E.g., if there is the same
    initial system in :samp:`OLD_H5_FILE11` and :samp:`OLD_H5_FILE12`, the one
    in :samp:`OLD_H5_FILE11` will be discarded and only the one in
    :samp:`OLD_H5_FILE12` will end up in :samp:`NEW_H5_FILE1`.

.. _pipeline_step3:

Step 3: Calculating extra values from stellar model time series and structure profile data
--------------------------------------------------

In this step we calculate extra quantities from the histories and profiles.
Those extra values are key parameters at He depletion, at onset of common
envelope evolution, and at core collapse.

Because some of the values may require a high precision in the data, we
recommend to use the data from the ORIGINAL compression to calculate them. But
the new values can be added to any :samp:`PSyGrid` object. Hence this step
requests three paths to be specified in :samp:`step_3.csv` besides the grid
type:

.. code-block::

    path_to_grid,grid_type,path_to_grid_ORIGINAL,path_to_processed_grid

.. table:: Description of required paths

    ======================  ===========
    Path                    Description
    ======================  ===========
    path_to_grid            path of the grid, which gets the values appended to it
    grid_type               type of the grid
    path_to_grid_ORIGINAL   path of the grid, where the values are calculated from
    path_to_processed_grid  path of the new grid (a copy of the one specified as :samp:`path_to_grid` with the appended values)
    ======================  ===========

.. note::
    This step use the path to the original MESA data as the unique identifier
    of each system in the :samp:`PSyGrid` object, thus the location of the MESA
    file cannot be changed between creating two :samp:`PSyGrid` objects of the
    same grid in :ref:`step1 <pipeline_step1>`. Similarly, the overlaying in
    :ref:`step2 <pipeline_step2>` needs to be the same, too. Therefore, we
    recommend to setup and run the pipeline with an
    :ref:`ini file <pipeline_ini>`.

.. _pipeline_step4:

Step 4: Training the interpolators
------------------------------------

To get interpolated data from our grids, in this step we train an interpolator
on the :samp:`PSyGrid` object. The file :samp:`step_4.csv` therefore has to
contain the following pieces of information: First, the grid containing the 
data, second, the grid type, third, the interpolation method (inlcuding whether 
the grid starts at RLO), and finally, the name of the interpolator object.

.. code-block::

    path_to_grid,interpolation_method,path_to_interpolator

.. note::
    The type of interpolator will be recognized from the name of the
    interpolator object. The syntax is :code:`IF_METHOD{_RLO}.pkl`. The
    :samp:`IF` stands for initial-final interpolator, the :samp:`METHOD` refers
    to the interpolator type. The grids starting at Roche-lobe overflow may be
    indicated in the name as well, but is not required.

.. table:: Currently supported interpolator types

    ==============  ===========
    :samp:`METHOD`  Description
    ==============  ===========
    linear          linear interpolation
    1NN             nearest neighbor
    ==============  ===========

.. _pipeline_stepF:

Step F: Exporting the data set
-----------------------------

After we have a complete data set, we would like to export it to be used for
the population synthesis. We now go to the final step, step F.
the final step even if more steps are introduced in the future. In
:samp:`step_F.csv`, there are again two paths required, a source and an export
path. The step will simply copy the source to the export location. Hence, here
the final :samp:`PSyGrid` objects and all the interpolator files are usually
addressed by this step.

.. code-block::

    path_to_grid,export_path

.. _pipeline_stepR:

Step R: Exporting a rerun
------------------------

Often, a grid will not successfully converge every binary on the first go. So 
we may need to export reruns which use modified conditions to fix 
non-converged models. This step is therefore only needed to build 
a new grid. Usually, one would run the steps to the point where the need of a 
fix arises. Additionally, before exporting a rerun, the logic for how to select 
a system to be included in the rerun and what should be changed needs to be 
implemented first.

For this step the :samp:`csv` file is called :samp:`rerun.csv` to avoid too much
confusion with other steps. It clearly has to run after step 1 and step 2, but it is 
not a usual step itself. It requires the path to a :samp:`PSyGrid` object to 
get the non-converged models from, the path to which the rerun should be stored (it creates 
the :samp:`grid.csv` and the :samp:`ini` file needed to
:ref:`setup a new run <mesa-grids-api>`), the grid type, the metallicity, the
type of the rerun specifying the logic and changes, and the cluster name.

.. code-block::

    path_to_grid,rerun_path,grid_type,rerun_metallicity,rerun_type,cluster

.. table:: Currently supported rerun types

    =====================  ==============  ===========
    :samp:`rerun_type`     Future version  Description
    =====================  ==============  ===========
    PISN                   default in v3+  Enables the MESA inlist commit, which stops MESA before getting dynamical to save a final profile there
    reverse_MT             default in v3+  Uses a MESA version with a bug fix, that the role of donor and accretor can switch during the simulation
    opacity_max            caution         Uses a fixed maximum opacity of 0.5 (this is only a last option change to get more stability)
    TPAGBwind              default in v3+  Enables the MESA inlist commit, which changes the wind during the TPAGB phase
    thermohaline_mixing    default in v3+  Uses thermohaline mixing in the inlist
    HeMB_MLTp_mesh         caution         Turns off magnetic braking for He stars; it uses less extreme parameters of the MLT++ (this can cause significant changes in the radius evolution of stars); it changes some more input values to change the resulation close to the surface
    more_mesh              workaround      Modifies the remeshing and allows for more cells in MESA
    conv_bdy_weight        caution         Disables the convective_bdy_weight where this caused segmentation faults (this avoids a bug in the old MESA version r11701)
    dedt_energy_eqn        caution         Enables MESA's dedt-form of the energy equation for numerical stability during rapid (superthermal) mass transfer
    dedt_hepulse           caution         Enables MESA's dedt-form of the energy equation for rapid mass transfer; at stripped HeZAMS, several MLT++ changes, v_flag and lnPgas_flag set to .true., and convective_bdy_weight disabled to help with stripped He star superadiabatic envelopes, pulsations, and WD cooling
    LBV_wind               default in v3+  Turns on LBV winds when crossing the Humphreys-Davidson limit as intended (due to a bug this was only applied after a retry); additionally, there are reruns `LBV_wind+thermohaline_mixing`, `LBV_wind+dedt_energy_eqn`, which combine the two rerun types. Any additional changes to these reruns are described here as LBV_wind+rerun_type
    no_age_limit           default in v3+  Allows low mass stars to evolve beyond the age of the universe, which is needed for grids where we jump on past ZAMS; additionally, there are reruns `no_age_limit+thermohaline_mixing` and `no_age_limit+dedt_energy_eqn`, which combine the two rerun types
    LBV_wind+dedt          caution         Enables MESA's dedt-form of the energy equation for numerical stability during rapid (superthermal) mass transfer and sets lnPgas_flag to .true. for numerical stability. Also disabled convective_bdy_weight as a degenerate core is forming (as probed by the central Coulomb coupling parameter) to avoid segmentation faults.
    LBV_wind+hepulse       caution         Contains the LBV_wind+dedt_energy_rerun; additionally, at stripped HeZAMS, the thresholds to trigger MLT++ are relaxed, and several timestep controls limiting the allowed variation of lgTeff and (cell-wise) T, as well as controls limiting the allowed variation of donor envelope mass are relaxed during mass transfer to improve convergence during envelope stripping. Also removes stopping conditions for Hubble time and TAMS that would be enforced for models less massive than roughly G-type stars, relevant to single_* and CO_* grids.

    =====================  ==============  ===========

