.. _pipeline_ini:

##############################################################
Documentation of the ini file for the post-processing pipeline
##############################################################

Aim of the ini file
===================

Using an ini file should help to keep an overview on large grid repositories
and ensures that all will be setup the same way.

There is a setup-pipeline script takes one argument:

.. code-block:: bash

    setup-pipeline PATH_TO_INI

The content of the ini file is desciped :ref:`below <pipeline_ini_sections>`.
It will create for each step, plot or check two files:

1. \*.csv
2. \*.slurm

Additionally, it will prepare a :samp:`logs` directory (it may creates a save
of some old logs) and create a shell script :samp:`run_pipeline.sh` to submit
all :ref:`tasks <pipeline>`. Hence, you can run all
:ref:`steps <pipeline_steps>` with simply running this script.

.. code-block:: bash

    ./run_pipeline.sh

.. note::
    Currently, the user needs to take care of having a POSYDON_data directly
    which includes the tables for the core-collape prescriptions him-/herself
    in the working directory.

.. _pipeline_ini_sections:

Sections in the ini file
========================

Account and slurm settings
--------------------------

First we specify some details for the slurm jobs. This is similar to the
:ref:`[slurm] section for running the grids <inifile_slurm>`. Here is an
example to show the supported key words:

.. code-block:: ini

    [account]
        ACCOUNT = 'fragkos'
        PARTITION = 'public-cpu'
        WALLTIME = '23:00:00'
        MAILTYPE = 'ALL'
        EMAIL = 'matthias.kruckow@unige.ch'

General pipeline settings
-------------------------

The next sections deals with the general inforamtion about the pipeline. First
it needs to know where the grids are located. The `PATH` specifies, where you
would like to get the pipeline files being created. The `VERBOSE` option will
be used for the creation of the pipeline files and during the run of the
pipeline.

Finially, we have switches to turn on (`True`) and off (`False`) individual
:ref:`steps <pipeline_steps>`.

.. code-block:: ini

    [pipeline setup]
        PATH_TO_GRIDS = '/srv/beegfs/scratch/shares/astro/posydon/POSYDON_GRIDS_v2/'
        VERSION = '' # 'v2' in quest and '' in yggdrasil
        PATH = '.' # working dir
        VERBOSE = True
        
        # steps
        CREATE_GRID_SLICES = True
        COMBINE_GRID_SLICES = True
        CALCULATE_EXTRA_VALUES = True
        TRAIN_INTERPOLATORS = True
        EXPORT_DATASET = True
        # rerun step
        RERUN = False

Step sections
-------------

The path of each grid will be joint as
:samp:`PATH_TO_GRIDS/VERSION/GRID_TYPE/METALLICITY/GRID_SLICE`. The
corresponding h5 files will have names according to
:samp:`PATH_TO_GRIDS/VERSION/GRID_TYPE/METALLICITY/COMPRESSION/GRID_SLICE.h5`.
All sections have common keywords:

.. table:: Common keywords of steps

    ==================  ===========
    Keyword             Description
    ==================  ===========
    GRID_TYPES          a list of grid types; the looped :samp:`GRID_TYPE` is used in the path name
    METALLICITIES       a list of lists of the metallicities of the grids; the looped :samp:`METALLICITY` is used in the path name; the outer list allows you to have different lists for each grid type
    GRID_SLICES         a list of lists of the grid slices; the looped :samp:`GRID_SLICE` is used in the path name; the outer list allows you to have different lists for each grid type
    COMPRESSIONS        a list of lists of compression types
    DROP_MISSING_FILES  boolean to igrnore missing files
    CREATE_PLOTS        a list of plots to make; this will be done independently whether the step is active or not, to make no plots put there an empty list or comment out such a line
    DO_CHECKS           a list of checks to perform; this will be done independently whether the step is active or not, to make no checks put there an empty list or comment out such a line
    ==================  ===========

Some :ref:`steps <pipeline_steps>` have more keywords, which are specific to
that step:

.. table:: Step specific keywords

    ====  ============================  ===========
    Step  Keyword                       Description
    ====  ============================  ===========
       1  STOP_BEFORE_CARBON_DEPLETION  indicating, whether high mass HMS stars should get their history croped short before carbon depletion (1) or not (0)
       2  GRID_SLICES                   for this step, we have 3 layers of lists: the outermost is still the grid type, the inner most is still the grid slice, the middle layer is the combined grid
       2  GRIDS_COMBINED                a list of lists of combined grids; the outermost list is again refering to grid type; this is used as name for the new combined grid instead of `GRID_SLICE`
       4  INTERPOLATION_METHODS         a list of the interpolator types which are trained
       4  CONTROL_GRIDS                 a list of lists of control grids for the `GRID_SLICES`; it need to have the same number of entries as the `GRID_SLICES`, to specify no control grid use an empty string
       R  RERUN_TYPE                    a defined rerun type
    ====  ============================  ===========

Here is an example of all the :ref:`steps <pipeline_steps>`:

.. code-block:: ini

    #CREATE_GRID_SLICES
    [step_1]
        # e.g. ['CO-HMS_RLO','CO-HeMS','HMS-HMS']
        GRID_TYPES = ['CO-HMS_RLO', 'CO-HeMS', 'HMS-HMS']
        # e.g. ['2e+00_Zsun', '1e+00_Zsun', '4.5e-01_Zsun', '2e-01_Zsun', '1e-01_Zsun', '1e-02_Zsun', '1e-03_Zsun', '1e-04_Zsun']
        METALLICITIES = [# CO-HMS_RLO
                         ['2e+00_Zsun', '1e+00_Zsun', '4.5e-01_Zsun', '2e-01_Zsun', '1e-01_Zsun', '1e-02_Zsun', '1e-03_Zsun', '1e-04_Zsun'],
                         # CO-HeMS
                         ['2e+00_Zsun', '1e+00_Zsun', '4.5e-01_Zsun', '2e-01_Zsun', '1e-01_Zsun', '1e-02_Zsun', '1e-03_Zsun', '1e-04_Zsun'],
                         # HMS-HMS
                         ['2e+00_Zsun', '1e+00_Zsun', '4.5e-01_Zsun', '2e-01_Zsun', '1e-01_Zsun', '1e-02_Zsun', '1e-03_Zsun', '1e-04_Zsun']
                        ]
        GRID_SLICES = [# CO-HMS_RLO
                       ['grid_low_res_0', 'grid_low_res_1', 'grid_low_res_2', 'rerun_PISN_grid_low_res_combined', 'rerun_TPAGBwind_grid_low_res_combined',
                        'grid_random_1', 'rerun_PISN_grid_random_combined', 'rerun_TPAGBwind_grid_random_combined'],
                       # CO-HeMS
                       ['grid_low_res_0', 'grid_low_res_1', 'grid_low_res_2', 'rerun_PISN_grid_low_res_combined',
                        'grid_random_1', 'grid_random_rerun', 'rerun_PISN_grid_random_combined'],
                       # HMS-HMS
                       ['grid_low_res_0', 'grid_low_res_1', 'grid_low_res_2', 'grid_low_res_3', 'grid_low_res_4', 'grid_low_res_5', 'rerun_PISN_grid_low_res_combined', 'rerun_reverse_MT_grid_low_res_combined', 'rerun_TPAGBwind_grid_low_res_combined',
                        'grid_random_1', 'rerun_PISN_grid_random_combined', 'rerun_reverse_MT_grid_random_combined', 'rerun_TPAGBwind_grid_random_combined']
                      ]
        COMPRESSIONS = [# CO-HMS_RLO
                        ['LITE', 'ORIGINAL', 'LITE_RLO', 'ORIGINAL_RLO'],
                        # CO-HeMS
                        ['LITE', 'ORIGINAL', 'LITE_RLO', 'ORIGINAL_RLO'],
                        # HMS-HMS
                        ['LITE', 'ORIGINAL']
                       ]
        DROP_MISSING_FILES = True
        # EXTRA PARAMETERS
        # only applied to HMS grids
        STOP_BEFORE_CARBON_DEPLETION = 1
        # supported plots: e.g. 'combined_TF12', 'termination_flag_1', 'termination_flag_2', 'termination_flag_3', 'termination_flag_4', and any quantity valid for a Z-plotting
        CREATE_PLOTS = []
        # supported checks: e.g. 'failure_rate'
        DO_CHECKS = []
    
    #COMBINE_GRID_SLICES
    [step_2]
        GRID_TYPES = ['CO-HMS_RLO', 'CO-HeMS', 'HMS-HMS']
        METALLICITIES = [# CO-HMS_RLO
                         ['2e+00_Zsun', '1e+00_Zsun', '4.5e-01_Zsun', '2e-01_Zsun', '1e-01_Zsun', '1e-02_Zsun', '1e-03_Zsun', '1e-04_Zsun'],
                         # CO-HeMS
                         ['2e+00_Zsun', '1e+00_Zsun', '4.5e-01_Zsun', '2e-01_Zsun', '1e-01_Zsun', '1e-02_Zsun', '1e-03_Zsun', '1e-04_Zsun'],
                         # HMS-HMS
                         ['2e+00_Zsun', '1e+00_Zsun', '4.5e-01_Zsun', '2e-01_Zsun', '1e-01_Zsun', '1e-02_Zsun', '1e-03_Zsun', '1e-04_Zsun']
                        ]
        GRID_SLICES = [# CO-HMS_RLO
                       [['grid_low_res_0', 'grid_low_res_1', 'grid_low_res_2'],
                        ['grid_low_res_0', 'grid_low_res_1', 'grid_low_res_2', 'rerun_PISN_grid_low_res_combined'],
                        ['grid_low_res_0', 'grid_low_res_1', 'grid_low_res_2', 'rerun_PISN_grid_low_res_combined', 'rerun_TPAGBwind_grid_low_res_combined'],
                        ['grid_random_1'],
                        ['grid_random_1', 'rerun_PISN_grid_random_combined'],
                        ['grid_random_1', 'rerun_PISN_grid_random_combined', 'rerun_TPAGBwind_grid_random_combined']],
                       # CO-HeMS
                       [['grid_low_res_0', 'grid_low_res_1', 'grid_low_res_2'],
                        ['grid_low_res_0', 'grid_low_res_1', 'grid_low_res_2', 'rerun_PISN_grid_low_res_combined'],
                        ['grid_random_1', 'grid_random_rerun'],
                        ['grid_random_1', 'grid_random_rerun', 'rerun_PISN_grid_random_combined']],
                       # HMS-HMS
                       [['grid_low_res_0', 'grid_low_res_1', 'grid_low_res_2', 'grid_low_res_3', 'grid_low_res_4', 'grid_low_res_5'],
                        ['grid_low_res_0', 'grid_low_res_1', 'grid_low_res_2', 'grid_low_res_3', 'grid_low_res_4', 'grid_low_res_5', 'rerun_PISN_grid_low_res_combined'],
                        ['grid_low_res_0', 'grid_low_res_1', 'grid_low_res_2', 'grid_low_res_3', 'grid_low_res_4', 'grid_low_res_5', 'rerun_PISN_grid_low_res_combined', 'rerun_reverse_MT_grid_low_res_combined'],
                        ['grid_low_res_0', 'grid_low_res_1', 'grid_low_res_2', 'grid_low_res_3', 'grid_low_res_4', 'grid_low_res_5', 'rerun_PISN_grid_low_res_combined', 'rerun_reverse_MT_grid_low_res_combined', 'rerun_TPAGBwind_grid_low_res_combined'],
                        ['grid_random_1'],
                        ['grid_random_1', 'rerun_PISN_grid_random_combined'],
                        ['grid_random_1', 'rerun_PISN_grid_random_combined', 'rerun_reverse_MT_grid_random_combined'],
                        ['grid_random_1', 'rerun_PISN_grid_random_combined', 'rerun_reverse_MT_grid_random_combined', 'rerun_TPAGBwind_grid_random_combined']]
                      ]
        GRIDS_COMBINED = [# CO-HMS_RLO
                          ['grid_low_res_combined', 'grid_low_res_combined_rerun1_PISN', 'grid_low_res_combined_rerun3_TPAGBwind',
                           'grid_random_combined', 'grid_random_combined_rerun1_PISN', 'grid_random_combined_rerun3_TPAGBwind'],
                          # CO-HeMS
                          ['grid_low_res_combined', 'grid_low_res_combined_rerun1_PISN',
                           'grid_random_combined', 'grid_random_combined_rerun1_PISN'],
                          # HMS-HMS
                          ['grid_low_res_combined', 'grid_low_res_combined_rerun1_PISN', 'grid_low_res_combined_rerun2_reverse_MT', 'grid_low_res_combined_rerun3_TPAGBwind',
                           'grid_random_combined', 'grid_random_combined_rerun1_PISN', 'grid_random_combined_rerun2_reverse_MT', 'grid_random_combined_rerun3_TPAGBwind']
                         ]
        COMPRESSIONS = [# CO-HMS_RLO
                        ['LITE', 'ORIGINAL', 'LITE_RLO', 'ORIGINAL_RLO'],
                        # CO-HeMS
                        ['LITE', 'ORIGINAL', 'LITE_RLO', 'ORIGINAL_RLO'],
                        # HMS-HMS
                        ['LITE', 'ORIGINAL']
                       ]
        DROP_MISSING_FILES = True
        # supported plots: e.g. 'combined_TF12', 'termination_flag_1', 'termination_flag_2', 'termination_flag_3', 'termination_flag_4', and any quantity valid for a Z-plotting
        CREATE_PLOTS = ['PLOT_AFTER_COMBINE']
        # supported checks: e.g. 'failure_rate'
        DO_CHECKS = ['CHECK_AFTER_COMBINE']
    
    #CALCULATE_EXTRA_VALUES
    [step_3]
        GRID_TYPES = ['CO-HMS_RLO', 'CO-HeMS', 'HMS-HMS']
        METALLICITIES = [# CO-HMS_RLO
                         ['2e+00_Zsun', '1e+00_Zsun', '4.5e-01_Zsun', '2e-01_Zsun', '1e-01_Zsun', '1e-02_Zsun', '1e-03_Zsun', '1e-04_Zsun'],
                         # CO-HeMS
                         ['2e+00_Zsun', '1e+00_Zsun', '4.5e-01_Zsun', '2e-01_Zsun', '1e-01_Zsun', '1e-02_Zsun', '1e-03_Zsun', '1e-04_Zsun'],
                         # HMS-HMS
                         ['2e+00_Zsun', '1e+00_Zsun', '4.5e-01_Zsun', '2e-01_Zsun', '1e-01_Zsun', '1e-02_Zsun', '1e-03_Zsun', '1e-04_Zsun']
                        ]
        GRID_SLICES = [# CO-HMS_RLO
                       ['grid_low_res_combined_rerun3_TPAGBwind', 'grid_random_combined_rerun3_TPAGBwind'],
                       # CO-HeMS
                       ['grid_low_res_combined_rerun1_PISN', 'grid_random_combined_rerun1_PISN'],
                       # HMS-HMS
                       ['grid_low_res_combined_rerun3_TPAGBwind', 'grid_random_combined_rerun3_TPAGBwind']
                      ]
        COMPRESSIONS = [# CO-HMS_RLO
                        ['LITE', 'LITE_RLO'],
                        # CO-HeMS
                        ['LITE', 'LITE_RLO'],
                        # HMS-HMS
                        ['LITE']
                       ]
        DROP_MISSING_FILES = True
        # supported plots: e.g. 'combined_TF12', 'termination_flag_1', 'termination_flag_2', 'termination_flag_3', 'termination_flag_4', and any quantity valid for a Z-plotting
        CREATE_PLOTS = ['PLOT_AFTER_EXTRA']
        # supported checks: e.g. 'failure_rate', 'CO_TYPE', 'SN_TYPE'
        DO_CHECKS = ['CHECK_AFTER_EXTRA']
    
    #TRAIN_INTERPOLATORS
    [step_4]
        GRID_TYPES = ['CO-HMS_RLO', 'CO-HeMS', 'HMS-HMS']
        METALLICITIES = [# CO-HMS_RLO
                         ['2e+00_Zsun', '1e+00_Zsun', '4.5e-01_Zsun', '2e-01_Zsun', '1e-01_Zsun', '1e-02_Zsun', '1e-03_Zsun', '1e-04_Zsun'],
                         # CO-HeMS
                         ['2e+00_Zsun', '1e+00_Zsun', '4.5e-01_Zsun', '2e-01_Zsun', '1e-01_Zsun', '1e-02_Zsun', '1e-03_Zsun', '1e-04_Zsun'],
                         # HMS-HMS
                         ['2e+00_Zsun', '1e+00_Zsun', '4.5e-01_Zsun', '2e-01_Zsun', '1e-01_Zsun', '1e-02_Zsun', '1e-03_Zsun', '1e-04_Zsun']
                        ]
        GRID_SLICES = [# CO-HMS_RLO
                       ['grid_low_res_combined_rerun3_TPAGBwind_processed'],
                       # CO-HeMS
                       ['grid_low_res_combined_rerun1_PISN_processed'],
                       # HMS-HMS
                       ['grid_low_res_combined_rerun3_TPAGBwind_processed']
                      ]
        INTERPOLATION_METHODS = ["linear","1NN"]
        COMPRESSIONS = [# CO-HMS_RLO
                        ['LITE_RLO'],
                        # CO-HeMS
                        ['LITE', 'LITE_RLO'],
                        # HMS-HMS
                        ['LITE']
                       ]
        CONTROL_GRIDS = [# CO-HMS_RLO
                         ['grid_random_combined_rerun3_TPAGBwind_processed'],
                         # CO-HeMS
                         ['grid_random_combined_rerun1_PISN_processed'],
                         # HMS-HMS
                         ['grid_random_combined_rerun3_TPAGBwind_processed']
                        ]
        DROP_MISSING_FILES = True
        # supported plots: e.g. 'combined_TF12', 'termination_flag_1', 'termination_flag_2', 'termination_flag_3', 'termination_flag_4', and any quantity valid for a Z-plotting
        CREATE_PLOTS = ['PLOT_AFTER_TRAINING']
        # supported checks: e.g. 'failure_rate'
        DO_CHECKS = ['CHECK_AFTER_TRAINING']
    
    #EXPORT_DATASET
    [step_9]
        GRID_TYPES = ['CO-HMS_RLO', 'CO-HeMS', 'HMS-HMS']
        METALLICITIES = [# CO-HMS_RLO
                         ['2e+00_Zsun', '1e+00_Zsun', '4.5e-01_Zsun', '2e-01_Zsun', '1e-01_Zsun', '1e-02_Zsun', '1e-03_Zsun', '1e-04_Zsun'],
                         # CO-HeMS
                         ['2e+00_Zsun', '1e+00_Zsun', '4.5e-01_Zsun', '2e-01_Zsun', '1e-01_Zsun', '1e-02_Zsun', '1e-03_Zsun', '1e-04_Zsun'],
                         # HMS-HMS
                         ['2e+00_Zsun', '1e+00_Zsun', '4.5e-01_Zsun', '2e-01_Zsun', '1e-01_Zsun', '1e-02_Zsun', '1e-03_Zsun', '1e-04_Zsun']
                        ]
        GRID_SLICES = [# CO-HMS_RLO
                       ['grid_low_res_combined_rerun3_TPAGBwind_processed'],
                       # CO-HeMS
                       ['grid_low_res_combined_rerun1_PISN_processed'],
                       # HMS-HMS
                       ['grid_low_res_combined_rerun3_TPAGBwind_processed']
                      ]
        COMPRESSIONS = [# CO-HMS_RLO
                        ['LITE_RLO'],
                        # CO-HeMS
                        ['LITE', 'LITE_RLO'],
                        # HMS-HMS
                        ['LITE']
                       ]
        DROP_MISSING_FILES = True
    
    #EXPORT_RERUNS
    [rerun]
        GRID_TYPES = ['CO-HMS_RLO', 'CO-HeMS', 'HMS-HMS']
        METALLICITIES = [# CO-HMS_RLO
                         ['2e+00_Zsun', '1e+00_Zsun', '4.5e-01_Zsun', '2e-01_Zsun', '1e-01_Zsun', '1e-02_Zsun', '1e-03_Zsun', '1e-04_Zsun'],
                         # CO-HeMS
                         ['2e+00_Zsun', '1e+00_Zsun', '4.5e-01_Zsun', '2e-01_Zsun', '1e-01_Zsun', '1e-02_Zsun', '1e-03_Zsun', '1e-04_Zsun'],
                         # HMS-HMS
                         ['2e+00_Zsun', '1e+00_Zsun', '4.5e-01_Zsun', '2e-01_Zsun', '1e-01_Zsun', '1e-02_Zsun', '1e-03_Zsun', '1e-04_Zsun']
                        ]
        GRID_SLICES = [# CO-HMS_RLO
                       ['grid_low_res_combined_rerun3_TPAGBwind','grid_random_combined_rerun3_TPAGBwind'],
                       # CO-HeMS
                       ['grid_low_res_combined_rerun1_PISN','grid_random_combined_rerun1_PISN'],
                       # HMS-HMS
                       ['grid_low_res_combined_rerun3_TPAGBwind','grid_random_combined_rerun3_TPAGBwind']
                      ]
        COMPRESSIONS = [# CO-HMS_RLO
                        ['LITE'],
                        # CO-HeMS
                        ['LITE'],
                        # HMS-HMS
                        ['LITE']
                       ]
        DROP_MISSING_FILES = True
        # example reruns are 'PISN', 'reverse_MT', 'TPAGBwind', 'opacity_max'
        RERUN_TYPE = 'opacity_max' 

There are some predefined shortcuts for lists of :ref:`plots <pipeline_plots>`
and :ref:`checks <pipeline_checks>`:

.. table:: Plot sets

    =====================  =====
    Set name               plots
    =====================  =====
    'PLOT_AFTER_CREATE'    
    'PLOT_AFTER_COMBINE'   'combined_TF12', 'termination_flag_1', 'termination_flag_2', 'termination_flag_3', 'termination_flag_4', 'rl_relative_overflow_1', 'rl_relative_overflow_2', 'lg_mtransfer_rate'
    'PLOT_AFTER_EXTRA'     'S1_MODEL01_CO_type', 'S1_MODEL01_SN_type', 'S1_MODEL01_mass', 'S1_MODEL01_spin', 'S1_MODEL01_m_disk_radiated', 'S1_MODEL02_CO_type', 'S1_MODEL02_SN_type', 'S1_MODEL02_mass', 'S1_MODEL02_spin', 'S1_MODEL02_m_disk_radiated', 'S1_MODEL03_CO_type', 'S1_MODEL03_SN_type', 'S1_MODEL03_mass', 'S1_MODEL03_spin', 'S1_MODEL03_m_disk_radiated', 'S1_MODEL04_CO_type', 'S1_MODEL04_SN_type', 'S1_MODEL04_mass', 'S1_MODEL04_spin', 'S1_MODEL04_m_disk_radiated', 'S1_MODEL05_CO_type', 'S1_MODEL05_SN_type', 'S1_MODEL05_mass', 'S1_MODEL05_spin', 'S1_MODEL05_m_disk_radiated', 'S1_MODEL06_CO_type', 'S1_MODEL06_SN_type', 'S1_MODEL06_mass', 'S1_MODEL06_spin', 'S1_MODEL06_m_disk_radiated', 'S1_MODEL07_CO_type', 'S1_MODEL07_SN_type', 'S1_MODEL07_mass', 'S1_MODEL07_spin', 'S1_MODEL07_m_disk_radiated', 'S1_MODEL08_CO_type', 'S1_MODEL08_SN_type', 'S1_MODEL08_mass', 'S1_MODEL08_spin', 'S1_MODEL08_m_disk_radiated', 'S1_MODEL09_CO_type', 'S1_MODEL09_SN_type', 'S1_MODEL09_mass', 'S1_MODEL09_spin', 'S1_MODEL09_m_disk_radiated', 'S1_MODEL10_CO_type', 'S1_MODEL10_SN_type', 'S1_MODEL10_mass', 'S1_MODEL10_spin', 'S1_MODEL10_m_disk_radiated'
    'PLOT_AFTER_TRAINING'  'INTERP_ERROR_age', 'INTERP_ERROR_star_1_mass', 'INTERP_ERROR_star_2_mass', 'INTERP_ERROR_period_days', 'INTERP_ERROR_S1_co_core_mass', 'INTERP_ERROR_S1_co_core_radius', 'INTERP_ERROR_S1_he_core_mass', 'INTERP_ERROR_S1_he_core_radius', 'INTERP_ERROR_S1_center_h1', 'INTERP_ERROR_S1_center_he4', 'INTERP_ERROR_S1_surface_h1', 'INTERP_ERROR_S1_surface_he4', 'INTERP_ERROR_S1_surf_avg_omega_div_omega_crit', 'INTERP_ERROR_S1_log_Teff', 'INTERP_ERROR_S1_log_L', 'INTERP_ERROR_S1_log_R', 'INTERP_ERROR_S1_spin_parameter', 'INTERP_ERROR_S1_lambda_CE_10cent', 'INTERP_ERROR_S2_co_core_mass', 'INTERP_ERROR_S2_co_core_radius', 'INTERP_ERROR_S2_he_core_mass', 'INTERP_ERROR_S2_he_core_radius', 'INTERP_ERROR_S2_center_h1', 'INTERP_ERROR_S2_center_he4', 'INTERP_ERROR_S2_surface_h1', 'INTERP_ERROR_S2_surface_he4', 'INTERP_ERROR_S2_surf_avg_omega_div_omega_crit', 'INTERP_ERROR_S2_log_Teff', 'INTERP_ERROR_S2_log_L', 'INTERP_ERROR_S2_log_R', 'INTERP_ERROR_S2_spin_parameter', 'INTERP_ERROR_S2_lambda_CE_10cent', 'INTERP_ERROR_S1_MODEL01_mass', 'INTERP_ERROR_S1_MODEL01_spin', 'INTERP_ERROR_S1_MODEL01_m_disk_radiated', 'INTERP_ERROR_S1_MODEL05_mass', 'INTERP_ERROR_S1_MODEL05_spin', 'INTERP_ERROR_S1_MODEL05_m_disk_radiated', 'INTERP_ERROR_S1_MODEL06_mass', 'INTERP_ERROR_S1_MODEL06_spin', 'INTERP_ERROR_S1_MODEL06_m_disk_radiated', 'INTERP_ERROR_S1_MODEL10_mass', 'INTERP_ERROR_S1_MODEL10_spin', 'INTERP_ERROR_S1_MODEL10_m_disk_radiated'
    =====================  =====

.. table:: Check sets

    ======================  ======
    Set name                checks
    ======================  ======
    'CHECK_AFTER_CREATE'    
    'CHECK_AFTER_COMBINE'   'failure_rate'
    'CHECK_AFTER_EXTRA'     'CO_type', 'SN_type'
    'CHECK_AFTER_TRAINING'  
    ======================  ======
