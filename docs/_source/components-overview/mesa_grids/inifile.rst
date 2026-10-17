.. _inifile:

#########################################
ini file documentation for POSYDON GRIDS
#########################################

How to write a configuration file
==================================

Grids are set up with the ``posydon-grid`` command (``posydon-grid types``,
``init``, ``explain``, ``setup``; ``posydon-grid run`` runs a grid). The
older ``posydon-setup-grid`` still works but is deprecated: it prints the
``posydon-grid`` command to use instead.

``posydon-grid setup`` cannot run without a
configuration file. This :samp:`ini` file needs to specify all aspects of the MESA
simulation that are static, i.e. those parts that are shared amongst every MESA
simulation in your grid. The design of the :samp:`ini` file file is such that it will
start with the MESA defaults for every possible MESA parameter, stack on top of
that the POSYDON suggested defaults, and then finally, if supplied, stack user
specific parameters on top of both the MESA defaults and the POSYDON defaults.
This scheme requires a high-performance computing environment with SLURM
installed.

Below we describe each section of the :samp:`ini` file and what each value corresponds
to. In general, lines with capitalized quantities surrounded by {} are
quantities that need to be adjusted for your specific case.

The default :samp:`ini` file is in POSYDON-MESA-INLISTS:
`grid_params.ini <https://github.com/POSYDON-code/POSYDON-MESA-INLISTS/blob/v2_sources/templates/grid_params.ini>`_

POSYDON no longer ships the ``grid_params/`` folder. Ini files still name
its files as ``${posydon_github_root}/grid_params/...`` (MESA defaults, ZAMS
models, column lists); when such a path does not exist, ``posydon-grid``
reads the same file from POSYDON-MESA-INLISTS at the commit POSYDON pins
(``posydon/grids/setup/support.py``) and says so in a warning.

The sections below the next one describe this (legacy) format. A shorter,
new format is also accepted; it is described first.

.. _inifile_new:

The new format
==============

An :samp:`ini` file with a ``[grid]`` section is read in the new format.
Write a starter file for a grid type with

.. code-block:: bash

    posydon-grid types
    posydon-grid init CO-HeMS grid.ini

and fill in the ``<...>`` values. A complete example:

.. code-block:: ini

    [grid]
    type = CO-HeMS              ; see --list-grid-types
    mode = fixed                ; fixed | dynamic (needs psycris_inifile)
    points = grid.csv           ; one row per MESA run
    run_directory = .

    [scenario]                  ; which POSYDON inlists to use
    release = v2.0
    rerun = matthias_lbv_thermo_dedt_noAgeLimit   ; optional
    commit = 8482d7d            ; or ref = <branch or tag>; one is required

    [user_inlists]              ; your own MESA inlists, on top
    star1_controls = my_inlist1

    [output]
    history_interval = 1
    star1 = history, final_model
    star2 = none
    binary = history

    [submission]
    backend = slurm_array       ; shell | slurm_array | slurm_mpi
    partition = {PARTITION}
    account = {ACCOUNT}
    walltime = 2-00:00:00
    email = {EMAIL}
    cpus_per_task = 4

Relative paths are relative to the :samp:`ini` file.

**The scenario** selects the POSYDON inlists. ``release`` is a folder of
`POSYDON-MESA-INLISTS <https://github.com/POSYDON-code/POSYDON-MESA-INLISTS>`_
(``r11701/releases/<release>/``) and ``rerun`` one of its physics variants
(``reruns/<rerun>/``, which holds only the keys and files it changes). There
is no default version of the inlist repository: give ``ref`` (a branch or
tag) or ``commit``. ``source = path`` with ``inlist_repo = <folder>`` uses a
local copy as it is, without git. ``mesa`` names the MESA version if a
release exists for several.

**Layers.** For every inlist the values are stacked, a later layer winning
key by key:

1. MESA's defaults of the release's MESA installation
   (``r11701/mesa/<installation>/``), for every inlist ``[mesa_inlists]``
   does not list files for;
2. the release's files for the grid type (and the rerun's on top of them);
3. your ``[user_inlists]``;
4. the keys the setup writes itself (wiring, output options, ZAMS models).

``[user_inlists]`` takes ``binary``, ``star1``, ``star2`` (both namelists),
a single namelist such as ``star1_controls`` or ``binary_job``,
``star1_formation_step0`` (on top of that formation step), and whole files
for ``history_columns``, ``binary_history_columns``, ``profile_columns``,
``run_binary_extras`` and ``run_star_extras``. Several files for one key are
separated by commas and applied in order. In a single-star grid, which
runs star 1 as its formation steps, ``star1`` files apply to every step,
and ``[output] star1`` applies to its last step, which runs the star to
its end (``final_model`` saves ``final_star1.mod``).

A ZAMS model set in ``[grid]`` (``zams``, ``zams_star1``, ``zams_star2``)
is used instead of the release's. ``[submission]`` settings left out
(``account``, ``partition``, ``walltime``, ``email``) are left out of the
SLURM scripts, so the cluster's defaults apply.

A whole file from a higher layer replaces the lower one: your
``run_star_extras`` replaces POSYDON's wherever a star runs (in the binary
and in the formation steps). Fortran sources are copied under the name MESA
compiles (``run_star_extras.f``, ``run_binary_extras.f``, ...), so your file
can have any name. MESA's makefiles and stock sources are filled in for
every ``[extras]`` key the file does not set.

To see where every value comes from without writing anything, run

.. code-block:: bash

    posydon-grid explain grid.ini

It lists, for each inlist, the keys POSYDON or you set, with the layer and
file of each value and the values it overrides.

**Command line.** A new-format file is the only input:

.. code-block:: bash

    posydon-grid setup grid.ini

It says the mode, run directory and backend itself (``[grid] mode``,
``run_directory``, ``[submission] backend``), so ``--grid-type``,
``--run-directory`` and ``--submission-type`` are refused for it. A legacy
file cannot say these, so for it they are given on the command line:

.. code-block:: bash

    posydon-grid setup grid_params.ini --grid-type fixed --submission-type slurm

**The run directory** gets a ``setup_record.ini`` with the :samp:`ini` file,
the POSYDON version, ``$MESA_DIR``, the scenario and the exact commit of the
inlist repository it used.

.. _inifile_slurm:


[slurm]
-------

This section designates all the SLURM-related parameters.

.. code-block:: ini

    [slurm]

    ; Number of nodes you would like to request
    number_of_nodes=1

    ; This is the number of processes that will be running MESA simulations
    number_of_mpi_tasks=1

    ; This is the number of cores you will let each MESA simulation use
    number_of_cpus_per_task=4

    ; run as job array
    ; cannot be false if running with --run-type sample
    job_array=True

    ; user name
    user={USERNAME}

    ; partition
    partition={PARTITION}

    ; account
    ; If running on geneva cluster set this to "default"
    account={ACCOUNT}

    ; wall-time
    walltime='2-00:00:00'

    ; work-directory: place, where MESA writes during runtime
    work_dir={WORK_DIRECTORY}

    ;email
    email={YOUR_EMAIL_ADDRESS}

    ; new group to be set with write permission
    ; if empty string no changes on group and permission
    ; group on yggdrasil: GL_S_Astro_POSYDON
    ; group on Quest: b1119
    newgroup=''

.. table:: general settings

    =======================  ===========
    Setting name             Description
    =======================  ===========
    number_of_nodes          The number of nodes each job will request
    number_of_mpi_tasks      (outdated) The number of parallel processes ganerates by MPI (recommended to not change)
    number_of_cpus_per_task  The number of CPUs each job will request (in best this should align with the number of threads set for MESA)
    job_array                Set to `True` to run a SLURM job array
    user                     User name
    partition                SLURM partition to run the jobs on
    account                  SLURM account to run the jobs on
    walltime                 Maximum time for each run, otherwise it will get cancelled by SLURM to avoid never ending MESA runs going on too short time steps
    work_dir                 The path to the place, where the data will be writting during runtime (it is recommended to use fast local node storage here), afterwards the data will be copied/moved to the current directory (uncommend this line to write always to the final location)
    email                    Your email address to receive notifications from SLURM
    newgroup                 The name of the owning group all the files should get (leave empty if no changes are needed here)
    =======================  ===========


[mesa_inlists]
--------------

This section designates all the basic MESA-specific parameters.

.. code-block:: ini

    [mesa_inlists]
    ; inlist types: binary_control,binary_job,star_job,star_control

    ; Please leave MESA defaults lines alone, please fill in the path to
    ; you local clone of POSYDON
    posydon_github_root={PATH_TO_POSYDON_DIRECTORY}

    ; There are a number of ways to build the physics of your MESA grid
    ; the first way is to point the sections below to your own MESA inlists,
    ; and/or the POSYDON default inlists (versions of which can be found in the
    ; following repo: https://github.com/POSYDON-code/POSYDON-MESA-INLISTS)
    ; you can also supply a scenario using syntax such as below, and the setup script
    ; will automatically find the inlists from POSYDON
    ; you want to use based on the git tag/commit and the scenario
    ; (in this case you are simulating MS-MS binaries)
    ; NOTE: You can use the scenario logic below and *still* supply your own local
    ; user mesa inlists that will overwrite or tweak some of the physics associated
    ; with the scenario.

    scenario = ['posydon', 'master-9ddb61bb0c482399fa5a41dd22fde41ccd8175d9', 'CO-H_star']

    ; zams_filename if a zams_filename is supplied this supercedes any star1 or star2 formation inlists
    ; and skips to running the binary with this pre-computed zams model.
    zams_filename = ${posydon_github_root}/grid_params/POSYDON-MESA-INLISTS/r11701/ZAMS_models/zams_z0.0142m2_y0.2703.data

    ; single_star_grid, this boolean, when True, will take the inlist1 from the binary mesa inlist section
    ; and run in a single star grid configuration
    single_star_grid = False

    ; Star1 formation - star1_job
    ; star1_formation_job_mesa_defaults = ${posydon_github_root}/grid_params/defaults/r11701/star/star_job.defaults
    ; star1_formation_job_posydon_defaults = ${user_template_root}/inlist1
    ; star1_formation_job_user = None

    ; Star2 formation - star2_job
    ; star2_formation_job_mesa_defaults = ${posydon_github_root}/grid_params/defaults/r11701/star/star_job.defaults
    ; star2_formation_job_posydon_defaults = ${user_template_root}/inlist2
    ; star2_formation_job_user = None

    ; Star1 formation - star1_control
    ; star1_formation_controls_mesa_defaults = ${posydon_github_root}/grid_params/defaults/r11701/star/controls.defaults
    ; star1_formation_controls_posydon_defaults = ${user_template_root}/inlist1
    ; star1_formation_controls_user = None

    ; Star2 formation - star2_control
    ; star2_formation_controls_mesa_defaults = ${posydon_github_root}/grid_params/defaults/r11701/star/controls.defaults
    ; star2_formation_controls_posydon_defaults = ${user_template_root}/inlist2
    ; star2_formation_controls_user = None

    ; binary_control
    binary_controls_mesa_defaults = ${posydon_github_root}/grid_params/defaults/r11701/binary/binary_controls.defaults
    ; binary_controls_posydon_defaults = ${posydon_github_root}/grid_params/POSYDON-MESA-INLISTS/r11701/default_common_inlists/binary/inlist_project
    ; binary_controls_user = ${user_template_root}/binary/inlist_project

    ; binary_job
    binary_job_mesa_defaults = ${posydon_github_root}/grid_params/defaults/r11701/binary/binary_job.defaults
    ; binary_job_posydon_defaults = ${posydon_github_root}/grid_params/POSYDON-MESA-INLISTS/r11701/default_common_inlists/binary/inlist_project
    ; binary_job_user = ${user_template_root}/binary/inlist_project

    ; star1_job
    star1_job_mesa_defaults = ${posydon_github_root}/grid_params/defaults/r11701/star/star_job.defaults
    ; star1_job_posydon_defaults = ${posydon_github_root}/grid_params/POSYDON-MESA-INLISTS/r11701/default_common_inlists/binary/inlist1
    ; star1_job_user =  ${user_template_root}/binary/inlist1

    ; star1_control
    star1_controls_mesa_defaults = ${posydon_github_root}/grid_params/defaults/r11701/star/controls.defaults
    ; star1_controls_posydon_defaults = ${posydon_github_root}/grid_params/POSYDON-MESA-INLISTS/r11701/default_common_inlists/binary/inlist1
    ; star1_controls_user = ${user_template_root}/binary/inlist1

    ; star2_job
    star2_job_mesa_defaults = ${posydon_github_root}/grid_params/defaults/r11701/star/star_job.defaults
    ; star2_job_posydon_defaults = ${user_template_root}/binary/inlist2
    ; star2_job_user = None

    ; star2_control
    star2_controls_mesa_defaults = ${posydon_github_root}/grid_params/defaults/r11701/star/controls.defaults
    ; star2_controls_posydon_defaults = ${user_template_root}/binary/inlist2
    ; star2_controls_user = None

    ; star history columns
    star_history_columns = ${posydon_github_root}/grid_params/POSYDON-MESA-INLISTS/r11701/default_common_inlists/history_columns.list

    ; binary history columns
    binary_history_columns = ${posydon_github_root}/grid_params/POSYDON-MESA-INLISTS/r11701/default_common_inlists/binary_history_columns.list

    ; profile columns
    profile_columns = ${posydon_github_root}/grid_params/POSYDON-MESA-INLISTS/r11701/default_common_inlists/profile_columns.list

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;;;;; MESA OUTPUT CONTROLS ;;;;;;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    ; controls how often MESA prints out the history of the evolution
    history_interval = 1

    ; Save binary history (history file will be named: )
    binary_history = True

    ; save history of star1
    history_star1 = True
    ;save final profile of star1
    final_profile_star1 = False
    ; save final model of star1
    final_model_star1 = True

    ; save history of star2
    history_star2 = False
    ; save profile of star2
    final_profile_star2 = False
    ; save final model of star2
    final_model_star2 = False

.. table:: settings for the MESA inlist

    =========================================  ===========
    Setting name                               Description
    =========================================  ===========
    posydon_github_root                        The path to your used POSYDON version
    scenario                                   List containing multiple information: 1) the source ('posydon' or 'user'(for future use)), 2) the git commit (the branch and full git hash for the inlist submodule separated by a dash), 3) the systems type ('HMS-HMS', 'CO-H_star', 'CO-He_star')
    zams_filename                              The location of the file containing the ZAMS models
    single_star_grid                           Flag to indicate single star or binary evolution
    star1_formation_job_mesa_defaults          (outdated) Path to the MESA job section defaults to form star 1
    star1_formation_job_posydon_defaults       (outdated) Path to the MESA job section inlist of POSYDON to form star 1
    star1_formation_job_user                   (outdated) Path to the MESA job section inlist of the user to form star 1
    star2_formation_job_mesa_defaults          (outdated) Path to the MESA job section defaults to form star 2
    star2_formation_job_posydon_defaults       (outdated) Path to the MESA job section inlist of POSYDON to form star 2
    star2_formation_job_user                   (outdated) Path to the MESA job section inlist of the user to form star 2
    star1_formation_controls_mesa_defaults     (outdated) Path to the MESA controls section defaults to form star 1
    star1_formation_controls_posydon_defaults  (outdated) Path to the MESA controls section inlist of POSYDON to form star 1
    star1_formation_controls_user              (outdated) Path to the MESA controls section inlist of the user to form star 1
    star2_formation_controls_mesa_defaults     (outdated) Path to the MESA controls section defaults to form star 2
    star2_formation_controls_posydon_defaults  (outdated) Path to the MESA controls section inlist of POSYDON to form star 2
    star2_formation_controls_user              (outdated) Path to the MESA controls section inlist of the user to form star 2
    binary_controls_mesa_defaults              (outdated) Path to the MESA controls section defaults to evolve the binary
    binary_controls_posydon_defaults           (outdated) Path to the MESA controls section inlist of POSYDON to evolve the binary
    binary_controls_user                       (outdated) Path to the MESA controls section inlist of the user to evolve the binary
    binary_job_mesa_defaults                   (outdated) Path to the MESA job section defaults to evolve the binary
    binary_job_posydon_defaults                (outdated) Path to the MESA job section inlist of POSYDON to evolve the binary
    binary_job_user                            (outdated) Path to the MESA job section inlist of the user to evolve the binary
    star1_job_mesa_defaults                    (outdated) Path to the MESA job section defaults to evolve star 1
    star1_job_posydon_defaults                 (outdated) Path to the MESA job section inlist of POSYDON to evolve star 1
    star1_job_user                             (outdated) Path to the MESA job section inlist of the user to evolve star 1
    star1_controls_mesa_defaults               (outdated) Path to the MESA controls section defaults to evolve star 1
    star1_controls_posydon_defaults            (outdated) Path to the MESA controls section inlist of POSYDON to evolve star 1
    star1_controls_user                        (outdated) Path to the MESA controls section inlist of the user to evolve star 1
    star2_job_mesa_defaults                    (outdated) Path to the MESA job section defaults to evolve star 2
    star2_job_posydon_defaults                 (outdated) Path to the MESA job section inlist of POSYDON to evolve star 2
    star2_job_user                             (outdated) Path to the MESA job section inlist of the user to evolve star 2
    star2_controls_mesa_defaults               (outdated) Path to the MESA controls section defaults to evolve star 2
    star2_controls_posydon_defaults            (outdated) Path to the MESA controls section inlist of POSYDON to evolve star 2
    star2_controls_user                        (outdated) Path to the MESA controls section inlist of the user to evolve star 2
    star_history_columns                       (outdated) Path to the history columns list of the stars
    binary_history_columns                     (outdated) Path to the history columns list of the binary
    profile_columns                            (outdated) Path to the profile columns list to write the final stellar profile
    history_interval                           Interval how often MESA will add a model to the star's histories
    binary_history                             Interval how often MESA will add a model to the binary history
    history_star1                              Flag, whether the history of star 1 should be saved
    final_profile_star1                        (outdated, done in :samp:`run_star_extras.f`) Flag, whether the final profil of star 1 should be saved
    final_model_star1                          Flag, whether the final model of star 1 should be saved
    history_star2                              Flag, whether the history of star 2 should be saved
    final_profile_star2                        (outdated, done in :samp:`run_star_extras.f`) Flag, whether the final profil of star 2 should be saved
    final_model_star2                          Flag, whether the final model of star 2 should be saved
    =========================================  ===========


[mesa_extras]
-------------

This section designates all the parameters for MESA makefiles and fortran files.

.. code-block:: ini

    [mesa_extras]
    ; path to MESA makefile for executable binary and star
    makefile_binary = ${MESA_DIR}/binary/work/make/makefile
    makefile_star = ${MESA_DIR}/star/work/make/makefile

    ; N.B. Normally system_type will determine which extras file gets used.
    ; posydon has a set of approved extras files for given types of systems
    ; and it will use these extra files by default but you may supply your own
    ; if you wish.

    ; user specified binary extra
    mesa_binary_extras = ${MESA_DIR}/binary/work/src/run_binary_extras.f
    ; user_binary_extras = ${mesa_inlists:posydon_github_root}/grid_params/POSYDON-MESA-INLISTS/r11701/default_common_inlists/binary/src/run_binary_extras.f

    ; user specified star extra - these go into the binary/src/ directory
    mesa_star_binary_extras = ${MESA_DIR}/binary/work/src/run_star_extras.f
    ; user_star_binary_extras =${mesa_inlists:posydon_github_root}/grid_params/POSYDON-MESA-INLISTS/r11701/default_common_inlists/binary/src/run_star_extras.f

    ; user specified star extras - these are for single star formation (e.g., pre-MS evolution)
    mesa_star1_extras = ${MESA_DIR}/star/work/src/run_star_extras.f
    ; user_star1_extras = ${mesa_inlists:posydon_github_root}/grid_params/POSYDON-MESA-INLISTS/r11701/default_common_inlists/binary/src/run_star_extras.f

    mesa_star2_extras = ${MESA_DIR}/star/work/src/run_star_extras.f
    ; user_star2_extras = ${mesa_inlists:posydon_github_root}/grid_params/POSYDON-MESA-INLISTS/r11701/default_common_inlists/binary/src/run_star_extras.f

    ; binary_run.f
    binary_run = ${MESA_DIR}/binary/work/src/binary_run.f

    ; star_run.f
    star_run = ${MESA_DIR}/star/work/src/run.f

.. table:: settings for the MESA extras

    =======================  ===========
    Setting name             Description
    =======================  ===========
    makefile_binary          Path to the make file of MESA's binary module
    makefile_star            Path to the make file of MESA's star module
    mesa_binary_extras       Path to MESA's binary module default :samp:`run_binary_extras.f`
    user_binary_extras       Path to the users/POSYDON :samp:`run_binary_extras.f`
    mesa_star_binary_extras  Path to MESA's binary module default :samp:`run_star_extras.f`
    user_star_binary_extras  Path to the users/POSYDON :samp:`run_star_extras.f`
    mesa_star1_extras        Path to MESA's star module default :samp:`run_star_extras.f`
    user_star1_extras        Path to the users/POSYDON :samp:`run_star_extras.f`
    mesa_star2_extras        Path to MESA's star module default :samp:`run_star_extras.f`
    user_star2_extras        Path to the users/POSYDON :samp:`run_star_extras.f`
    binary_run               Path to MESA's binary module :samp:`binary_run.f`
    star_run                 Path to MESA's star module :samp:`run.f`
    =======================  ===========


[run_parameters]
----------------

This section designates the run parameters for a grid.

.. code-block:: ini

    [run_parameters]
    ; If running posydon-run-grid with option --grid-type fixed
    ; then the grid is a file with all the different samples you would like to
    ; run MESA on.
    ; If posydon-make-grid is run with --grid-type dynamic, then grid is
    ; a file of pre-run MESA simulations from which you will generate new samples to
    ; run MESA on (i.e. generate grid points on the fly).

    grid = {PATH_TO_GRID}

The :samp:`grid` specifies where to find the csv file to read the runs from.
