.. _posydon_grid:

########################
The posydon-grid command
########################

``posydon-grid`` sets up grids of MESA runs from an :samp:`ini` file. A grid
is set up in four steps:

.. code-block:: bash

    posydon-grid init HMS-HMS grid.ini     # 1. write a starter ini file
    # 2. edit grid.ini, write grid.csv
    posydon-grid explain grid.ini          # 3. check where each inlist value comes from
    posydon-grid setup grid.ini            # 4. write the run directory

and then run, with the script ``setup`` wrote (``run_grid.sh`` on a SLURM
cluster, ``grid_command.sh`` without SLURM). The
:ref:`how-to guides <MESA-grids>` go through this for several kinds of
grid; :ref:`inifile` describes the :samp:`ini` file.


What you need
=============

- **MESA.** ``$MESA_DIR`` and ``$MESASDK_ROOT`` must be set, and the MESA SDK
  loaded, before ``posydon-grid setup``: it compiles MESA into the run
  directory, and writes both paths into the run scripts. The POSYDON v2.0
  grids use a modified MESA r11701.
- **The POSYDON inlists**: a clone of `POSYDON-MESA-INLISTS
  <https://github.com/POSYDON-code/POSYDON-MESA-INLISTS>`_ on its
  ``v2_sources`` branch, named by ``[scenario] inlist_repo`` (default:
  ``~/.posydon_mesa_inlists``):

  .. code-block:: bash

      git clone --branch v2_sources https://github.com/POSYDON-code/POSYDON-MESA-INLISTS.git ~/POSYDON-MESA-INLISTS

  ``posydon-grid`` reads the clone as it is and never runs git on it:
  updating it (``git pull``), or pinning it to a commit, is up to you.

``posydon-grid types``, ``init`` and ``explain`` need neither MESA nor
SLURM, so you can prepare and check an :samp:`ini` file on a laptop.


Commands
========

``posydon-grid types``
    Lists the grid types: ``HMS-HMS``, ``CO-HMS_RLO``, ``CO-HeMS``,
    ``HeMS-HMS``, ``single_HMS`` and ``single_HeMS``.

``posydon-grid init GRID_TYPE INIFILE``
    Writes a starter :samp:`ini` file with every option of the grid type,
    each with a short note. Options left at their default are commented out
    (``; key = value``); remove the ``;`` to set one, and fill in the
    ``<...>`` values. It does not overwrite an existing file.

``posydon-grid explain INIFILE``
    Shows how the inlists will be built, and writes nothing. For every
    namelist (``binary_controls``, ``star1_controls``, ...,
    ``star1_formation_step0_controls``, ...) it lists each key that POSYDON
    or you set, with the layer and file of its value, and below it the
    values it overrides::

        max_age = -1  [posydon_inlists] r11701/releases/v2.0/reruns/matthias_no_age_limit/common/binary/inlist1
            over 1.4e10  [posydon_inlists] r11701/releases/v2.0/common/binary/inlist1
            over 1d36  [mesa_inlists] r11701/mesa/default/defaults/star/controls.defaults

    The output is long; save it and search it
    (``posydon-grid explain grid.ini > explain.txt``). Keys left at MESA's
    default are not listed. The keys the setup writes itself (the ZAMS
    model, the output options, the wiring between the inlists) are added on
    top when the inlists are written.

``posydon-grid setup INIFILE``
    Writes the run directory (``[grid] run_directory``, relative to the
    :samp:`ini` file): it compiles MESA, writes the inlists, copies the
    column lists and Fortran sources, writes the run scripts and the record.
    Running it again on the same run directory replaces the files it writes,
    with a warning for each, and leaves everything else (such as the output
    of earlier runs) in place.

``posydon-grid run ...``
    Runs ``posydon-run-grid``, which runs the grid points, with the same
    arguments. The run scripts call ``posydon-run-grid`` themselves, so you
    rarely need it.

A new-format :samp:`ini` file says everything ``setup`` needs. The options
``--grid-type``, ``--run-directory`` and ``--submission-type`` are for
legacy :samp:`ini` files only (below).

A mistake in the :samp:`ini` file or the environment (a wrong value, a
missing file, ``$MESASDK_ROOT`` not set, ``init`` onto an existing file) is
printed as one line, such as ``posydon-grid: error: [submission] backend:
...``, and the command exits with status 2. A Python traceback means a bug
in ``posydon-grid``; please report it.


The layers
==========

Every inlist is stacked from layers; a later layer wins key by key:

1. ``[mesa_inlists]``: MESA's defaults, of the MESA installation the release
   names (``r11701/mesa/<installation>/defaults/`` in the inlist
   repository), unless your ``[mesa_inlists]`` section gives other files;
2. ``[posydon_inlists]``: the release's inlists for the grid type
   (``r11701/releases/<release>/``), and the rerun's changes on top
   (``reruns/<rerun>/``), as ``[scenario]`` selects them;
3. ``[user_inlists]``: your own inlists, from the ``[user_inlists]``
   section;
4. the keys the setup writes itself.

Whole files (column lists, Fortran sources) are not stacked: the highest
layer that has one wins. ``posydon-grid explain`` and the record use the
layer names above.

.. note::

    MESA stops with a namelist error, when the run starts, at a key it does
    not know. So ``posydon-grid explain`` and ``setup`` warn about every
    key of your files (``[user_inlists]`` and ``[mesa_inlists]``) that is
    not a key of the namelist it goes into: a misspelled key
    (``max_modle_number`` in ``&controls``), or a key in the wrong namelist
    (``max_model_number`` in ``&star_job``). The warning names the key, the
    namelist and the file. It does not stop the setup.

    The keys of a namelist are those of MESA's ``.defaults`` file for it, of
    the MESA version the grid uses: the installation the release names
    (``r11701/mesa/<installation>/defaults/``), so a release for another
    MESA version is checked against that version. This holds also where
    your ``[mesa_inlists]`` replaces MESA's defaults: your file is checked
    against the installation's. Without a ``[scenario] release`` (or for a
    legacy :samp:`ini` file), the reference is the ``.defaults`` files of
    your ``[mesa_inlists]`` (legacy: ``*_mesa_defaults``), and a namelist
    without one is not checked. Names are compared without case and
    without an array index: ``x_logical_ctrl(5)`` is ``x_logical_ctrl``.
    The release's own files are not checked by the setup.


The run directory
=================

``posydon-grid setup`` writes:

=============================== ============================================
``binary/``, ``star1/``,        MESA's work folders: the compiled executables,
``star2/``                      the inlists (``binary/inlist_project``,
                                ``binary/inlist1``, ``binary/inlist2``;
                                ``star1/inlist_step0``, ... for formation
                                steps and single stars) and ``src/``
``column_lists/``               the history and profile columns
``mk``                          the script that compiled MESA
``setup_record.ini``            what the run directory was set up from
=============================== ============================================

and, depending on ``[submission] backend``:

================ ===============================================================
``slurm_array``  ``job_array_grid_submit.slurm`` (one array task per grid
                 point), ``cleanup.slurm`` (compresses the output afterwards)
                 and ``run_grid.sh``, which submits both. Run
                 ``cd run && ./run_grid.sh``.
``slurm_mpi``    ``mpi_grid_submit.slurm``, ``cleanup.slurm`` and
                 ``run_grid.sh``: one MPI job, for dynamic grids.
``shell``        ``grid_command.sh``, which runs the points one after the other
                 on the machine it is started on (``job_array = true``), or as
                 one ``mpirun`` for a dynamic grid (``job_array = false``).
                 Run ``cd run && ./grid_command.sh``.
================ ===============================================================

Each grid point runs in its own folder in the run directory, named after
its parameters and ending in ``_grid_index_<row>``; MESA's screen output
is in ``out.txt`` there.


The record
==========

``setup_record.ini`` says what the run directory was set up from, so a grid
can be traced, compared and set up again:

- ``[setup]``: the :samp:`ini` file, its format, the grid type and mode, when,
  on which host and with which command;
- ``[posydon]``: the POSYDON that ran: version, location, commit, whether its
  code had local changes (``dirty``), and Python;
- ``[mesa]``: ``$MESA_DIR`` and ``$MESASDK_ROOT``;
- ``[scenario]``: the release and rerun, and the inlist repository as it was
  read: folder, remote, commit, whether the files used differed from that
  commit, and a checksum of those files (``inlist_fingerprint``). Two grids
  with the same fingerprint read identical files;
- ``[inlists]``: for every written namelist, the files stacked into it,
  lowest layer first;
- ``[whole_files]``: for every column list and Fortran source, the file it
  was copied from.

To set up the same grid again from the same files, keep the :samp:`ini` file
and check out the inlist repository at ``inlist_commit``.


Legacy ini files
================

:samp:`ini` files of the earlier format (with ``[slurm]``,
``[mesa_inlists]``, ``[mesa_extras]`` and ``[run_parameters]`` sections) are
still read. They cannot say the grid mode, run directory or backend, so
these are given on the command line:

.. code-block:: bash

    posydon-grid setup grid_params.ini --grid-type fixed --submission-type slurm

The legacy format checks out the inlist repository at the commit its
``scenario`` names. ``posydon-setup-grid``, the earlier command, still works
for them but is deprecated: it prints the ``posydon-grid`` command to use
instead. For new grids, use the new format: ``posydon-grid init``.
