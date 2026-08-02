.. _config-reference:

Config Reference
================

POSYDON is driven by ``ini`` files that select which physical assumptions to
use, how many binaries to evolve, how much memory/disk to use (via ``dump_rate``),
where to write the output, and so on. There are three families of ``ini`` files,
each documented on its own page. This page is the central index that points you
to the correct one.

.. list-table::
   :header-rows: 1

   * - Family
     - What it configures
     - Page
   * - Population
     - The population-synthesis run: number of binaries, seeds, metallicity,
       ``dump_rate``, output directories. Loaded by a ``PopulationRunner``.
     - :ref:`population-parameters <pop-params-guide>`
   * - Post-processing
     - The pipeline that turns raw MESA output into a ``PSyGrid``: which steps
       to run, their order, inputs and outputs.
     - :ref:`pipeline-ini <pipeline_ini>`
   * - MESA
     - The parameters for the MESA simulations and the MESA grid completion
       API: inlists, directories, scheduling.
     - :ref:`mesa-inifile <inifile>`

How the three families relate
-----------------------------

The three ``ini`` families correspond to three stages of the workflow
(see :ref:`workflow`):

1. run MESA simulations on the grids: :ref:`MESA inifile <inifile>`;
2. turn the grids into a ``PSyGrid``: the :ref:`pipeline inifile <pipeline_ini>`;
3. run the final population synthesis: the
   :ref:`population parameters <pop-params-guide>` file.

A quick example
---------------------------------

The default population configuration ships with the code. Make a copy, edit it,
and hand it to a ``PopulationRunner``:

.. code-block:: python

   from posydon.popsyn.synthetic_population import PopulationRunner
   poprun = PopulationRunner('my_population.ini', verbose=True)
   poprun.evolve()

See each page linked above for the full list of keys and their defaults.

Command-line executables
--------------------------

The bin executables (``posydon_setup_grid``, ``posydon_run_grid``,
``posydon_setup_pipeline``, ``posydon_run_pipeline``, ``posydon_popsyn``) read
these same ``ini`` files from the command line. See the
:doc:`bin executables <../api_reference/bin>` reference.

.. seealso::

   * :ref:`pop-params-guide` - the population run configuration.
   * :ref:`pipeline_ini` - the post-processing configuration.
   * :ref:`inifile` - the MESA inifile command line.
