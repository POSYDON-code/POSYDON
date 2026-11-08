.. _job_schedulers:

##############
Job schedulers
##############

POSYDON's command line tools write job scripts for a cluster's job scheduler:
``posydon-popsyn`` for population synthesis, ``posydon-setup-grid`` for MESA
grids and ``posydon-setup-pipeline`` for post-processing. Two schedulers are
supported:

- ``slurm``: SLURM (``sbatch``). This is the default.
- ``pbs``: PBS Pro (``qsub``), e.g. on ALCF Crux.

The tools describe their jobs independently of the scheduler, so switching
schedulers is a single setting:

.. list-table::
   :header-rows: 1
   :widths: 25 40 35

   * - Tool
     - Select the scheduler
     - Scheduler options
   * - ``posydon-popsyn setup``
     - ``--scheduler pbs``
     - ``--scheduler_option KEY=VALUE`` (repeatable)
   * - ``posydon-setup-grid``
     - ``--submission-type pbs``
     - keys in the ``[slurm]``/``[pbs]``/``[scheduler]`` ini section
   * - ``posydon-setup-pipeline``
     - ``SCHEDULER = 'pbs'`` in ``[account]``
     - upper case keys in ``[account]``, e.g. ``SELECT``

The general job settings (account, queue/partition, walltime, email) keep
their usual names for every scheduler. Generated files use the scheduler's
extension, e.g. ``step_1.slurm`` or ``step_1.pbs``, and the submission
scripts (``slurm_submit.sh``/``pbs_submit.sh``, ``run_grid.sh`` and
``run_pipeline.sh``) chain the jobs with the scheduler's dependencies.

``posydon-popsyn setup`` saves the job settings of a run in
``popsyn_jobs.json``. ``posydon-popsyn check`` uses this file to create rescue
jobs for the same scheduler. Arguments given to ``check`` override the saved
settings.


SLURM
=====

Each array task is a SLURM array task, started with ``srun``.

.. list-table:: SLURM options
   :header-rows: 1
   :widths: 25 75

   * - Option
     - Description
   * - ``exclude``
     - Comma separated list of nodes to exclude. For ``posydon-popsyn`` this
       is also available as ``--exclude``.
   * - ``launcher``
     - Command prefix starting a task (default: ``srun``).


PBS (e.g. ALCF Crux)
====================

Many PBS sites, including ALCF, schedule **whole nodes**: jobs never share a
node, and only a few jobs per project may be queued or running. POSYDON's jobs
are mostly arrays of small tasks (e.g. single-core population synthesis
batches). One PBS job per task would occupy a full node for every task, so
POSYDON *packs* array tasks onto nodes instead:

- One PBS job requests
  :math:`\lceil N_\mathrm{tasks}/\mathrm{tasks\_per\_node} \rceil` nodes, or
  at most ``max_nodes``.
- ``mpiexec`` starts the POSYDON task runner
  (``python -m posydon.CLI.schedulers.task_runner``) in each slot.
- Each runner runs the same command a SLURM array task would run, once for
  each array index assigned to it. It sets ``POSYDON_ARRAY_TASK_ID``,
  ``POSYDON_ARRAY_TASK_COUNT``, ``POSYDON_ARRAY_TASK_MIN`` and
  ``POSYDON_ARRAY_JOB_ID``, and writes one log file per task, as on SLURM.
- If there are more tasks than slots, each slot runs its tasks one after
  another. Choose the walltime accordingly.
- The job fails if any task failed. This means jobs depending on it with
  ``afterok`` (e.g. the merge job of a population) only start when all tasks
  succeeded.

The job scripts export the environment of ``qsub`` (``-V``), as SLURM does
by default, so submit them from a shell with your POSYDON environment
activated. The scripts change into the submission directory when they start.

.. list-table:: PBS options
   :header-rows: 1
   :widths: 25 75

   * - Option
     - Description
   * - ``select``
     - Additional resources of each chunk in ``-l select=N:...``, e.g.
       ``system=crux``.
   * - ``filesystems``
     - Value of ``-l filesystems=``, e.g. ``home:eagle`` (required at ALCF).
   * - ``tasks_per_node``
     - Array tasks per node. The default is the number that fits on a node,
       given the cores and memory per task (e.g. 28 popsyn tasks with
       ``--mem_per_cpu 9G``, 16 MESA runs with 4 cores and 4 GB per core).
   * - ``max_nodes``
     - Maximum number of nodes of an array job.
   * - ``cores_per_node``, ``memory_per_node``
     - Size of a node, used for the default ``tasks_per_node`` (default: 128
       and ``256G``, an ALCF Crux node).
   * - ``launcher``
     - Launcher template with the fields ``{ntasks}``, ``{ppn}``, ``{depth}``
       and ``{nodes}``. The default is ALCF's Cray PALS ``mpiexec``:
       ``mpiexec -n {ntasks} --ppn {ppn} --depth {depth} --cpu-bind depth``.
       For Open MPI, e.g., use ``mpirun -np {ntasks} --map-by ppr:{ppn}:node``.
   * - ``place``
     - Value of ``-l place=`` (default: ``scatter``).
   * - ``export_env``
     - Export the environment of ``qsub`` (default: ``True``).

ALCF Crux examples
------------------

Population synthesis:

.. code-block:: bash

    posydon-popsyn setup population_params.ini --job_array 200 \
        --scheduler pbs --queue workq-route --account <PROJECT> \
        --scheduler_option select=system=crux \
        --scheduler_option filesystems=home:eagle \
        --walltime 12:00:00 --merge_walltime 01:00:00 --mem_per_cpu 9G
    ./pbs_submit.sh

MESA grid (``posydon-setup-grid --submission-type pbs --inifile grid.ini``),
in the grid ini file:

.. code-block:: ini

    [pbs]
    number_of_nodes=1
    number_of_mpi_tasks=1
    number_of_cpus_per_task=4
    job_array=True
    partition=workq-route
    account={PROJECT}
    walltime='24:00:00'
    email={YOUR_EMAIL_ADDRESS}
    ; PBS options
    select=system=crux
    filesystems=home:eagle
    max_nodes=8

Post-processing pipeline, in the ``[account]`` section:

.. code-block:: ini

    [account]
        ACCOUNT = '<PROJECT>'
        PARTITION = 'workq-route'
        WALLTIME = '12:00:00'
        MAILTYPE = 'FAIL'
        EMAIL = 'me@university.edu'
        SCHEDULER = 'pbs'
        SELECT = 'system=crux'
        FILESYSTEMS = 'home:eagle'

.. note::

    On ALCF Crux, ``workq`` jobs are limited to 24 hours, so MESA runs that
    need longer will not finish within a single job.


Adding a scheduler
==================

Jobs are described by :class:`posydon.CLI.schedulers.base.JobSpec` and
:class:`~posydon.CLI.schedulers.base.ArraySpec`. Commands and log file names
can use the placeholders ``{task_id}``, ``{job_id}`` and ``{job_end}``.

To support another scheduler:

1. Subclass :class:`posydon.CLI.schedulers.base.Scheduler`.
2. Set ``name``, ``extension``, ``submit_program`` and ``options``.
3. Implement ``directives``, ``body``, ``submit_line`` and ``mpi_launcher``.
4. Add the class to ``SCHEDULERS`` in ``posydon/CLI/schedulers/__init__.py``.

All tools then accept the new scheduler name.
