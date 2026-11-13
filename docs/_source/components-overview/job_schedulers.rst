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
``run_pipeline.sh``) chain the jobs with the scheduler's dependencies. They
stop with an error message as soon as a submission fails, so no job is
submitted without the jobs it depends on.

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
- If there are more tasks than slots, a slot that becomes free takes the
  next task that no other slot has taken yet (a work queue). All slots stay
  busy until the tasks run out, even if the tasks' run times differ a lot,
  as for MESA runs. Choose the walltime accordingly (see
  :ref:`large_grids`).
- The slots are spread evenly over the cores of a node, e.g. 56
  single-core tasks get 2 cores each, so both sockets of a node (and their
  memory) are used.
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
       ``{depth}`` is the number of cores per slot: the cores of a node
       divided evenly between the slots, but at least the cores per task.
       For Open MPI, e.g., use ``mpirun -np {ntasks} --map-by ppr:{ppn}:node``.
   * - ``place``
     - Value of ``-l place=`` (default: ``scatter``).
   * - ``export_env``
     - Export the environment of ``qsub`` (default: ``True``).
   * - ``rerunnable``
     - ``True`` writes ``#PBS -r y``, ``False`` writes ``#PBS -r n``; not set
       by default. A rerunnable job is restarted **from the beginning**, e.g.
       after being preempted (see below).
   * - ``task_stats``
     - Directory in which each array task's start, end, CPU time and peak
       memory are recorded (see :ref:`task_stats`). Off by default.
   * - ``task_order``
     - ``dynamic`` (default): with more tasks than slots, a free slot takes
       the next task nobody has taken yet. ``static``: slot *r* runs the tasks
       *r*, *r* + slots, ... assigned when the job starts.
   * - ``queue_dir``
     - Directory of the work queues (default: ``.task_queue`` in the
       submission directory). It must be on a file system shared by all
       nodes, e.g. Eagle. Each run of a job uses its own subdirectory with a
       single small file; the directory can be deleted when no job is
       running.

Make sure ``mpiexec`` is the right one
--------------------------------------

The default ``launcher`` uses the first ``mpiexec`` in your ``PATH``. Installing
``mpi4py`` with conda usually installs conda's own MPI (e.g. Open MPI) into
the environment, whose ``mpiexec`` then comes before the system's. On ALCF
systems the job fails with ``An unrecognized option was included on the
mpiexec command line: --ppn``. Check with ``which -a mpiexec`` and, if needed,
give the system's launcher with its full path, e.g. on Crux:

.. code-block:: bash

    --scheduler_option "launcher=/opt/cray/pals/1.8/bin/mpiexec -n {ntasks} --ppn {ppn} --depth {depth} --cpu-bind depth"

(or ``launcher=...`` in the grid ini, ``LAUNCHER = '...'`` in the pipeline's
``[account]`` section). The path may change when the system is updated.


ALCF Crux queues
----------------

.. list-table::
   :header-rows: 1
   :widths: 18 12 12 58

   * - Queue
     - Nodes per job
     - Walltime
     - Notes
   * - ``debug``
     - 1-4
     - 2 h
     - 1 running and 1 queued job per user, so dependent jobs (e.g. the popsyn
       merge job) can't be queued at the same time; submit them later (e.g.
       with ``posydon-popsyn check``)
   * - ``workq-route``
     - 1-176
     - 24 h
     - Routes to ``workq``: 20 jobs queued or running and 10 running per
       project; use it for production runs
   * - ``preemptable``
     - 1-10
     - 72 h
     - Up to 20 jobs per project; jobs can be **killed without warning** when
       jobs are submitted to ``demand``
   * - ``demand``
     - 1-64
     - 1 h
     - By request only

Check the current limits with ``qstat -Qf <queue>``.

Preemption and reruns
---------------------

A preempted job is killed as a whole: all its nodes and all array tasks
packed into it. With ``rerunnable=True`` PBS restarts the job script **from
the beginning**, so every task of the job is run again:

- MESA grids: set ``skip_finished_runs=True`` in the ``[run_parameters]``
  section of the grid ini. ``posydon-run-grid`` then skips grid points whose
  directory already contains a finished MESA run (one with a termination
  code), instead of deleting and rerunning them. Unfinished runs start from
  the beginning.
- Population synthesis: don't make the jobs rerunnable; after a kill,
  ``posydon-popsyn check`` creates a rescue job for the missing batches only.

Smaller jobs (fewer nodes) lose less when they are killed.


.. _large_grids:

Large MESA grids
----------------

A MESA grid larger than one job can hold runs as a series of jobs:

1. Set ``max_nodes`` to the job size you want, e.g. ``max_nodes=20`` with
   ``tasks_per_node=64`` gives 1280 slots. Without ``max_nodes``, a job
   would request a slot for every grid point.
2. Set ``skip_finished_runs=True`` in the ``[run_parameters]`` section.
3. Submit the grid job. Its slots work through the grid points in order
   until the walltime ends; runs still going are killed.
4. Submit the same job script again (e.g. ``qsub grid.pbs``). It skips the
   grid points that have finished (within seconds each) and continues with
   the killed and remaining ones. Repeat until all grid points are done.

Runs that need longer than the walltime never finish this way; run them with
a longer walltime (e.g. in ``preemptable``) or more threads per run.

Two jobs running at the same time on the same grid would not know about
each other's work and could run the same grid point twice, so submit the
next job only after the previous one ended (e.g. with
``qsub -W depend=afterany:<job ID> grid.pbs``).


.. _task_stats:

Per-task statistics
-------------------

The scheduler's own accounting can't show the resources of the individual
tasks of a packed job. With the PBS option ``task_stats=<directory>``, the
task runner records the start, end, exit code, CPU time and peak memory
(including child processes, e.g. MESA) of every task, one CSV file per rank.
Tasks that were still running when the job was killed have a start but no
end. Summarize them with:

.. code-block:: bash

    python -m posydon.CLI.schedulers.task_stats <directory> --cpus-per-task 4

which prints the numbers of finished, failed and unfinished tasks, the
distribution of run times and peak memory, and a suggested ``tasks_per_node``
for the measured peak memory.


ALCF Crux examples
------------------

Population synthesis, e.g. 1 million binaries per metallicity in 56 tasks
on one node (each task needs about 3.6 GB):

.. code-block:: bash

    posydon-popsyn setup population_params.ini --job_array 56 \
        --scheduler pbs --queue workq-route --account <PROJECT> \
        --scheduler_option select=system=crux \
        --scheduler_option filesystems=home:eagle \
        --scheduler_option tasks_per_node=56 \
        --walltime 04:00:00 --merge_walltime 06:00:00
    ./pbs_submit.sh

Run populations from a directory on Eagle rather than from your home
directory.

MESA grid (``posydon-setup-grid --submission-type pbs --inifile grid.ini``),
in the grid ini file:

.. code-block:: ini

    [pbs]
    number_of_nodes=1
    number_of_mpi_tasks=1
    number_of_cpus_per_task=2
    job_array=True
    partition=workq-route
    account={PROJECT}
    walltime='24:00:00'
    email={YOUR_EMAIL_ADDRESS}
    ; PBS options: 64 runs with 2 threads per node, 20 nodes per job
    select=system=crux
    filesystems=home:eagle
    tasks_per_node=64
    max_nodes=20
    ; optional: record run times and memory per grid point
    task_stats=stats

    [run_parameters]
    ...
    ; keep grid points that already finished (for a series of jobs, or
    ; rerunnable jobs)
    skip_finished_runs=True

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
    need longer will not finish within a single job. The ``preemptable`` queue
    allows 72 hours, with the risks described above.


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
