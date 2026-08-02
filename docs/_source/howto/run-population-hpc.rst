.. _run-hpc:

Running on an HPC
=================

A large population run is typically launched on an HPC facility with SLURM.
Each metallicity is handled by one job, and the binaries of that metallicity
are split across a job array. Getting the array size and the walltime right is
the key to running efficiently without either wasting resources or timing out.

.. note::

   A full worked example of an HPC population is given in ``pop_syn.ipynb``
   under ``tutorials-examples/population-synthesis/``.

Choosing the job-array size
---------------------------

Size the array so that each array element handles at least a few thousand
binaries. Every job has a fixed startup cost because it has to load the grids,
so an array that is too fine spends most of its time loading rather than
computing. As a rule of thumb, plan for at least about 1000 binaries per job.

Choosing the walltime
---------------------

Each binary takes roughly 1-2 seconds to evolve. A good starting estimate is::

    walltime  ~  (seconds per binary) x (binaries per job)

For example, 100,000 binaries split over 100 jobs means each job evolves 1,000
binaries, which takes about 33 minutes, so a walltime of about 45 minutes
(``00:45:00``) is reasonable.

Balancing array size and walltime
----------------------------------

The two parameters are linked:

* if the walltime is too long, increase the array size so each job finishes
  sooner and the overall population completes faster;
* if the walltime is too short, decrease the array size, because every job has
  a fixed setup cost independent of how many binaries it runs.

A minimal SLURM script looks like:

.. code-block:: bash

   #!/bin/bash
   #SBATCH --nodes=1
   #SBATCH --array=0-99
   #SBATCH --time=00:45:00
   #SBATCH --mem-per-cpu=10G

   source activate posydon-debug
   export PATH_TO_POSYDON=/path/to/posydon
   export PATH_TO_POSYDON_DATA=/path/to/posydon_data
   python run_population.py $SLURM_ARRAY_TASK_ID

Memory also matters per CPU, so read :ref:`memory-tuning` before finalizing the
``--mem-per-cpu`` request.

Related reading
---------------

* :ref:`memory-tuning` - how much RAM to request and ``dump_rate`` guidance.
* :ref:`fix-reruns` - what to do about binaries that could not be evolved.
* The :ref:`population parameters <pop-params-guide>` reference - how a
  population run is configured through an ``ini`` file.
