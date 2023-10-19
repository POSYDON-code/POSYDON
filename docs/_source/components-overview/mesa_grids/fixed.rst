.. _fixed_grid:

#####################
Run a fixed MESA grid
#####################

Default evolutionary parameters
===============================

POSYDON makes it easy to run a fixed grid of binaries with MESA within a high
performance computing environment. Note that these instructions assume you
are using slurm as your job scheduler. We just need to start by setting a few
environmental variables. We recommend you add these to your `.bashrc` or
`.bash_profile` so they are included at login. The filepaths of these need to be
edited for your own installation.

.. code-block::

    export MESA_DIR="/projects/b1119/mesa_sdk/mesa-r11701"
    export MESASDK_ROOT="/projects/b1119/mesa_sdk/mesasdk"
    source $MESASDK_ROOT/bin/mesasdk_init.sh

Next, you will want to create a working directory for your grid of binaries.
For this example, we will create a working directory called `example_grid`:

.. code-block::

    mkdir example_grid
    cd example_grid

This directory will contain all the binary runs and associated data generated
by the individual MESA runs as well as the functionality for slurm to run the
grid.

Now that we have our grid, we need to create two separate files: a file
containing the list of binaries we want to evolve and an .ini file that
contains all the parameters associated with running the grid. Let's start with
the list of binaries we want to evolve. For the case of this example, let's
evolve a series of high-mass X-ray binaries with neutron star accretors at a
series of initial binary periods. We'll use our favorite text editor to add the
following lines into a file we'll call `grid.csv`. As we'll see below, it can
be given any name, but the entries must be separated by commas.

.. code-block::

    m1,m2,initial_period_in_days
    10,1.4,1
    10,1.4,10
    10,1.4,100
    10,1.4,1000

POSYDON internally allows any variable to be included here, so long as it has
the same name as a variable in MESA. However, these three parameters ought to
always be specified, otherwise you will be evolving a binary with a default
orbital period and masses.

As for the .ini file, you can copy over the `example
<https://github.com/POSYDON-code/POSYDON/blob/development/grid_params/grid_params.ini>`_
from the :ref:`inifile` page, which provides a more detailed description of all
the entries. Make sure to carefully go through each of the separate entries and
adjust them for your particular needs. Place that .ini file (which we have
named `example_grid.ini`) in our example_grid directory.

Now, with our `grid.csv` and `example_grid.ini` files ready to go, we use a
POSYDON script `posydon-setup-grid` to generate all the necessary files to run
the grid using slurm. Note the different arguments here. We are designating
that: 1) we are using a fixed grid, 2) that our job scheduler is slurm, and 3)
the name of the .ini file. We may expand compatibility for other job
schedulers in the future, but for now slurm is the only implemented option.

.. code-block::

    posydon-setup-grid --grid-type fixed --submission-type slurm --inifile example_grid.ini

This will take a minute to run. You will note that it includes a pull from the
POSYDON GitHub repository to ensure the correct, designated version of the MESA
inlists are being used. It is possible that you will need to enter your
GitHub username and password during this step. Additionally, note that this
step takes some time to run, as several MESA files are being compiled.

Once finished, you will note that in our `example_grid` directory, we now have
several new files and new directories. One advantage of this implementation is
that all the binaries we run use a single version of the MESA executables,
rather than each binary having its own copy of the same files.

Finally, we are ready to submit our jobs with:

.. code-block::

    sbatch slurm_job_array_grid_submit.sh

Once the grid of runs is completed, we recommend you use our provided PSyGrid
functionality to interpret and collate the individual binary runs
(link for documentation: :mod:`posydon.grids.psygrid`.)


Non-default evolutionary parameters
===================================

If you want to change the binary evolution parameters so you are using
non-default options, we have constructed a hierarchy of MESA inlists. You can
provide a non-default option in your own user-provided MESA inlist, which is
explicitly linked in the inifile (make sure you uncomment the appropriate
lines). We additionally provide the capability to use your own
`run_star_extras.f` and `run_binary_extras.f` as well as provide your own list
of history and profile columns. See :ref:`inifile` for details.
