.. _MESA-grids:

Architecting MESA Simulation Grids with POSYDON
===============================================

You can use POSYDON framework to run your own, customized MESA simulations of single and binary systems. Here we will go through POSYDON's dedicated MESA configuration repository and infrastructure for one or even big grids of MESA runs.

The POSYDON MESA simulation parameters are stored in the `POSYDON-MESA-INLISTS <https://github.com/POSYDON-code/POSYDON-MESA-INLISTS>`_ Github repository, on its ``v2_sources`` branch.
Clone it once:

.. code-block:: bash

    git clone --branch v2_sources https://github.com/POSYDON-code/POSYDON-MESA-INLISTS.git ~/POSYDON-MESA-INLISTS


This repository contains the MESA inlists used by POSYDON to run single and binary star simulations.
At the time of writing, the POSYDON v2.0.0 code is compatible with MESA r11701 and the inlists are formatted as such.
Additionally, multiple changes were made to the MESA source code to fix bugs, including reverse mass transfer in binaries (where the secondary star fills its Roche lobe).

.. warning::

    Please contact us if you want to run MESA grids with this adapted version of MESA.

Each POSYDON release is a folder of this repository (``r11701/releases/v2.0/``), and the reruns of a release, the physics variants some of its grids were run again with, are folders inside it (``reruns/<rerun>/``). You choose a release and a rerun in your :samp:`ini` file; see :ref:`posydon_grid`.

To follow the next step of the tutorial, you will need to have a MESA version compatible with POSYDON installed on your machine. If you do not have MESA installed, please follow the instructions on the `MESA website <https://docs.mesastar.org/en/release-r23.05.1/>`_.

.. warning::

    The POSYDON v2.0.0 code is compatible with MESA r11701 inlists. Support might not be available for the latest MacOS version.


How-to guides
-------------

Each guide sets up one kind of grid with ``posydon-grid``, from an empty
folder to a run directory ready to submit. They were run as shown, with the
outputs kept.

.. toctree::
   :maxdepth: 1

   how-to/hms_hms_rerun
   how-to/single_hms_rerun
   how-to/single_hems_rerun
   how-to/stack_own_inlist
   how-to/shell_laptop

Ready-made :samp:`ini` files for every grid type are in the `examples folder
<https://github.com/POSYDON-code/POSYDON/tree/main/docs/_source/components-overview/mesa_grids/examples>`_.

Getting Started Tutorials (legacy ini format)
---------------------------------------------

.. note::

    The tutorials below use the legacy :samp:`ini` format (``grid_params.ini``,
    with ``[slurm]``, ``[mesa_inlists]``, ``[mesa_extras]`` and
    ``[run_parameters]``), which ``posydon-grid`` still reads. For new grids,
    follow the how-to guides above.

I. Running your first MESA simulation using POSYDON
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following Jupyter Notebook guides you through our dedicated MESA submission API, showcasing how to simulate a HMS-HMS binary systems.

The first version of this tutorial showcases how to run the simulation using SLURM on HPC.

.. toctree::

    1_hms_hms

In case you want to run the simulation locally, you can follow the second version of this tutorial.

.. toctree::

    laptop

To gather more information about the MESA simulation submission API ini file, please refer to the :ref:`MESA Grids API <mesa-grids-api>` section of the documentation.

II. Running your first MESA grid simulation using POSYDON
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following tutorial will show you how to run a grid of MESA simulations using POSYDON. The grid will consist of 100 simulations of HMS-HMS binary systems at 0.1Zsun metallicity in the mass ratio slice 0.7 each with different primary masses and orbital periods. The simulations will be run on a HPC using SLURM.

.. toctree::

    running_a_grid

Now that you have run your first grid, you can process, visualize and explore
the results using the POSYDON :ref:`post-processing pipeline
<processing-pipeline>`.

III. Running single stars with POSYDON
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This tutorials shows how to run a grid of single star simulations using POSYDON.

.. toctree::

   running_single_hms

We show how to export EEPs and make a PSyGrid object in this tutorial
(see the :ref:`PSyGrid <psygrid>` and :ref:`post-processing pipeline
<processing-pipeline>` pages for how the extracted data is organized).

Advanced Tutorials
------------------

Creating and Customizing MESA Grids layer by layer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Every inlist is stacked from MESA's defaults, the POSYDON release (and
rerun), and your own inlists. See the how-to
:doc:`how-to/stack_own_inlist`, and :ref:`posydon_grid` for the layers.

Running a dynamically sampled grid of MESA simulations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. warning::

    This feature is experimental. Please contact us if you encounter any issues.

Provided a sparse, rectilinearly sampled MESA grid, POSYDON allows to dynamically sample the parameter space to enhance the coverage of the parameter space to improve classification and interpolation accuracy.


Support & Feedback
------------------

Encountered obstacles or have questions? Ensure you consult our
:ref:`code usage FAQ <code-usage>` or drop by the
:ref:`contact information <contact_info>` page.
