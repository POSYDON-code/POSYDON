.. _MESA-grids:

Architecting MESA Simulation Grids with POSYDON
===============================================

Boost your MESA simulation grid prowess by diving into the intricacies of POSYDON's grid creation techniques and the dedicated MESA configuration repository. Harness the power of our API to craft precise and customized grids.

Getting Started Tutorials
-------------------------

The POSYDON MESA simulation parameters are stored in `POSYDON-MESA-INLISTS <https://github.com/POSYDON-code/POSYDON-MESA-INLISTS>`_ GitHub submodule. To obtain the submodule, clone the POSYDON repository with the following command in your terminal:

.. code-block:: bash

    cd $PATH_TO_POSYDON
    git submodule init
    git submodule update grid_params/POSYDON-MESA-INLISTS/

Similarly to the POSYDON code repository, there is a `main` branch that points to the latest stable version of the POSYDON-MESA-INLISTS repository associated to a POSYDON code release and a `development` branch that contains the latest stable version of our fiducial MESA configurations, run the following command to get the lates development version of the sumodule:

.. code-block:: bash

    cd grid_params/POSYDON-MESA-INLISTS/
    git checkout development

To follow the next step of the tutorial, you will need to have MESA installed on your machine. If you do not have MESA installed, please follow the instructions on the `MESA website <https://docs.mesastar.org/en/release-r23.05.1/>`_.

.. warning::

    The POSYDON v2.0.0 code is compatible with MESA r11701. Support might not be available for the latest MacOS version.


Running your first MESA simulation using POSYDON
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Immerse yourself in our Jupyter Notebook which guides you through our dedicated MESA submission API, showcasing how to simulate a HMS-HMS binary systems.

The first version of this tutorial showcases how to run the simulation using SLURM on HPC. 

.. toctree::
    
    1_hms_hms

In case you want to run the simulation locally, you can follow the second version of this tutorial.

.. toctree::
   
    laptop

To gather more information about the MESA simulation submission API ini file, please refer to the :ref:`MESA Grids API <mesa-grids-api>` section of the documentation.

Running your first MESA grid simulation using POSYDON
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following tutorial will show you how to run a grid of MESA simulations using POSYDON. The grid will consist of 100 simulations of HMS-HMS binary systems at 0.1Zsun metallicity in the mass ratio slice 0.7 each with different primary masses and orbital periods. The simulations will be run on a HPC using SLURM.

.. toctree::
    
    running_a_grid

Now that you have run your first grid, you can process, visualize and explore the results using the POSYDON post-processing pipeline API. TODO: link

Rimulating single stars with POSYDON
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This tutorials shows how to run a grid of single star simulations using POSYDON.

.. toctree::
   
   running_single_hms

We show how to export EEPs and make a PSyGrid object in this tutorial TODO: link the advanced tutorial of step 2.

Advanced Tutorials
------------------

Creating and Customizing MESA Grids layer by layer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Step-by-step guide on how to leverage our MESA simulation submission API for creating and customizing MESA grids to suit your needs.

Dive deeper into advanced techniques for grid creation, by layering MESA and POSYDON default inlist parameters and customizing your inlists variables on top to your needs.

.. toctree::
    
    notebook: Creating and Customizing MESA Grids layer by layer

Running a dynamically sampled grid of MESA simulations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. warning::

    This feature is experimental. Please contact us if you encounter any issues.

Provided a sparsed, rectilinearly sampled MESA grid, POSYDON allows to dynamically sample the parameter space to henence the coverage of the paramter space to imporve classification and interpolation accuracy.

.. toctree::
   
   dynamic


Support & Feedback
------------------

Encountered obstacles or have questions? Ensure you consult our [FAQ](link-to-FAQ) or drop by the [contact-information.rst](contact-information.rst) page.

