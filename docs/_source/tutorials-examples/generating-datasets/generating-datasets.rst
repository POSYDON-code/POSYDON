.. _generating-datasets:

Crafting the Core Datasets for POSYDON
======================================

Discover the intricate details behind POSYDON's dataset generation, including diving deep into the vast world of MESA datasets, learning the art of downsampling, and mastering the controls through our Processing Pipeline API.

Getting Started Tutorials
-------------------------

I. Generating a POSYDON PSyGrid Dataset
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Dive into our Jupyter Notebook that will show you how to create a MESA PSyGrid dataset from scratch.

.. toctree::
    
    just_step_1

To learn more about the PSyGrid object or the Processing Pipeline API ini file, check out the :ref:`API Documentation <processing-pipeline>`.


II. The Full POSYDON Processing Pipeline Experience 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

From A to Z, this tutorial shows you how to process, concatenate, down sample, plot grids and check failure rates, train interpolation objects and export the POSYDON dataset for population synthesis.

.. toctree::

    run_full_pipeline

Congratulations! You've now mastered the POSYDON Processing Pipeline.


.. _plot_1d:

III. 1D Plotting Functionalities for POSYDON PSyGrids
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Unlock the power of 1D plotting functionalities for POSYDON PSyGrids. This tutorial shows you how to easily visualize your single and binary star tracks leveraging the `plot1D` method of the PSyGrid object.

.. toctree::

    plot_1D


.. _plot_2d:

IV. 2D Plotting Functionalities for POSYDON PSyGrids
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Unlock the power of 2D plotting functionalities for POSYDON PSyGrids. This tutorial shows you how to easily visualize your single and binary star tracks leveraging the `plot2D` method of the PSyGrid object.

.. toctree:: 

    plot_2D

Advanced Tutorials
------------------

Export MESA Simulation Points to Rerun Using the Processing Pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Learn how to export MESA simulation points to rerun using the Processing Pipeline.

.. toctree::
    
    step_rerun


Export Single Star PSyGrid Datasets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Learn how to export single star PSyGrid datasets. This is an advanced tutorial because it requires knowledge about the EEP code.

.. note::
    
    POSYDON v2.0.0 does not embed a python interface to compute EEPs but relies on the Fortran code of Aaron Dotter (2016). Future POSYDON code releases will include a python interface to compute EEPs and embed the export of single stellar grid into the post-processing pipeline.

.. toctree::

    processing_single_hms
