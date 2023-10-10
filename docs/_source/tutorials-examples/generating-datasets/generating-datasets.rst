.. _generating-datasets:

Crafting the Core Datasets for POSYDON
======================================

Discover the intricate details behind POSYDON's dataset generation, including diving deep into the vast world of MESA datasets, learning the art of downsampling, and mastering the controls through our Processing Pipeline API.

Getting Started Tutorials
-------------------------

Generating a POSYDON PSyGrid Dataset
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Dive into our Jupyter Notebook that will show you how to create a MESA PSyGrid dataset from scratch.

.. toctree::
    
    just_step_1

To learn more about the PSyGrid object or the Processing Pipeline API ini file, check out the [API Documentation](api.rst).


The fullt POSYDON Processing Pipeline experince 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

From A to Z, this tutorial shows you how to process, concatenate, downsample, plot the grids and check failure rate, train the interpolation object and export the POSYDON dataset for population synthesis.

.. toctree::

    run_full_piepeline

Congratulations! You now master the POSYDON Processing Pipeline. To learn more about the PSyGrid object or the Processing Pipeline API ini file, check out the [API Documentation](api.rst) and the indepth components overview.

Advanced Tutorials
------------------

Export MESA Simulation Points to Rerun Using the Processing Pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Learn how to export MESA simulation points to rerun using the Processing Pipeline.

.. toctree::
    
    step_rerun


Export Single Star PSyGrid Datasets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Learn how to export single star PSyGrid datasets. This is an advaced tuotorial because it requires knowledge about the EEP code.

.. node::
    
    POSYDON v2.0.0 does not embed a python interface to compute EEPs but relyes on the Fortran code of Aaron Dotter (2016). Feature POSYDON code releases will include a python interface to compute EEPs and embed the export of single stellar grid into the post-processig pipeline.

.. toctree::

    processing_single_hms

Support & Feedback
------------------

Facing challenges or have queries? Make sure to check our [FAQ](link-to-FAQ) or visit the [contact-information.rst](contact-information.rst) page.

