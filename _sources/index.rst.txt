.. _home:

POSYDON
===================================

**POSYDON** combines functionality to run large grids of MESA binary
simulations, efficiently store those simulations, interpolation over those
grids to synthesize binary populations, along with a series of visualization
routines to interpret the results.

This documentation is designed to provide an introduction to using POSYDON,
which is described in a `paper <https://arxiv.org/abs/2202.05892>`_ currently
under review. There are many excellent textbooks for learning about binary
evolution.

To report a bug or download the bleeding edge version of POSYDON, please visit
our `GitHub page <https://github.com/POSYDON-code/POSYDON>`_.


*****************
Table of Contents
*****************

.. toctree::
   :maxdepth: 1
   :caption: User guide

   install/install
   data/data
   api

.. toctree::
   :maxdepth: 1
   :caption: Binary Star Grids

   run_mesa_grids/inifile
   run_mesa_grids/fixed/fixed
   PSyGrid/PSyGrid
..
   run_mesa_grids/dynamic/dynamic

.. toctree::
   :maxdepth: 1
   :caption: Generating Populations

   pop_synth/POSYDON_populations
   pop_synth/pop_synth
   pop_synth/custom_flow
   SingleStar/SingleStar
   BinaryStar/BinaryStar
   advanced_pop_options/custom_hooks
   advanced_pop_options/debugging_binaries
   interpolation/IF-Interpolator/IFInterpolator

.. toctree::
   :maxdepth: 1
   :caption: Plotting

   visualization/plot1D/plot1D
   visualization/plot2D/plot2D
   visualization/VHD/VHD
   interpolation/violinplot

.. toctree::
   :maxdepth: 1
   :caption: About
   :hidden:

   about/licence
   about/attribution


License & Attribution
---------------------

Copyright (c) 2022, Tassos Fragos

If you make use of POSYDON in your work, please cite our paper
(`arXiv <https://arxiv.org/abs/2202.05892>`_,
`ADS <https://ui.adsabs.harvard.edu/abs/2022arXiv220205892F/abstract>`_).
