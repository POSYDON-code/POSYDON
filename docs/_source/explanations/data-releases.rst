.. _data-releases:

################
Data Releases
################

POSYDON shares the precomputed (down-sampled) MESA simulation grids, together with their
post-processed products (interpolators and classifiers), so that you can run
population synthesis without having to run MESA yourself. These shared datasets
are called data releases.

Data Release 1 (DR1)
--------------------

DR1 was the first published dataset accompanying POSYDON v1. It contains a
fixed grid of binary-system MESA simulations. Precise storage and content
details of DR1 are described in the original POSYDON paper, `Fragos et
al. (2023) <https://ui.adsabs.harvard.edu/abs/2023ApJS..264...45F/abstract>`_.

Because the code evolved substantially between v1 and v2, DR1 is not the default
data used by current POSYDON versions. The default download is DR2.

Data Release 2 (DR2, default)
------------------------------------

DR2 accompanies POSYDON v2 and is the default dataset used by the code. It
contains fixed grids of single- and binary-star MESA simulations processed
through the POSYDON :ref:`processing pipeline <processing-pipeline>`, together
with the trained interpolators and classifiers used for population synthesis.

Metallicities
~~~~~~~~~~~~~~

Unlike DR1, DR2 contains simulations over **eight different metallicities**:
::

    1e-04_Zsun  1e-03_Zsun  ...  (solar Zsun and above)

so that populations can be constructed over a wide range of star formation
environments. Each metallicity is a separate download, which is convenient if
you only need a subset (e.g. only the solar-metallicity grid to follow the
*Getting Started* examples).

Storage requirements
~~~~~~~~~~~~~~~~~~~~~~

The amount of disk space you need depends on the subset you download. The table
below gives the canonical figures used throughout these documentation, with a
single consistent set of numbers:

===================== =================
What you download      Approximate size
===================== =================
Full DR2 dataset       ~ 110 GB
Each DR2 metallicity   ~ 12 GB
DR1 dataset            ~ 40 GB
===================== =================

The exact footprint also depends on how many run types (``HMS-HMS``, ``He-He``,
and so on; see :ref:`grid-types`) and metallicities you include.

How to download
-----------------

The simplest way to obtain the data is with the ``get-posydon-data``
command-line tool. For example, to download just the solar-metallicity DR2
data used by the quick-start guides:

.. code-block:: bash

   get-posydon-data DR2_1Zsun

And to download the full DR2 dataset:

.. code-block:: bash

   get-posydon-data DR2
See the ``get-posydon-data`` reference page (``api_reference/get_posydon_data.rst``)
for the full list of options.

.. seealso::

   * :ref:`grid-types` for what a grid contains and how DR2 is organized.
   * :ref:`installation-guide` for how data is located on disk (the
     ``PATH_TO_POSYDON_DATA`` variable).
