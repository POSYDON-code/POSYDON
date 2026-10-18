.. _installation-guide:

Installation Guide
******************

.. contents:: On this page
   :local:
   :depth: 1

Before you start
================

- **Operating system**: Linux (officially tested, including HPC clusters) or
  macOS (supported on a best-effort basis).
- **conda >= 23.1.0** (Anaconda, Miniconda or Miniforge). Older versions can
  take hours to resolve POSYDON's dependencies; see
  :ref:`installation-issues` if ``conda --version`` reports something older.
- **Python 3.11**, which is installed into the conda environment for you.
- **Disk space** for the data you download: about 12 GB for a single
  metallicity, about 110 GB for the full DR2 dataset (see
  :ref:`data-releases`).
- **Memory**: at least 8 GB of RAM, 16 GB recommended. Population runs are
  parallelized with one process per CPU, so a multi-core machine or HPC
  cluster helps.

Quick install
=============

The installation script creates a conda environment, installs POSYDON, sets
the required environment variables and asks which data to download:

.. code-block:: bash

    bash <(curl -sSL https://raw.githubusercontent.com/POSYDON-code/POSYDON/main/install.sh)

Then continue with :ref:`verify-installation`. If you prefer to do each step
yourself, follow one of the manual installations below.

Manual install with conda
=========================

1. **Create and activate an environment**

   .. code-block:: bash

       conda create --name posydon_env python=3.11 -y
       conda activate posydon_env
       conda config --add channels conda-forge

2. **Install POSYDON**

   .. code-block:: bash

       conda install -c posydon posydon

3. **Set the environment variables**

   ``PATH_TO_POSYDON`` points to the directory that contains the installed
   ``posydon`` package, and ``PATH_TO_POSYDON_DATA`` is where the data will be
   stored. Choose a data directory with enough free space; keeping it
   separate from the code is recommended.

   .. code-block:: bash

       conda env config vars set PATH_TO_POSYDON=$(python -c "import os, posydon; print(os.path.dirname(os.path.dirname(posydon.__file__)))")
       conda env config vars set PATH_TO_POSYDON_DATA=/path/where/you/want/to/store/data
       conda activate posydon_env

   ``conda env config vars`` stores the variables in the environment, so they
   are set automatically every time you run ``conda activate posydon_env``.
   Re-activating the environment (the last line) applies them now.

4. **Download the data**: see :ref:`download-data`.

Development install from GitHub
===============================

Use this if you want the latest development version, or plan to modify the
code or :ref:`contribute <how-to-contribute>`.

1. **Create and activate an environment**, as in step 1 above.

2. **Clone the repository and install it in editable mode**

   .. code-block:: bash

       git clone https://github.com/POSYDON-code/POSYDON.git
       cd POSYDON
       pip install -e .

   The repository is cloned into the current directory, so move to where you
   want it first. On Apple Silicon Macs, see
   :ref:`the macOS notes <installation-issues-macos>` if ``pip`` fails to
   build ``tables``.

3. **Set the environment variables**, from inside the ``POSYDON`` directory:

   .. code-block:: bash

       conda env config vars set PATH_TO_POSYDON=$(pwd)
       conda env config vars set PATH_TO_POSYDON_DATA=/path/where/you/want/to/store/data
       conda activate posydon_env

   If you are not using conda, ``export`` the same two variables in your
   ``~/.bashrc`` (or your shell's equivalent) instead.

4. **Download the data**: see :ref:`download-data`.

.. _download-data:

Download the data
=================

POSYDON needs precomputed MESA grids, interpolators and auxiliary files to run.
The ``get-posydon-data`` command downloads them into ``PATH_TO_POSYDON_DATA``.

For the :ref:`getting-started examples <first-population>` and most tutorials,
the solar-metallicity set (about 12 GB, including the auxiliary data) is
enough:

.. code-block:: bash

    get-posydon-data DR2_1Zsun

To download all eight DR2 metallicities (about 110 GB):

.. code-block:: bash

    get-posydon-data DR2

Run ``get-posydon-data -l complete`` to list all datasets, and
``get-posydon-data -h`` for the other options. The data can also be downloaded
manually from the `POSYDON community on Zenodo
<https://zenodo.org/communities/posydon/>`_. See :ref:`data-releases` for what
each release contains.

.. _verify-installation:

Verify the installation
=======================

With the environment activated, check that POSYDON imports and both
environment variables are set:

.. code-block:: bash

    python -c "import posydon; print(posydon.__version__)"
    echo $PATH_TO_POSYDON $PATH_TO_POSYDON_DATA

You are now ready for :ref:`your first population <first-population>`. If
something went wrong, see :ref:`installation-issues`.

Optional extras
===============

Some features need additional packages. From a GitHub checkout, install them
with ``pip install ".[<extra>]"`` in the ``POSYDON`` directory; with a conda
install, install the listed package directly.

.. list-table::
   :header-rows: 1
   :widths: 10 30 60

   * - Extra
     - Package
     - Needed for
   * - ``hpc``
     - ``mpi4py``
     - Running MESA grids or parallel population runs on HPC clusters. We
       recommend ``conda install mpi4py``, as installing it with ``pip`` often
       fails.
   * - ``ml``
     - ``tensorflow``
     - Active learning and profile interpolation. Not needed for the provided
       initial-final interpolators.
   * - ``vis``
     - ``PyQt5``
     - Experimental visualization tools.
   * - ``doc``
     - Sphinx and related packages
     - :ref:`Building the documentation locally <building-docs>`.

Running the tutorials in Jupyter
--------------------------------

The tutorials are Jupyter notebooks. To run them, install JupyterLab in your
environment and register the environment as a kernel:

.. code-block:: bash

    conda install -c conda-forge jupyterlab ipykernel
    python -m ipykernel install --user --name=posydon_env
    jupyter lab
