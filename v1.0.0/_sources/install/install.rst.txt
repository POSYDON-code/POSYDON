.. _install:

############
Installation
############

=========================================
Installing the latest release from GitHub
=========================================


Creating a conda environment
----------------------------

We recommend to install the Python distribution Anaconda and to install POSYDON
in a virtual environment. After creating the environment (you can choose any
name, e.g., `posydon`, or `posydon_env`) like this:

.. code-block::

    conda create -n posydon python=3.7

Make sure you agree to any questions about installing required libraries. To
proceed with the installation, we will need to activate the environment:

.. code-block::

    conda activate posydon

Cloning POSYDON
---------------
Clone the repository in a local directory, e.g. `/home/POSYDON/`, with

.. code-block::

    git clone https://github.com/POSYDON-code/POSYDON.git


The directory will contain the following structure:

.. code-block::

    posydon/
    README.md
    setup.py
    etc.

Installing POSYDON
------------------
Exporting the global path
~~~~~~~~~~~~~~~~~~~~~~~~~
Export the path to the cloned POSYDON code (you can add this line to your
.bashrc or .bash_profile.), e.g.

.. code-block::

    export PATH_TO_POSYDON=/home/POSYDON/

Installing the package
~~~~~~~~~~~~~~~~~~~~~~
From the cloned POSYDON directory execute the commands to install POSYDON and
the `mpi4py` dependency

.. code-block::

    pip install -e .
    conda install mpi4py


Downloading POSYDON data
------------------------
OPTION 1: Zenodo (latest official version)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Export the path to where you want to clone the data, e.g. `/home/`, and
download the data with the following commands

.. code-block::

    export PATH_TO_POSYDON_DATA=/home/
    get-posydon-data


OPTION 2: git LFS submodule (latest development version)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Inside the cloned POSYDON repository run the following commands to
install and initialise git LFS, and download the data.

.. code-block::

    export PATH_TO_POSYDON_DATA=$PATH_TO_POSYDON/data/
    conda install git-lfs
    git lfs install
    git submodule init
    git submodule update data/POSYDON_data


=========================================
Installing POSYDON documentations modules
=========================================

These modules are needed in order to compile the documentation

.. code-block::

    pip install -e .[doc]

To compile the documentation and open the html page do the following

.. code-block::

    cd docs/
    make html
    open _build/html/index.html


Installation Notes/FAQ
----------------------

.. note::

    USING IPYTHON OR JUPYTER-NOTEBOOKS WITH POSYDON ENVIRONMENT

    Please note that using the global instance of the conda jupyter-notebook
    or ipython will most likely fail when trying to use posydon.
    PLEASE explicitly install both into the posydon environment with either

    ``conda install jupyter ipython``

    ``pip install jupyter ipython``
