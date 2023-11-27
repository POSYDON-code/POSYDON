.. _install:

############
Installation
############

The easiest way to start with POSYDON is to install it from Anaconda
(Automatic Installation), and call an executable script that downloads the grid
data automatically (Automatic Data Download).

Alternatively, if you are interested in modifying the code, you can install
POSYDON manually through `github` (Manual Installation). You can still use the
automatic script to download the data, however, if you are about to modify the
grid data, you can get them using the `git-lfs` module (Manual Data Download).

Below we describe both ways.

==========================================================
Recommended way (Automatic installation and Data Download)
==========================================================

Installing POSYDON v1.0.0 from Anaconda (Automatic Installation)
----------------------------------------------------------------

We recommend to install the Python distribution Anaconda and to install POSYDON
in a virtual environment. Specifically, we recommend using the installation we
have provided which can be accessed through conda-forge. On Linux, the new
conda environment can be created (we have named our environment posydon-example,
but you can choose any name), the conda-forge channel added, and the required
library installation can all be completed in one line:

.. code-block::

    conda create --name posydon-example -c posydon -c conda-forge posydon

On Mac (or if you have problems with the above command on Linux), these steps
likely need to be separately run:

.. code-block::

    conda create -n posydon-conda python=3.7
    conda activate posydon-conda
    conda config --add channels conda-forge
    conda config --add channels posydon
    conda config --set channel_priority false
    conda install posydon

Now, you can source the environment with

.. code-block::

    conda activate posydon


Getting the POSYDON data (Automatic Data Download)
--------------------------------------------------
Export the path to where you want to clone the data, e.g. `/home/`, and
download the data from ZENODO with the following commands

.. code-block::

    export PATH_TO_POSYDON_DATA=/home/
    get-posydon-data

(Note that you can add the export command to your .bashrc or .bash_profile.)



==========================================================
Installing POSYDON v.1.0 from GitHub (Manual Installation)
==========================================================


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

    git checkout main
    pip install -e .
    conda install mpi4py


Downloading POSYDON data (after Manual Installation via `github`)
-----------------------------------------------------------------
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


Installing POSYDON documentations modules
-----------------------------------------

In the case of manual installation you can also alter and build the
documentation. These modules are needed in order to compile the documentation

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
