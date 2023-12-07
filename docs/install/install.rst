.. _install:

############
Installation
############

The easiest way to start with POSYDON is to install it from Anaconda
(Recommended way), and call an executable script that downloads the grid
data automatically.

Alternatively, if you are interested in modifying the code, you can install
POSYDON manually through `github` (Manual Installation). Below we describe both
ways.

=========================================================
Installing POSYDON v1.0.0 from Anaconda (Recommended way)
=========================================================

Installation
------------

We recommend using Anaconda and to install POSYDON in a virtual environment.
We have created a package which can be accessed through conda-forge. On Linux,
the new conda environment can be created (we have named our environment
posydon-example, but you can choose any name), the conda-forge channel added,
and the required library installation can all be completed in one line using
the terminal or command line:

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

Now, you can activate the environment with

.. code-block::

    conda activate posydon


Downloading the POSYDON data
----------------------------
Because the data is large, ~10 GB, it must be downloaded
with an explicit command. Export the path to where you want
to clone the data, e.g. `/home/`, and download the data from
ZENODO with the following commands

.. code-block::

    export PATH_TO_POSYDON_DATA=/home/
    get-posydon-data

(Note that you can add the export command to your .bashrc or .bash_profile.)



==========================================================
Installing POSYDON v.1.0 from GitHub (Manual Installation)
==========================================================

Creating a conda environment
----------------------------

As above, we recommend using Anaconda to install POSYDON in a virtual
environment. After creating the environment (you can choose any name, e.g.,
`posydon`, or `posydon_env`) like this:

.. code-block::

    conda create -n posydon python=3.7

Make sure you agree to any questions about installing required libraries. To
proceed with the installation, you will need to activate the environment:

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
the `mpi4py` dependency (first make sure you are on the main branch!)

.. code-block::

    git checkout main
    pip install -e .
    conda install mpi4py


Getting the POSYDON data
------------------------
Export the path to where you want to clone the data, e.g. `/home/`, and
download the data from ZENODO with the following commands

.. code-block::

    export PATH_TO_POSYDON_DATA=/home/
    get-posydon-data

(Note that you can add the export command to your .bashrc or .bash_profile.)


Installing POSYDON documentation modules
----------------------------------------

In the case of manual installation you can also alter and build the
documentation. These modules are needed in order to compile the documentation

.. code-block::

    pip install -e .[doc]

To compile the documentation and open the html page use the following commands

.. code-block::

    cd docs/
    make html
    open _build/html/index.html


======================
Installation Notes/FAQ
======================

.. note::

    USING IPYTHON OR JUPYTER-NOTEBOOKS WITH POSYDON ENVIRONMENT

    Please note that using the global instance of the conda jupyter-notebook
    or ipython will most likely fail when trying to use posydon.
    PLEASE explicitly install both into the posydon environment with either

    ``conda install jupyter ipython``

    ``pip install jupyter ipython``
