.. _installation-guide:

Installation Guide
******************

This guide will help you install POSYDON step-by-step. We recommend using Anaconda, a package manager, to manage the installation and dependencies, ensuring a smooth setup process.

.. contents:: Table of Contents
   :local:

Installing POSYDON
==================

Anaconda (Recommended)
----------------------

1. **Install Anaconda**

    If you haven't already, download and install Anaconda from `Anaconda's official website <https://www.anaconda.com/products/distribution>`_.

2. **Create a New Conda Environment**

    Open your terminal or Anaconda prompt and create a new environment called ``posydon_env`` using Python 3.11 (the ``-y`` automatically answers all confirmations with yes):

    .. code-block:: bash

        conda create --name posydon_env python=3.11 -y

    Activate the new environment:

    .. code-block:: bash

        conda activate posydon_env

    As many of the required versions of packages are available on ``conda-forge`` you may need to add this channel to conda if you have not already:

    .. code-block:: bash

        conda config --add channels conda-forge

3. **Install POSYDON**

    With the environment activated, install the POSYDON code from the POSYDON channel using (by default the installation will take place in the current directory, hence please navigate to your desired location first or use the ``-p`` option to specify a path):

    .. code-block:: bash

        conda install -c posydon posydon

    .. note:: 
        In the next step you will need to know where this package was installed. You can find it by running the following 
        commands in Python:

        .. code-block:: python

            import posydon
            print(posydon.__file__)

.. _posydon-env:

4. **Set Environment Variables**

    Export the required paths (please change the location names accordingly to your installation from the previous step):

    .. code-block:: bash

        conda env config vars set PATH_TO_POSYDON=/path/to/your/posydon/installation
        conda env config vars set PATH_TO_POSYDON_DATA=/path/where/you/want/to/store/data

    The path for the data location is up to you, but keeping the data separate from the code is recommended for organization.

    .. note:: 
        So that your new `conda` environment is loaded whenever you open a new terminal, you can add the following line to 
        your :code:`~/.bashrc` or :code:`~/.bash_profile` or your shell equivalent:
        
        .. code-block:: bash
            
            conda activate posydon_env

5. **Download the Dataset**

    You can use POSYDON's built-in API command (the downloaded data will be downloaded to the directory specified by :code:`PATH_TO_POSYDON_DATA`):

    .. warning::
        Executing the following command ``get-posydon-data`` will download the full DR2 data set. This includes data for 
        all eight of the available metallicities, plus auxillary data. This may be more than you want, the data for each 
        metallicity requires about 10 GB of disk space. In total, the DR2  dataset requires 103 GB of disk space.

    .. code-block:: bash

        get-posydon-data

    You may use :code:`get-posydon-data -h` to see all the options for this command, which allows you to list all the datasets and download the ones of your choice.
    If you would rather download something more minimal to run populations, you will need at least the auxillary data:

    .. code-block:: bash

        get-posydon-data auxiliary

    and data for at least one metallicity (the examples and tutorials typically use solar):

    .. code-block:: bash

        get-posydon-data DR2_1Zsun

    Alternatively, you can manually download the datasets from Zenodo. You can find the POSYDON datasets on the `POSYDON community <https://zenodo.org/communities/posydon/>`_ on Zenodo.

.. _stable-version:

GitHub (Stable Version)
----------------------------

For users interested in the latest features and developments, you can install POSYDON directly from its GitHub repository:

1. **Clone the Repository**

    In your terminal or command prompt execute one of the following command to clone the repo with the ``https`` protocol:

    .. warning:: 
        By default, the repository will be placed in the current directory, so navigate to your desired location before proceeding.

    .. code-block:: bash

        git clone https://github.com/POSYDON-code/POSYDON.git

    For the ``ssh`` protocol:

    .. code-block:: bash

        git clone git@github.com:POSYDON-code/POSYDON.git

    Or for the Github CLI:

    .. code-block:: bash

        gh repo clone POSYDON-code/POSYDON

2. **Install the Stable Version**

    .. warning::
        If you are installing POSYDON on a Mac with Apple M1 or M2 chips, you should first install ``hdf5`` and ``pytables`` through ``conda``
        with ``conda install hdf5 pytables``, before following the instructions below.

    Navigate to the cloned repository's directory:

    .. code-block:: bash

        cd POSYDON

    Install the software as an editable package using ``pip``:

    .. code-block:: bash

        pip install -e .

3. **Set Environment Variables**

    Next, export the required paths as environment variables. From the ``POSYDON`` directory that you just performed the installation in, 
    you can execute ``pwd`` if you are unsure of the path name. Please change the location names accordingly below to your installation 
    path:

    .. code-block:: bash

        export PATH_TO_POSYDON=/path/to/your/posydon/installation
        export PATH_TO_POSYDON_DATA=/path/where/you/want/to/store/data

    The path for the data location is up to you, but keeping the data separate 
    from the code is recommended for organization.

    .. note:: 
        You can add these lines to your :code:`~/.bashrc` or :code:`~/.bash_profile` or your shell equivalent to ensure the environment variables are set every time you open a new terminal.

4. **Download the Dataset**

    You can use POSYDON's built-in API command (the downloaded data will be downloaded to the directory specified by :code:`PATH_TO_POSYDON_DATA`):

    .. warning::
        Executing the following command ``get-posydon-data`` will download the full DR2 data set. This includes data for 
        all eight of the available metallicities, plus auxillary data. This may be more than you want, the data for each 
        metallicity requires about 10 GB of disk space. In total, the DR2  dataset requires 103 GB of disk space.

    .. code-block:: bash

        get-posydon-data

    You may use :code:`get-posydon-data -h` to see all the options for this command, which allows you to list all the datasets and download the ones of your choice.
    If you would rather download something more minimal to run populations, you will need at least the auxillary data:

    .. code-block:: bash

        get-posydon-data auxillary

    and data for at least one metallicity (the examples and tutorials typically use solar):

    .. code-block:: bash

        get-posydon-data DR2_1Zsun

    Alternatively, you can manually download the datasets from Zenodo. You can find the POSYDON datasets on the `POSYDON community <https://zenodo.org/communities/posydon/>`_ on Zenodo.

Installing Jupyter for Tutorials
=================================

Our tutorials are provided as Jupyter notebooks. If you want to run these notebooks interactively, you will need to have either Jupyter Lab or Jupyter Notebook installed.

1. **Using Anaconda (Recommended)**


    If you have already installed Anaconda as suggested earlier in the installation guide, installing Jupyter Lab or Notebook is straightforward:

    .. code-block:: bash

        conda install -c conda-forge jupyterlab

    Or, for the classic Jupyter Notebook:

    .. code-block:: bash

        conda install -c conda-forge notebook

2. **Alternatively, via pip**


    If you prefer using ``pip``, you can also install Jupyter Lab or Notebook using the following commands:

    .. code-block:: bash

        pip install jupyterlab

    Or, for the classic Jupyter Notebook:

    .. code-block:: bash

        pip install notebook

3. **After Installation**


    Once installed, you can start Jupyter Lab or Notebook by running:

    .. code-block:: bash

        jupyter lab

    Or:

    .. code-block:: bash

        jupyter notebook

    From the terminal or command prompt. This will open a browser window where you can navigate to the downloaded notebooks and run them interactively.

    .. note::
        Remember to navigate to the directory containing the Jupyter notebooks or you won't see them listed in the Jupyter interface.



Installing additional dependencies (Optional) 
=============================================

For some specific functionalities, you may need to install additional dependencies.
Below are the instructions for installing these dependencies and what they are used for.

Running grids of MESA models using POSYDON
------------------------------------------

If you are planning to create MESA grids using POSYDON on HPC facilities, it's essential to have ``mpi4py`` installed to take advantage of parallel computations.
You do not need to have ``mpi4py`` installed if you are only running population synthesis simulations.

1. **Install mpi4py via Anaconda (Recommended)**:

    .. code-block:: bash

        conda install mpi4py

2. **Alternatively, via pip**:

    .. code-block:: bash

        pip install ".[hpc]"


.. warning::
    Users have reported issues when trying to install ``mpi4py`` via pip. If you encounter any issues, try installing ``mpi4py`` through 
    Anaconda. If you cannot solve the issue, please refer to the :ref:`Troubleshooting Guide <installation-issues>` or seek support from 
    the community or developers, see the :ref:`contact us <contact_info>` page.


Documentation Building
------------------------

If you're interested in building the POSYDON documentation locally:

1. **Install Documentation Modules**:

    Navigate to your POSYDON directory and install the required documentation modules:

    .. code-block:: bash

        pip install ".[doc]"

    .. warning::
        If you are installing POSYDON on a Mac with Apple M1 or M2 chips, you should install ``pandoc`` separately through brew with 
        ``brew install pandoc``.

2. **Compile the Documentation**:

    Once you have the required modules installed, you can build the documentation using Sphinx:

    .. code-block:: bash

        cd docs
        make html

    This command will generate the HTML documentation in the ``_build/html`` directory within the ``docs`` folder.

3. **Open the Compiled Documentation**:

    After successfully building the documentation, you can view it in your preferred browser. Navigate to the build directory and open the ``index.html``:

    .. code-block:: bash

        open _build/html/index.html

    .. note::
        The ``open`` command works on macOS. If you're using a different OS, you might need to open the ``index.html`` using your file manager or use a different command.



Machine Learning Dependencies
---------------------------------------

For users who wish to utilize POSYDON's latest machine learning features. 
This is specifically used in the active learning module and profile interpolation.
You do not require these dependencies if you are using the provided Initial-Final interpolators.

1. **Navigate to your POSYDON directory** (where the ``setup.py`` is located) and run:

    .. code-block:: bash

        pip install ".[ml]"


Installing Experimental Visualization Libraries
-----------------------------------------------

POSYDON provides experimental visualization libraries to enhance the experience of data analysis and results visualization. While these libraries offer advanced features, please note that they might still be in development and could be subject to changes.

To install these experimental visualization libraries

1. **Navigate to your POSYDON directory** (where the `setup.py` is located) and run:

    .. code-block:: bash
   
        pip install ".[vis]"

    After installing these libraries, you can access various visualization tools and features integrated within POSYDON. Ensure to consult the documentation or any guides associated with these features for their optimal usage.

    .. note::
        As these are experimental features, feedback, and bug reports regarding the visualization tools are highly appreciated. It will aid the development and optimization of these features for future stable releases.
