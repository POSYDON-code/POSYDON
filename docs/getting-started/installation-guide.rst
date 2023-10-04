.. _installation-guide:

Installation Guide
------------------

This guide will help you step by step in installing POSYDON. We recommend using Anaconda, a package manager, to manage the installation and dependencies, ensuring a smooth setup process.

.. contents:: Table of Contents
   :local:

Using Anaconda (Recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **Install Anaconda**

    If you haven't already, download and install Anaconda from `Anaconda's official website <https://www.anaconda.com/products/distribution>`_.

2. **Create a New Conda Environment**

    Open your terminal or Anaconda prompt and create a new environment using Python 3.11:

    .. code-block:: bash

        conda create --name posydon_env python=3.11 -y

    Activate the new environment:

    .. code-block:: bash

        conda activate posydon_env

3. **Install POSYDON**

    .. warning::
        This documentation describes the POSYDON v2.0.0 code which is not yet available on Anaconda. Please use the development version for now. See :ref:`Using the Development Version <dev-version>` for more details.

    With the environment activated, install POSYDON using:

    .. code-block:: bash

        conda install -c posydon posydon

.. _posydon-env:

4. **Set Environment Variables**

    Export the required paths:

    .. code-block:: bash

        export PATH_TO_POSYDON=/path/to/your/posydon/installation
        export PATH_TO_POSYDON_DATA=/path/where/you/want/to/store/data

    .. note:: 
        You can add these lines to your ``~/.bashrc`` or ``~/.bash_profile`` or your shell equivalent to ensure the environment variables are set every time you open a new terminal.

5. **Download the Dataset**

    .. warning::
        The POSYDON v2.0.0 dataset is not yet available on Zenodo. The above instructions currently point to the POSYDON v1.0.0 dataset release. Please refer to the development version of the dataset available on Northwestern and UNIGE HPC facilities for now. To have access to latest pre-release dataset (230914) you must be a POSYDON core developer, please refer to the #developers Slack channel.

    You can use POSYDON's built-in API command:

    .. code-block:: bash

        get-posydon-data

    Alternatively, you can manually download the dataset from Zenodo using the provided `link <https://zenodo.org/record/6384235>`_. (TODO: update link to v2)

.. _dev-version:

Using the Development Version
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For users interested in the latest features and developments, you can install POSYDON directly from its GitHub repository:

1. **Clone the Repository**

    In your terminal or command prompt:

    .. code-block:: bash

        git clone https://github.com/POSYDON-code/POSYDON.git

2. **Install the Development Version**

    .. warning::
        If you are installing POSYDON on a Mac with Apple M1 or M2 chips, you should first install `hdf5` and `pytables` through conda with `conda install hdf5 pytables`, before following the instractions below.

    Navigate to the cloned repository's directory:

    .. code-block:: bash

        cd POSYDON

    Install the software as an editable package using:

    .. code-block:: bash

        pip install -e .

3. **Set Environment Variables and Download Data**

    Refer back to the recommended installation steps, starting from :ref:`point 4 <posydon-env>`, to set the necessary environment variables and download the required dataset.


Using POSYDON on HPC Facilities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you are planning to run POSYDON's population synthesis on a High-Performance Computing (HPC) facility, it's essential to have `mpi4py` installed to enable parallel computations. 

1. **Install mpi4py via Anaconda (Recommended)**:

    .. code-block:: bash

        conda install mpi4py

2. **Alternatively, via pip**:

    .. code-block:: bash

        pip install ".[hpc]"

    .. warning::
        Users have reported issues when trying to install `mpi4py` via pip. If you encounter any issues, try installing `mpi4py` through Anaconda. If you cannot solve the issue, please refer to the :ref:`Troubleshooting Guide <installation-issues>` or seek support from the community or developers, see the :ref:`contact us <contact-info>` page.

Machine Learning Modules Installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For users who wish to utilize POSYDON's latest machine learning features:

1. **Navigate to your POSYDON directory** (where the `setup.py` is located) and run:

    .. code-block:: bash

        pip install ".[ml]"


Installing Experimental Visualization Libraries
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

POSYDON provides experimental visualization libraries to enhance the experience of data analysis and results visualization. While these libraries offer advanced features, please note that they might still be in development and could be subject to changes.

To install these experimental visualization libraries

1. **Navigate to your POSYDON directory** (where the `setup.py` is located) and run:

    .. code-block:: bash
   
        pip install ".[vis]"

    After installing these libraries, you can access various visualization tools and features integrated within POSYDON. Ensure to consult the documentation or any guides associated with these features for their optimal usage.

    .. note::
        As these are experimental features, feedback, and bug reports regarding the visualization tools are highly appreciated. It will aid the development and optimization of these features for future stable releases.


Documentation Installation & Compilation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you're interested in building the POSYDON documentation locally:

1. **Install Documentation Modules**:

    Navigate to your POSYDON directory and install the required documentation modules:

    .. code-block:: bash

        pip install ".[doc]"

2. **Compile the Documentation**:

    Once you have the required modules installed, you can build the documentation using Sphinx:

    .. code-block:: bash

        cd docs
        make html

3. **Install Pandoc via Anaconda**

    .. warning::
        If you are installing POSYDON on a Mac with Apple M1 or M2 chips, you should install `pandoc` through brew with `brew install pandoc`.

    .. code-block:: bash

        conda install pandoc

4. **Open the Compiled Documentation**:

    After successfully building the documentation, you can view it in your preferred browser. Navigate to the build directory and open the `index.html`:

    .. code-block:: bash

        open _build/html/index.html

    .. note::
        The `open` command works on macOS. If you're using a different OS, you might need to open the `index.html` using your file manager or use a different command.


Installing Jupyter for Tutorials
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The tutorials in the "Getting Started" guide are provided as Jupyter notebooks. If you want to run these notebooks interactively, you will need to have either Jupyter Lab or Jupyter Notebook installed.

1. **Using Anaconda (Recommended)**


    If you have already installed Anaconda as suggested earlier in the installation guide, installing Jupyter Lab or Notebook is straightforward:

    .. code-block:: bash

        conda install -c conda-forge jupyterlab

    Or, for the classic Jupyter Notebook:

    .. code-block:: bash

        conda install -c conda-forge notebook

2. **Alternatively, via pip**


    If you prefer using `pip`, you can also install Jupyter Lab or Notebook using the following commands:

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


Additional Notes
~~~~~~~~~~~~~~~~~

- After installation, ensure you verify the setup by following our :ref:`Verification Guide <verification>`.
- Always ensure you activate the `posydon_env` environment before running POSYDON.
- If you encounter issues during the installation, consult our :ref:`Troubleshooting Guide <installation-issues>` or seek support from the community or developers, see the :ref:`contact us <contact-info>` page.

