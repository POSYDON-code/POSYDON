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

        conda create --name posydon_env python=3.11

    Activate the new environment:

    .. code-block:: bash

        conda activate posydon_env

3. **Install POSYDON**

    With the environment activated, install POSYDON using:

    .. code-block:: bash

        conda install -c posydon posydon

4. **Set Environment Variables**

    Export the required paths:

    .. code-block:: bash

        export PATH_TO_POSYDON=/path/to/your/posydon/installation
        export PATH_TO_POSYDON_DATA=/path/where/you/want/to/store/data

5. **Download the Dataset**

    You can use POSYDON's built-in API command:

    .. code-block:: bash

        get-posydon-data

    Alternatively, you can manually download the dataset from Zenodo using the provided `link <https://zenodo.org/record/6384235>`_. (TODO: update link to v2)

Using the Development Version
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For users interested in the latest features and developments, you can install POSYDON directly from its GitHub repository:

1. **Clone the Repository**

    In your terminal or command prompt:

    .. code-block:: bash

        git clone https://github.com/POSYDON-code/POSYDON.git

2. **Install the Development Version**

    Navigate to the cloned repository's directory:

    .. code-block:: bash

        cd POSYDON

    Install the software using:

    .. code-block:: bash

        pip install .

3. **Set Environment Variables and Download Data**

    Refer back to the recommended installation steps, starting from point 4, to set the necessary environment variables and download the required dataset.

Using POSYDON on HPC Facilities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you are planning to run POSYDON's population synthesis on a High-Performance Computing (HPC) facility, it's essential to have `mpi4py` installed to enable parallel computations. 

1. **Install mpi4py via Anaconda**:

    .. code-block:: bash

        conda install -c anaconda mpi4py

2. **Alternatively, via pip**:

    .. code-block:: bash

        pip install mpi4py

Machine Learning Modules Installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For users who wish to utilize POSYDON's latest machine learning features:

1. **Navigate to your POSYDON directory** (where the `setup.py` is located) and run:

    .. code-block:: bash

        pip install .[ml]


Installing Experimental Visualization Libraries
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
POSYDON provides experimental visualization libraries to enhance the experience of data analysis and results visualization. While these libraries offer advanced features, please note that they might still be in development and could be subject to changes.

To install these experimental visualization libraries, you can use the following pip command:

.. code-block:: bash

   pip install .[vis]

After installing these libraries, you can access various visualization tools and features integrated within POSYDON. Ensure to consult the documentation or any guides associated with these features for their optimal usage.

.. note::
   As these are experimental features, feedback, and bug reports regarding the visualization tools are highly appreciated. It will aid the development and optimization of these features for future stable releases.


Documentation Installation & Compilation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you're interested in building the POSYDON documentation locally:

1. **Install Documentation Modules**:

    Navigate to your POSYDON directory and install the required documentation modules:

    .. code-block:: bash

        pip install .[doc]

2. **Compile the Documentation**:

    Once you have the required modules installed, you can build the documentation using Sphinx:

    .. code-block:: bash

        cd docs
        make html

3. **Open the Compiled Documentation**:

    After successfully building the documentation, you can view it in your preferred browser. Navigate to the build directory and open the `index.html`:

    .. code-block:: bash

        open _build/html/index.html

    Note: The `open` command works on macOS. If you're using a different OS, you might need to open the `index.html` using your file manager or use a different command.


Additional Notes
~~~~~~~~~~~~~~~~~

- After installation, ensure you verify the setup by following our :ref:`Verification Guide <verification>`.
- Always ensure you activate the `posydon_env` environment before running POSYDON.
- If you encounter issues during the installation, consult our :ref:`Troubleshooting Guide <installation-issues>` or seek support from the community or developers, see the :ref:`contact us <contact-info>` page.

