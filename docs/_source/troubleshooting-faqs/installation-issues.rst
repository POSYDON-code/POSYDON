.. _installation-issues:

Common Installation Issues
--------------------------

From time to time, users might encounter issues during the installation of POSYDON. This page aims to address common installation problems and offer solutions. If your problem isn't covered here, please `report the issue <https://github.com/POSYDON-code/POSYDON/issues>`_ so we can assist you and possibly update this page for the benefit of others.

1. **Slow ``conda`` solving or installation taking hours:**

    **Conda Version Requirements**: POSYDON has a complex dependency tree, and older conda versions (especially those prior to v23.1.0) use very slow dependency solvers that can take hours or may never complete the installation process. Modern conda versions (>= 23.1.0) include the libmamba solver by default, which resolves dependencies efficiently in seconds to minutes.

    **Check your conda version and solver:**

    .. code-block:: bash

        conda --version
        conda config --show solver

    **Solutions:**

    a. **Update conda** (Recommended): If you have administrative access or can install conda locally, we strongly recommend updating to the latest conda version (2025 or later), which includes the fast libmamba solver by default:

       .. code-block:: bash

           # Update conda in your base environment
           conda update -n base conda

           # Or install a fresh conda distribution from https://www.anaconda.com/download

    b. **Install and configure the libmamba solver** (For existing conda installations):

       If you cannot update conda but have version 4.12 or later, you can install and configure the libmamba solver:

       .. code-block:: bash

           # Install the libmamba solver
           conda install -n base conda-libmamba-solver

           # Set libmamba as the default solver
           conda config --set solver libmamba

           # Verify the configuration
           conda config --show solver

    c. **Use mamba** (Alternative): Another option is to use the ``mamba`` package manager, which is a drop-in replacement for ``conda`` but is much faster (`click here for more details <https://www.anaconda.com/blog/a-faster-conda-for-a-growing-community>`_). However, this approach has not been fully vetted with POSYDON.

    .. note::
        If you're on an HPC cluster with an old system-wide conda installation (e.g., conda 2021.11), you may need to install a recent conda version locally in your home directory rather than using the system version.

    ``conda`` can also be slow and sometimes gets stuck on "Verifying transaction" or "Executing transaction", especially when installing packages on a cluster, as it creates many small files which can be difficult for HPC clusters to handle. The libmamba solver helps with this issue as well.

2. **Failed Dependencies**:
    - **Description**: Sometimes, certain dependencies might fail to install or conflict with pre-existing ones.
    - **Solution**: Install POSYDON in a fresh conda environment. Some of the required package versions are only available on ``conda-forge``, so make sure it is added as a channel with ``conda config --add channels conda-forge``. If a single dependency still fails, try installing it separately with ``conda`` before installing POSYDON.

    We try our best to keep the dependencies updated and compatible. However, if you encounter issues, please let us know!

3. **Environment Variables Not Set**:
    - **Description**: POSYDON cannot find its data, or complains that ``PATH_TO_POSYDON`` or ``PATH_TO_POSYDON_DATA`` is not set.
    - **Solution**: Check the variables with ``echo $PATH_TO_POSYDON $PATH_TO_POSYDON_DATA``. If they are empty, set them as described in the :ref:`installation guide <installation-guide>`. Variables set with ``conda env config vars set`` only take effect after re-activating the environment.

4. **Insufficient Storage Space**:
    - **Description**: The data download fails or fills up your disk.
    - **Solution**: Download only the metallicities you need (about 12 GB each) instead of the full DR2 dataset (about 110 GB), see :ref:`download-data`, or point ``PATH_TO_POSYDON_DATA`` to a larger disk.

.. _installation-issues-macos:

5. **Installation Fails on Apple Silicon Macs**:
    - **Description**: On Macs with Apple M-series chips, ``pip install -e .`` can fail while building ``tables`` (PyTables) or ``h5py``.
    - **Solution**: Install the HDF5 libraries from conda first, then repeat the ``pip`` installation:

      .. code-block:: bash

          conda install hdf5 pytables

6. **mpi4py Fails to Install**:
    - **Description**: ``pip install ".[hpc]"`` or ``pip install mpi4py`` fails to compile.
    - **Solution**: Install it with ``conda install mpi4py`` instead. On HPC clusters, you may need to load the cluster's MPI module first; ask your cluster's support team which MPI implementation to use.

Still Facing Issues?
~~~~~~~~~~~~~~~~~~~~

If you've tried the solutions above and are still encountering problems, please `open an issue <https://github.com/POSYDON-code/POSYDON/issues>`_ or :ref:`contact us <contact_info>`. Include your operating system, the output of ``conda --version``, the commands you ran and the full error message.
