.. _installation-issues:

Common Installation Issues
-------------------

From time to time, users might encounter issues during the installation of POSYDON. This page aims to address common installation problems and offer solutions. If your problem isn't covered here, please `report the issue <https://github.com/POSYDON-code/POSYDON/issues>`_ so we can assist you and possibly update this page for the benefit of others.

1. **Slow `conda` solving or installation taking hours:**

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

    c. **Use mamba** (Alternative): Another option is to use the `mamba` package manager, which is a drop-in replacement for `conda` but is much faster (`click here for more details <https://www.anaconda.com/blog/a-faster-conda-for-a-growing-community>`_). However, this approach has not been fully vetted with POSYDON.

    .. note::
        If you're on an HPC cluster with an old system-wide conda installation (e.g., conda 2021.11), you may need to install a recent conda version locally in your home directory rather than using the system version.

    `conda` can also be slow and sometimes gets stuck on "Verifying transaction" or "Executing transaction", especially when installing packages on a cluster, as it creates many small files which can be difficult for HPC clusters to handle. The libmamba solver helps with this issue as well.

2. **Failed Dependencies**:
    - **Description**: Sometimes, certain dependencies might fail to install or conflict with pre-existing ones.
    - **Solution**: Try installing the failed dependencies separately using ``pip`` or ``conda`` before installing POSYDON. If using ``conda``, consider creating a fresh environment specifically for POSYDON. Additionally, some of the package versions are only available on ``conda-forge``. Please ensure that ``conda-forge`` has been added to conda as a channel using the command ``conda config --add channels conda-forge``.

    We try our best to keep the dependencies updated and compatible. However, if you encounter issues, please let us know!

3. **Failed to Set Environment Variables**:
    - **Description**: After installation, you might forget to set the ``PATH_TO_POSYDON`` and ``PATH_TO_POSYDON_DATA`` environment variables.
    - **Solution**: Refer back to the :ref:`installation guide <installation-guide>` to ensure all post-installation steps are followed.

4. **Insufficient Storage Space**:
    - **Description**: POSYDON requires around 140GB of free storage space for the down-sampled MESA simulation library and related files.
    - **Solution**: Ensure you have enough space on your installation drive. Delete unnecessary files or consider using a larger storage solution.

Troubleshooting Tips
~~~~~~~~~~~~~~~~~~~~

1. **Check Error Messages**: Always read the error messages during the installation. They often provide hints or exact reasons for the failure.

2. **Use a Virtual Environment**: Installing POSYDON in a virtual environment (e.g., via ``conda`` or ``venv``) can prevent conflicts with system-wide packages.

3. **Update Your Package Managers**: Ensure that both ``pip`` and ``conda`` (if you use them) are updated to their latest versions.

4. **System Compatibility**: Double-check the compatibility of POSYDON with your operating system and Python version.

5. **Check Online Forums**: Platforms like Stack Overflow often have threads related to common Python package installation issues.

Still Facing Issues?
~~~~~~~~~~~~~~~~~~~~

If you've tried the solutions above and are still encountering problems, please :ref:`contact us <contact_info>` with as much detail as possible. We're here to help!

Thank you for your patience, and we appreciate your interest in POSYDON!
