.. _installation-issues:

Installation Issues
-------------------

From time to time, users might encounter issues during the installation of POSYDON. This page aims to address common installation problems and offer solutions. If your problem isn't covered here, please `report the issue <https://github.com/POSYDON-code/POSYDON/issues>`_ so we can assist you and possibly update this page for the benefit of others.

Common Installation Issues
~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **Failed Dependencies**:
    - **Description**: Sometimes, certain dependencies might fail to install or conflict with pre-existing ones.
    - **Solution**: Try installing the failed dependencies separately using ``pip`` or ``conda`` before installing POSYDON. If using ``conda``, consider creating a fresh environment specifically for POSYDON.

2. **Insufficient Storage Space**:
    - **Description**: POSYDON requires around 40GB of free storage space for the lite MESA simulation library and related files.
    - **Solution**: Ensure you have enough space on your installation drive. Delete unnecessary files or consider using a larger storage solution.

3. **Proxy or Network Issues**:
    - **Description**: Installation might fail if you're behind a strict network proxy or firewall.
    - **Solution**: If you're using ``conda``, you can set proxy settings using the ``--proxy`` flag. For pip, you can use the ``--proxy`` flag as well.

4. **Failed to Set Environment Variables**:
    - **Description**: After installation, you might forget to set the ``PATH_TO_POSYDON`` and ``PATH_TO_POSYDON_DATA`` environment variables.
    - **Solution**: Refer back to the [installation guide](link-to-installation-guide) to ensure all post-installation steps are followed.

5. **Error with Experimental Visualization Libraries**:
    - **Description**: Issues after installing the experimental visualization libraries with ``pip install ".[vis]"``.
    - **Solution**: These libraries are experimental. Please ensure you have all required system dependencies. Consider reaching out to our support channels for more guidance or troubleshooting.

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
