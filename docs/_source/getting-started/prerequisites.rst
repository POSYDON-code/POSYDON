Prerequisites
-------------

Before proceeding with the installation and usage of POSYDON, there are certain prerequisites that must be met. This page outlines all the necessary requirements for a smooth experience with the tool.

Operating System
~~~~~~~~~~~~~~~~

POSYDON has been tested on the following operating systems:

- Linux (Ubuntu 20.04 and newer, CentOS 7, etc.)
- macOS (10.14 and newer)
<<<<<<< HEAD
- Windows 10
=======
>>>>>>> eirini_CE_fix

(Note: Adapt and expand the list based on the actual compatibility of POSYDON)

Software Dependencies
~~~~~~~~~~~~~~~~~~~~~

1. Python 3.11
<<<<<<< HEAD
2. [List other software dependencies here, if any.]
=======
2. More dependencies are defined in our setup file and should be automatically applied during the installation. To allow you having different versions of the packages we require we strongly encourage to use a dedicated conda environment for POSYDON.
>>>>>>> eirini_CE_fix

(Note: Expand on any specific versions or configurations, if necessary.)

Hardware Requirements
~~~~~~~~~~~~~~~~~~~~~

**Storage**:

<<<<<<< HEAD
In order to run a population synthesis model with POSYDON, you must ensure that you have approximately **40GB** of free storage space. This is crucial for downloading the lite MESA simulation library, interpolation objects, and other auxiliary files used by the code.
=======
In order to run a population synthesis model with POSYDON, you must ensure that you have enough free storage space. It depends on the dataset how much free space you need:

- DR1: approximately **40GB**
- DR2: approximately **140GB** (the default in the :ref:`installation <installation-guide>`)
- one of the eight metallicities in DR2: approximately **20GB**

This is crucial for downloading the lite MESA simulation library, interpolation objects, and other auxiliary files used by the code.
>>>>>>> eirini_CE_fix

Memory:

- Minimum: 8GB RAM
- Recommended: 16GB RAM

(Note: Adjust memory requirements based on your knowledge of POSYDON's demands.)

CPU:

- Multi-core processor recommended for optimal performance, e.g. HPC cluster.

(Note: Include more specific CPU details if necessary.)

Network:

- A stable internet connection for downloading necessary libraries and datasets.

Installation Steps
~~~~~~~~~~~~~~~~~~

Proceed to :ref:`Installation Guide <installation-guide>` to understand the steps required to get POSYDON up and running.

Final Notes
~~~~~~~~~~~

Ensure all the above prerequisites are met before initiating the setup for POSYDON to ensure a seamless and hassle-free experience.

If you encounter any issues, please refer to the :ref:`Troubleshooting Guide <installation-issues>` or :ref:`contact us <contact_info>`.
