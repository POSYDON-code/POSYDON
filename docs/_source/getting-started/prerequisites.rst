Prerequisites
-------------

Before proceeding with the installation and usage of POSYDON, there are certain prerequisites that must be met. This page outlines all the necessary requirements for a smooth experience with the tool.

Operating System
~~~~~~~~~~~~~~~~

POSYDON has been tested on the following operating systems:

- Linux (Ubuntu 20.04 and newer, CentOS 7, etc.)
- macOS (10.14 and newer)

POSYDON is developed and tested on Linux (and in particular on HPC clusters),
which is the platform POSYDON is officially tested against. We aim to keep
macOS working, but if you hit an issue on macOS please see the
:ref:`troubleshooting guide <installation-issues>`.

Software Dependencies
~~~~~~~~~~~~~~~~~~~~~

1. Python 3.11
2. More dependencies are defined in our setup file (``pyproject.toml``) and are
   installed automatically during the installation. To keep different package
   versions from interfering with each other, we strongly encourage you to use a
   dedicated conda environment for POSYDON.

Hardware Requirements
~~~~~~~~~~~~~~~~~~~~~

**Storage**: The exact footprint depends on which data release and how many
metallicities you download. See :ref:`data-releases` for the canonical numbers.

In order to run a population synthesis model with POSYDON, you must ensure that you have enough free storage space. It depends on the dataset how much free space you need:

- DR2: approximately **140GB** (the default in the :ref:`installation <installation-guide>`)
- one of the eight metallicities in DR2: approximately **20GB**
- DR1: approximately **40GB**

This is crucial for downloading the lite MESA simulation library, interpolation objects, and other auxiliary files used by the code.

Memory:

- Minimum: 8GB RAM
- Recommended: 16GB RAM

The amount of memory a population run needs per CPU, and the related
``dump_rate`` setting, are covered in the :ref:`memory tuning <memory-tuning>`
how-to.

CPU:

- A multi-core processor is recommended for optimal performance (e.g. an HPC cluster), since a run is parallelized across CPUs (one process per CPU).

Installation Steps
~~~~~~~~~~~~~~~~~~

Proceed to :ref:`Installation Guide <installation-guide>` to understand the steps required to get POSYDON up and running.

Final Notes
~~~~~~~~~~~

Ensure all the above prerequisites are met before initiating the setup for POSYDON to ensure a seamless and hassle-free experience.

If you encounter any issues, please refer to the :ref:`Troubleshooting Guide <installation-issues>` or :ref:`contact us <contact_info>`.
