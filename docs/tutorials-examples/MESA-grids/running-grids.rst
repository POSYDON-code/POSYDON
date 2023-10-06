.. _MESA-grids:

Architecting MESA Simulation Grids with POSYDON
===============================================

Boost your MESA simulation grid prowess by diving into the intricacies of POSYDON's grid creation techniques and the dedicated MESA configuration repository. Harness the power of our API to craft precise and customized grids.

Getting Started Tutorials
-------------------------

The POSYDON MESA simulation parameters are stored in `POSYDON-MESA-INLISTS <https://github.com/POSYDON-code/POSYDON-MESA-INLISTS>`_ GitHub submodule. To obtain the submodule, clone the POSYDON repository with the following command in your terminal:

.. code-block:: bash

    cd $PATH_TO_POSYDON
    git submodule init
    git submodule update grid_params/POSYDON-MESA-INLISTS/

Similarly to the POSYDON code repository, there is a `main` branch that points to the latest stable version of the POSYDON-MESA-INLISTS repository associated to a POSYDON code release and a `development` branch that contains the latest stable version of our fiducial MESA configurations, run the following command to get the lates development version of the sumodule:

.. code-block:: bash

    cd grid_params/POSYDON-MESA-INLISTS/
    git checkout development

To follow the next step of the tutorial, you will need to have MESA installed on your machine. If you do not have MESA installed, please follow the instructions on the `MESA website <https://docs.mesastar.org/en/release-r23.05.1/>`_.

.. warning::

    The POSYDON v2.0.0 code is compatible with MESA r11701. Support might not be available for the latest MacOS version.


MESA Configuration Repository
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Immerse yourself in our Jupyter Notebook which guides you through our dedicated repository, showcasing the vast array of MESA configurations available.

.. toctree::
    
    notebook: Exploring the MESA Configuration Repository

Harnessing the API Interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Step-by-step guide on how to leverage our API for creating and customizing MESA grids to suit your needs.

.. toctree::
    
    notebook: Mastering Grid Creation with the API

Advanced Tutorials
------------------

Optimizing Grid Creation
~~~~~~~~~~~~~~~~~~~~~~~~

Dive deeper into advanced techniques for grid creation, ensuring efficiency and precision.
.. toctree::
    
    notebook: Advanced Grid Crafting Techniques

Automation with the API
~~~~~~~~~~~~~~~~~~~~~~~

Explore automation possibilities with the API, streamlining your grid creation processes.

.. toctree::
   
   notebook: API-driven Grid Automation

Deep Dive: Behind the Scene Customizations
------------------------------------------

Custom Grid Configurations
~~~~~~~~~~~~~~~~~~~~~~~~~~

Learn the art of crafting specialized grid configurations tailored for specific simulations.

.. toctree::

   notebook: Tailoring Custom Grid Configurations

Advanced API Features
~~~~~~~~~~~~~~~~~~~~~

Delve into the more intricate functionalities of the API, unlocking extensive customization options for MESA grid creation.

.. toctree::

    notebook: Diving Deeper into API Features

Grid Integrity and Validation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Ensure the grids you create are of the highest standard by performing dedicated validation checks.

.. toctree::

    notebook: Grid Integrity and Assurance

Support & Feedback
------------------

Encountered obstacles or have questions? Ensure you consult our [FAQ](link-to-FAQ) or drop by the [contact-information.rst](contact-information.rst) page.

