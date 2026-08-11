==============================
POSYDON Documentation
==============================


.. image:: posydon_logo.*
  :width: 800



Welcome to POSYDON!
-------------------

Welcome to the official documentation of **POSYDON** - *POpulation SYnthesis with Detailed binary-evolution simulatiONs*. Whether you're a new user, a contributor, or just exploring, we're delighted to have you here!

About POSYDON
-------------

POSYDON is a next-generation single- and binary-star population synthesis tool. Our vision is to provide researchers and enthusiasts with a state-of-the-art platform to delve into the intricacies of stellar structures and binary evolutions using MESA. With full stellar structure modeling, advanced machine learning techniques, and a modular architecture, POSYDON stands at the forefront of astrophysical simulations.

To stay up-to-date with the latest news about POSYDON, check out our `official website <https://www.posydon.org>`_ for more details.


Using this Documentation
------------------------

We've tried to group the documentation into four distinct sections to help you navigate easily:

- **Tutorials:** Step-by-step guides to help you get started and explore POSYDON's capabilities.
- **User Guides:** In-depth explanations of POSYDON's features and functionalities.
- **Reference:** For developers and advanced users, detailing the classes, methods, and functions available in POSYDON.
- **Explanation:** Background and conceptual material on the POSYDON framework.


.. toctree::
   :maxdepth: 1
   :caption: How-to (Getting Started)

   getting-started/prerequisites
   getting-started/installation-guide
   getting-started/first-grids
   getting-started/first-population

.. toctree::
   :maxdepth: 1
   :caption: Tutorials

   tutorials-examples/population-synthesis/binary-pop-syn
   tutorials-examples/generating-datasets/generating-datasets
   tutorials-examples/MESA-grids/running-grids

.. toctree::
   :maxdepth: 1
   :titlesonly:
   :caption: POSYDON school material

   POSYDON School 2025 <https://github.com/POSYDON-code/POSYDON-2025-School-Labs>



.. toctree::
   :maxdepth: 1
   :titlesonly:
   :caption: Explanation

   explanation/grid-types
   explanation/data-releases
   components-overview/processing-pipeline
   components-overview/machine-learning-components
   components-overview/stellar-binary-simulation

.. toctree::
   :maxdepth: 1
   :caption: Reference

   components-overview/mesa-grids
   api_reference/posydon
   api_reference/bin


.. toctree::
   :maxdepth: 1
   :caption: Introduction

   introduction-acknowledgements/intro
   Collaborative Team <https://posydon.org/team.html>
   Publications <https://ui.adsabs.harvard.edu/public-libraries/ZZsD9bzLTzWnLV3hwyJxbA>


.. toctree::
   :maxdepth: 1
   :caption: Releases and Datasets

   Releases <https://github.com/POSYDON-code/POSYDON/releases>
   Datasets <https://zenodo.org/communities/posydon/>


.. toctree::
   :maxdepth: 1
   :caption: Support and Contact

   contact-support/contact-information
   troubleshooting-faqs/installation-issues
   troubleshooting-faqs/code-questions
   Report an Issue <https://github.com/POSYDON-code/POSYDON/issues/new/choose>
   contributing/how-to-contribute

Acknowledgments
---------------

**Team Members:** POSYDON is being led and developed by a dedicated team of astrophysicists and computer scientists. At the helm are Principal Investigators Tassos Fragos (Université de Genève) and Vicky Kalogera (Northwestern University), along with many talented individuals. You can read more about the team on the `Collaborative Team <https://posydon.org/team.html>`_ page.

**Funding Agencies:** The POSYDON project is supported primarily by two sources: the Swiss National Science Foundation Professorship grant (PI Fragos) and the Gordon and Betty Moore Foundation (PI Kalogera).

**Licensing:** POSYDON is licensed under "BSD 3-Clause". For full licensing details, please refer to our `LICENSE <https://github.com/POSYDON-code/POSYDON/blob/main/LICENSE.md>`_.
