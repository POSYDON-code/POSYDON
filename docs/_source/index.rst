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

How is this Documentation Structured?
-------------------------------------

- **Introductions** Discover what POSYDON is, its objectives, and its science scope.
- **Getting Started:** A quick guide to get POSYDON up and running on your machine.
- **User Guides:** Detailed guides on using various features.
- **Tutorials and Examples:** Step-by-step walk-throughs of typical use-cases.
- **In-Depth Components Overview:** Delve deep into the core components of POSYDON.
- **API Reference:** A comprehensive reference for developers.
- **Troubleshooting and FAQs:** Answers to common questions and problems.
- **Contributing to POSYDON:** Join our community and help us grow!
- **Release Notes and Changelog:** Stay updated with our version history.
- **Contact and Support:** Reach out to us with your queries or feedback.

.. - **POSYDON Workflow:** Learn about how to conduct simulations from start to finish.

.. toctree::
   :maxdepth: 1
   :caption: Introduction

   introduction-acknowledgements/intro
   Collaborative Team <https://posydon.org/team.html>
   Publications <https://ui.adsabs.harvard.edu/public-libraries/ZZsD9bzLTzWnLV3hwyJxbA>

.. toctree::
   :maxdepth: 1
   :caption: Getting Started

   getting-started/prerequisites
   getting-started/installation-guide

.. toctree::
   :maxdepth: 1
   :titlesonly:
   :caption: Tutorials and Examples

   tutorials-examples/population-synthesis/binary-pop-syn
   tutorials-examples/generating-datasets/generating-datasets
   tutorials-examples/MESA-grids/running-grids

.. toctree::
   :maxdepth: 1
   :titlesonly: 
   :caption: In-Depth Components Overview

   components-overview/mesa-grids
   components-overview/processing-pipeline
   components-overview/machine-learning-components
   components-overview/stellar-binary-simulation

.. .. toctree::
..    :maxdepth: 1
..    :caption: POSYDON Workflow

..    posydon-workflow/population-specification
..    posydon-workflow/stellar-evolution-choices
..    posydon-workflow/mesa-evolutionary-tracks
..    posydon-workflow/interacting-binary-grids
..    posydon-workflow/running-simulations

.. toctree::
   :maxdepth: 1
   :caption: API Reference

   api_reference/posydon
   api_reference/bin

.. toctree::
   :maxdepth: 1
   :caption: Troubleshooting and FAQs

   troubleshooting-faqs/installation-issues
   troubleshooting-faqs/code-questions


.. toctree::
   :maxdepth: 1
   :caption: Contributing to POSYDON

   contributing/how-to-contribute
   contributing/code-style-guidelines
   contributing/reporting-issues


.. toctree::
   :maxdepth: 1
   :caption: Releases and Datasets

   Releases <https://github.com/POSYDON-code/POSYDON/releases>
   Datasets <https://zenodo.org/communities/posydon/>
   

.. toctree::
   :maxdepth: 1
   :caption: Contact and Support

   contact-support/contact-information





Acknowledgments
---------------

**Team Members:** POSYDON is being led and developed by a dedicated team of astrophysicists and computer scientists. At the helm are Principal Investigators Tassos Fragos (Université de Genève) and Vicky Kalogera (Northwestern University), along with many talented individuals. You can read more about the team on the :ref:`Collaborative Team <https://posydon.org/team.html>`_ page.

**Funding Agencies:** The POSYDON project is supported primarily by two sources: the Swiss National Science Foundation Professorship grant (PI Fragos) and the Gordon and Betty Moore Foundation (PI Kalogera).

**Licensing:** POSYDON is licensed under "BSD 3-Clause". For full licensing details, please refer to our `LICENSE <https://github.com/POSYDON-code/POSYDON/blob/main/LICENSE.md>`_.
