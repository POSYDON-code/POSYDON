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

Using these Documentation
-------------------------

These documentation are organized according to the `Diátaxis <https://diataxis.fr/>`_ framework, which groups material into four distinct modes of use. As you work with POSYDON, you will move between them naturally:

- **Learn (Tutorials):** Step-by-step, goal-oriented lessons that let you experience POSYDON by building confidence through practice. Start here if you are new to the code and want to get a first simulation running.
- **Do (How-to):** Task-focused guides that solve a specific problem or accomplish a concrete objective (running on an HPC cluster, tuning memory, adding your own step, contributing). Use these when you know what you want to achieve but need the exact recipe.
- **Understand (Explanation):** Background and conceptual material that explains *why* POSYDON works the way it does. These pages help you form a mental model of grids, the post-processing pipeline, and population synthesis.
- **Look up (Reference):** Technical descriptions of the code and configuration: the API, command-line executables, and parameter files. Use these when you need an authoritative statement of what something is.

The navigation below groups every page under one of these four headings, so each page has an unambiguous home depending on what you are trying to do.

.. toctree::
   :maxdepth: 1
   :caption: Learn (Tutorials)

   tutorials-examples/population-synthesis/binary-pop-syn
   tutorials-examples/generating-datasets/generating-datasets
   tutorials-examples/MESA-grids/running-grids

.. toctree::
   :maxdepth: 1
   :caption: Do (How-to)

   getting-started/prerequisites
   getting-started/installation-guide
   getting-started/first-grids
   getting-started/first-population
   howto/run-population-hpc
   howto/analyze-population
   howto/memory-tuning
   howto/fix-reruns
   howto/add-step-user-module
   howto/contribute

.. toctree::
   :maxdepth: 1
   :caption: Understand (Explanation)

   explanation/workflow
   explanation/grid-types
   explanation/data-releases
   explanation/glossary
   components-overview/stellar-binary-simulation
   components-overview/processing-pipeline
   components-overview/machine-learning-components

.. toctree::
   :maxdepth: 1
   :caption: Look up (Reference)

   components-overview/mesa-grids
   reference/config
   api_reference/posydon
   api_reference/bin

.. toctree::
   :maxdepth: 1
   :titlesonly:
   :caption: POSYDON school material

   POSYDON School 2025 <https://github.com/POSYDON-code/POSYDON-2025-School-Labs>

.. toctree::
   :maxdepth: 1
   :caption: Releases and Datasets

   Releases <https://github.com/POSYDON-code/POSYDON/releases>
   Datasets <https://zenodo.org/communities/posydon/>

.. toctree::
   :maxdepth: 1
   :caption: Support and Contact

   introduction-acknowledgements/intro
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
