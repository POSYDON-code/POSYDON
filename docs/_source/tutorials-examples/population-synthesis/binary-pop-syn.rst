.. _binary-pop-syn:


Binary-Star Population Synthesis with POSYDON
==============================================

Embark on a journey into the fascinating world of binary-star population synthesis with POSYDON! 
Whether you are a budding astrophysicist or a seasoned researcher, our tutorials guide you every step of the way.
From simulating your first binaries to running large-scale populations and diving deep into data analysis, POSYDON offers you a seamless experience.

Getting Started Tutorials
-------------------------
    .. warning::
        The tutorials are designed to be followed in order, as some depend on the outputs of the previous tutorials.

These tutorials are created to get started with POSYDON on your local machine and on a HPC facility.

I. Your First Binary Simulations with POSYDON 🌠
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Simulate your first 10 binaries and experience the power of POSYDON.

.. toctree::
    :maxdepth: 1

    10_binaries_pop_syn

Now that you know how to run your first simulations, explore the different customization options reading the :ref:`population prameters <pop-params-guide>` documentation.


II. Large-Scale Population Synthesis on HPC Facilities 🚀
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this tutorial, you'll learn to run massive simulations of 1 million binaries across 8 different metallicities.

.. note:: 
    
    Ensure you've installed `mpi4py` to use this tutorial. If not, follow our :ref:`installation guide <installation-guide>`.

.. toctree::
    :maxdepth: 1

    pop_syn

III. Analyzing Merging DCO Populations: Rates & Observations 🔍
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Delve into the results of your simulations learning how to analyze the results of your simulations and compute rates and observational properties of merging double compact object (DCO) populations.

.. toctree::
    :maxdepth: 1

    bbh_analysis



Advanced Tutorials
------------------

Computing Long-Duration Gamma-Ray Bursts from Double Compact Object Populations 🌌
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this tutorial, you will learn how to compute long-duration gamma-ray bursts rates associated to the formation of merging binary black holes

.. toctree::
    :maxdepth: 1

    lgrb_pop_syn


Explore the assumption of star-formation history and compute DCO rates for one single metallicity 📖
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this tutorial, you will learn how to change the assumption of star-formation history and how to analysise and compute rates for a population synthesis model consisting of one single metallicity.

.. toctree::
    :maxdepth: 1

    one_met_pop_syn

X-ray binaries: computing the X-ray luminosity function 🩻 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TODO: bring v1 tutorial to v2 leveraging the SyntheticPopulation class

Advanced Visualization with Van den Heuvel Diagrams 🎨
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Dive deeper into POSYDON's visualization capabilities. This tutorial explores the intricacies of the Van den Heuvel diagrams, allowing users to intuitively visualize and analyze individual or multiple binaries from the population synthesis models. Discover multiple visualization modes, interactive features, and binary analysis tools to enhance your understanding of the synthesized data. Perfect for those looking to gain a richer visual perspective on their simulations.

.. toctree::
    :maxdepth: 1

    Van den Heuvel Diagrams <VHD>


Deep Dive: Behind the Scenes Customizations
-------------------------------------------

Dive into the nitty-gritty of POSYDON's inner workings. Customize, extend, and tweak the system to fit your unique requirements:

Customizing the Population Synthesis Class 🛠️
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   - Understand the foundational class driving our binary-star population synthesis. Learn to modify and adapt it to your needs.
   - 🔗 [Delve into class customizations](link-to-customization-notebook1)

(Add other in-depth customization tutorials here...)

Remember, mastering POSYDON is a journey, not a destination. As you progress through these tutorials, you'll uncover layers of capabilities and functionalities. Keep exploring and happy simulating!
