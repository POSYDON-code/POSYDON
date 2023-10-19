.. _intro:

Introduction, Objectives, and Scope
-----------------------------------

What is POSYDON?
~~~~~~~~~~~~~~~~

POSYDON, an acronym for **POpulation SYnthesis with Detailed binary-evolution simulatiONs**, is a next-generation single and binary-star population synthesis code. Designed to seamlessly incorporate full stellar structure and evolution modeling, it leverages the advanced capabilities of MESA for a comprehensive and accurate simulation of stellar systems.

Code Objectives
~~~~~~~~~~~~~~~~

The primary goal of POSYDON is to simulate stellar populations with unparalleled precision and detail. The modular nature of the code allows users to specify initial population properties and make informed choices regarding stellar evolution paths. By using MESA evolutionary tracks, POSYDON provides simulations for single, non-interacting, and interacting binaries organized in grids. With the integration of advanced machine-learning methods, the code performs classifications and a variety of interpolation calculations. Active learning guides the development of irregular grids, ensuring computational efficiency and robustness.

Science Scope
~~~~~~~~~~~~~

With its sophisticated algorithms and detailed simulations, POSYDON has broad science applicability. The version 2 of POSYDON is tailored towards the evolution of massive binary stars. This includes their intricate paths leading to the formation of neutron stars and black holes. Furthermore, the code is equipped to handle and analyze these evolutions across a vast range of metallicities, accommodating various stellar environments and conditions.

Citation Disclaimer
~~~~~~~~~~~~~~~~~~~

If you use POSYDON for your research and it contributes to a publication, we kindly request that you cite our foundational papers for this version:

- `POSYDON: A General-purpose Population Synthesis Code with Detailed Binary-evolution Simulations <https://ui.adsabs.harvard.edu/abs/2023ApJS..264...45F/abstract>`_, **Fragos et al. (2023)**

.. code-block::

        @ARTICLE{2023ApJS..264...45F,
            author = {{Fragos}, Tassos and {Andrews}, Jeff J. and {Bavera}, Simone S. and {Berry}, Christopher P.~L. and {Coughlin}, Scott and {Dotter}, Aaron and {Giri}, Prabin and {Kalogera}, Vicky and {Katsaggelos}, Aggelos and {Kovlakas}, Konstantinos and {Lalvani}, Shamal and {Misra}, Devina and {Srivastava}, Philipp M. and {Qin}, Ying and {Rocha}, Kyle A. and {Rom{\'a}n-Garza}, Jaime and {Serra}, Juan Gabriel and {Stahle}, Petter and {Sun}, Meng and {Teng}, Xu and {Trajcevski}, Goce and {Tran}, Nam Hai and {Xing}, Zepei and {Zapartas}, Emmanouil and {Zevin}, Michael},
            title = "{POSYDON: A General-purpose Population Synthesis Code with Detailed Binary-evolution Simulations}",
            journal = {\apjs},
            keywords = {Binary stars, Close binary stars, Compact binary stars, Interacting binary stars, X-ray binary stars, Compact objects, Stellar remnants, Black holes, Neutron stars, Gravitational wave sources, Stellar evolutionary models, Stellar populations, 154, 254, 283, 801, 1811, 288, 1627, 162, 1108, 677, 2046, 1622, Astrophysics - Solar and Stellar Astrophysics},
            year = 2023,
            month = feb,
            volume = {264},
            number = {2},
            eid = {45},
            pages = {45},
            doi = {10.3847/1538-4365/ac90c1},
            archivePrefix = {arXiv},
            eprint = {2202.05892},
            primaryClass = {astro-ph.SR},
            adsurl = {https://ui.adsabs.harvard.edu/abs/2023ApJS..264...45F},
            adsnote = {Provided by the SAO/NASA Astrophysics Data System}
        }


- POSYDON Version 2:  Population Synthesis across Metallicities with Detailed Binary-Evolution Simulation, Andrews et al. (in prep.)

Acknowledging these works ensures that the development team receives proper credit and supports the continued advancement and maintenance of the POSYDON software.
