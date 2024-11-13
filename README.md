# POSYDON

POSYDON is a next-generation single- and binary-star population synthesis framework, incorporating, fully self-consistent, state-of-the-art stellar structure and evolution modelling, using the [MESA](https://docs.mesastar.org) code, throughout the evolution of both binary components. This allows for a more accurate treatment of physical processes in stellar and binary evolution, including: realistic mass-transfer calculations and assessment of stability, internal angular-momentum transport and tides, stellar core sizes, mass-transfer rates, and orbital periods. The code is modular in many aspects and the user can specify initial population properties and adopt choices that determine how the evolution of a binary proceeds. Machine-learning methods are incorporated and applied on the grids of detailes single- and binary-star models for various classification and interpolation calculations, and the development of irregular grids guided by active learning, for computational efficiency.

POSYDON is being developed by a collaborative team of astrophysicists and computer scientists led by Principal Investigators Tassos Fragos (Université de Genève) and Vicky Kalogera (Northwestern University).  

In [Fragos et al. (2023)](https://ui.adsabs.harvard.edu/abs/2023ApJS..264...45F/abstract) and [Andrews et al. (2024)](https://ui.adsabs.harvard.edu/abs/2024arXiv241102376A/abstract),  we describe the detailed methodology and implementation of POSYDON, including the assumed physics of stellar and binary evolution, the extensive grids of detailed single- and binary-star models, the postprocessing, classification, and interpolation methods we developed for use with the grids, and the treatment of evolutionary phases that are not based on precalculated grids.


### For instruction on how to install and use POSYDON [see here](https://posydon.org/docs/index.html).
