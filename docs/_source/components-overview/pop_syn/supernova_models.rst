.. _supernova_models:

Supernova Models
================

POSYDON supports different supernova models for the formation of compact objects.
The supernova model can be set in the configuration file of the population synthesis model.

Overview
--------

The supernova treatment in POSYDON is determined by two key aspects:

1. **Core-collapse mechanism**: The physical prescription for computing remnant properties
2. **Computational approach**: How and when the collapse outcome is calculated

Core-Collapse Mechanisms
------------------------

POSYDON supports the following core-collapse mechanisms:

.. list-table:: Core-collapse mechanisms
   :header-rows: 1
   :widths: 20 10 70

   * - Mechanism
     - Engine
     - Description
   * - **Fryer+12-delayed**
     - ``None``
     - Classical rapid and delayed supernova engines from Fryer et al. (2012)
   * - **Fryer+12-rapid**
     - ``None``
     - Classical rapid and delayed supernova engines from Fryer et al. (2012)
   * - **Sukhbold+16-engine**
     - ``N20``, ``W20``, ``S19.8``, ``W15``, ``W18``
     - Uses pre-computed results from Sukhbold et al. (2016) based on neutrino-powered explosions
   * - **Patton&Sukhbold20-engine**
     -  ``N20``, ``W20``, ``S19.8``, ``W15``, ``W18``
     - Advanced engine combining Patton & Sukhbold (2020) results for realistic explosion landscapes
   * - **Maltsev+25-engine**
     - ``M16``
     - Maltsev et al. (2025) with updated explodability criteria and neutrino-driven physics. **Only valid M_CO < 10 Msun!**
   * - **Maltsev+25-MCO-rapid**
     - ``None``
     - Maltsev et al. (2025) rapid BPS prescription based on the CO core mass and metallicity (arXiv:2503.23856). Supports the ``optimistic``, ``balanced`` and ``pessimistic`` extrapolation modes via ``Maltsev25_MCO_extrapolation_mode`` (default ``balanced``).
   * - **Couch+20-engine**
     - ``"1.0", "1.2", "1.23", "1.25", "1.27", "1.3", "1.4"``
     - Simulations of turbulence-aided neutrino-driven core-collapse supernovae from Couch et al. (2020)
   * - **direct**
     - ``None``
     - Simplified prescriptions that directly collapse the pre-supernova to the baryonic mass
   * - **direct_he_core**
     - ``None``
     - Simplified prescriptions that directly collapse the pre-supernova to the baryonic mass of the helium core



Pulsational pair-instability supernova prescriptions
----------------------------------------------------

POSYDON also supports pulsational pair-instability supernova (PPISN) prescriptions, which can be set in the configuration file:

.. list-table:: (P)PISN prescriptions
   :header-rows: 1
   :widths: 20 10 70

   * - Prescription
     - Parameters
     - Description
   * - **Marchant+19**
     - ``None``
     - Marchant et al. (2019) prescription for PPISN and PISN mass loss from
     | Breivik et al. (2020).
   * - **Hendriks+23**
     - ``PISN_CO_shift`` and ``PPI_extra_mass_loss``
     - Hendriks et al. (2023) prescription for PPISN and PISN mass loss.
     | ``PISN_CO_shift`` shifts the CO core mass threshold for PPI onset
     | and ``PPI_extra_mass_loss`` adds extra mass loss during PPISN events.
     | PISN occurs if the remnant mass after PPI mass loss is less than 10 Msun.



Computational Approaches
------------------------

The outcome of core collapse can be computed using three different approaches, controlled by configuration flags:

Pre-Trained Interpolators
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block::

   use_interp_values = True
   use_profiles = True
   use_core_masses = True

This is the **default and recommended approach** for population synthesis.
The pre-trained interpolators use outcomes that were computed offline from the detailed MESA stellar profiles.
Because the on-the-fly calculations use the down-sampled profiles, the interpolators are generally more accurate and faster.


On-the-Fly Calculations with Downsampled Profiles
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block::

   use_interp_values = False
   use_profiles = True
   use_core_masses = False

This approach recalculates the core-collapse outcome during the simulation using stellar profiles. However, unlike the detailed profiles used to train the interpolators, these profiles are **downsampled** (reduced resolution) for computational efficiency.

**Advantages**: More flexible, can handle evolutionary scenarios not covered by pre-computed grids, allows for model exploration

**Disadvantages**: Slower than interpolators, potential accuracy loss due to profile downsampling

**When to use**: For detailed model comparisons, sensitivity studies, or cases requiring custom physics

Core-Masses Approach (Classical Population Synthesis)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block::

   use_interp_values = False
   use_profiles = False
   use_core_masses = True

This classical approach uses only the core masses at carbon depletion, without requiring detailed stellar profiles.

Pre-Defined Supernova Models
----------------------------

POSYDON provides a set of pre-defined supernova models (``SN_MODEL_v2_XX``) for which
the initial-final interpolator has been trained.
This requires the user to correctly set the supernova mechanism and engine in the configuration file of the population synthesis model.

Model Overview
^^^^^^^^^^^^^^

The following table lists all pre-defined supernova models and their key characteristics:

.. list-table:: Pre-defined SN Models
   :header-rows: 1
   :widths: 12 25 8 18 18

   * - Model
     - Mechanism
     - Engine
     - Conserve H-env
     - PPI Mass Loss

   * - SN_MODEL_v2_01
     - Fryer+12-delayed
     - \-
     - No
     - -20.0

   * - SN_MODEL_v2_02
     - Fryer+12-delayed
     - \-
     - Yes
     - -20.0

   * - SN_MODEL_v2_03
     - Fryer+12-delayed
     - \-
     - No
     - 0.0

   * - SN_MODEL_v2_04
     - Fryer+12-delayed
     - \-
     - Yes
     - 0.0

   * - SN_MODEL_v2_05
     - Fryer+12-rapid
     - \-
     - No
     - -20.0

   * - SN_MODEL_v2_06
     - Fryer+12-rapid
     - \-
     - Yes
     - -20.0

   * - SN_MODEL_v2_07
     - Fryer+12-rapid
     - \-
     - No
     - 0.0

   * - SN_MODEL_v2_08
     - Fryer+12-rapid
     - \-
     - Yes
     - 0.0

   * - SN_MODEL_v2_09
     - Sukhbold+16-engine
     - N20
     - No
     - -20.0

   * - SN_MODEL_v2_10
     - Sukhbold+16-engine
     - N20
     - Yes
     - -20.0

   * - SN_MODEL_v2_11
     - Sukhbold+16-engine
     - N20
     - No
     - 0.0

   * - SN_MODEL_v2_12
     - Sukhbold+16-engine
     - N20
     - Yes
     - 0.0

   * - SN_MODEL_v2_13
     - Patton&Sukhbold20-engine
     - N20
     - No
     - -20.0

   * - SN_MODEL_v2_14
     - Patton&Sukhbold20-engine
     - N20
     - Yes
     - -20.0

   * - SN_MODEL_v2_15
     - Patton&Sukhbold20-engine
     - N20
     - No
     - 0.0

   * - SN_MODEL_v2_16
     - Patton&Sukhbold20-engine
     - N20
     - Yes
     - 0.0

   * - SN_MODEL_v2_17
     - Sukhbold+16-engine
     - W20
     - No
     - -20.0

   * - SN_MODEL_v2_18
     - Sukhbold+16-engine
     - W20
     - Yes
     - -20.0

   * - SN_MODEL_v2_19
     - Sukhbold+16-engine
     - W20
     - No
     - 0.0

   * - SN_MODEL_v2_20
     - Sukhbold+16-engine
     - W20
     - Yes
     - 0.0

   * - SN_MODEL_v2_21
     - Patton&Sukhbold20-engine
     - W20
     - No
     - -20.0

   * - SN_MODEL_v2_22
     - Patton&Sukhbold20-engine
     - W20
     - Yes
     - -20.0

   * - SN_MODEL_v2_23
     - Patton&Sukhbold20-engine
     - W20
     - No
     - 0.0

   * - SN_MODEL_v2_24
     - Patton&Sukhbold20-engine
     - W20
     - Yes
     - 0.0

   * - SN_MODEL_v2_25
     - Maltsev+25-engine
     - M16
     - No
     - -20.0

   * - SN_MODEL_v2_26
     - Maltsev+25-engine
     - M16
     - Yes
     - -20.0

   * - SN_MODEL_v2_27
     - Maltsev+25-engine
     - M16
     - No
     - 0.0

   * - SN_MODEL_v2_28
     - Maltsev+25-engine
     - M16
     - Yes
     - 0.0

   * - SN_MODEL_v2_29
     - Maltsev+25-MCO-rapid
     - ``None`` (balanced)
     - No
     - -20.0

   * - SN_MODEL_v2_30
     - Maltsev+25-MCO-rapid
     - ``None`` (balanced)
     - Yes
     - -20.0

   * - SN_MODEL_v2_31
     - Maltsev+25-MCO-rapid
     - ``None`` (balanced)
     - No
     - 0.0

   * - SN_MODEL_v2_32
     - Maltsev+25-MCO-rapid
     - ``None`` (balanced)
     - Yes
     - 0.0

All pre-defined models use:

- ``use_interp_values = False`` (on-the-fly calculations, not pre-trained interpolators)
- ``use_profiles = True`` (detailed MESA profiles, downsampled during calculation)
- ``use_core_masses = False``

Configuration
^^^^^^^^^^^^^

To use a pre-defined model, specify it in your population synthesis configuration file:

.. code-block:: ini

   [supernova]
   model = SN_MODEL_v2_13

This will automatically load all parameters for that model. Alternatively, you can customize individual parameters:

.. code-block:: ini

   [supernova]
   mechanism = Patton&Sukhbold20-engine
   engine = N20
   ECSN = Tauris+15
   conserve_hydrogen_envelope = False
   PPI_extra_mass_loss = -20.0
   use_interp_values = True
   use_profiles = False
   use_core_masses = False
