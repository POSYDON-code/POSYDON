.. _population_reweighting:

Population Reweighting
======================

POSYDON allows you to reweight your populations after they have been evolved,
based on the initial conditions of the binaries and single stars. This is a
form of importance sampling: the weights are obtained from the ratio of the
PDF of the target population to the PDF of the simulated population, evaluated
at the initial conditions of each model.

The reweighting step is performed on a ``TransientPopulation`` object with the
``get_model_weights`` method. The implementation expects the following initial
properties to be available in the population data:

- ``S1_mass_i``: initial primary mass
- ``S2_mass_i``: initial secondary mass
- ``orbital_period_i``: initial orbital period
- ``state_i``: initial state, which is used to distinguish binaries from
  initially single stars

Supported initial-condition distributions
----------------------------------------

The reweighting implementation in ``posydon.popsyn.norm_pop`` currently supports
several distribution families. The table below summarizes the main options.

.. list-table:: Supported reweighting distributions
   :header-rows: 1
   :widths: 20 25 25 30

   * - Initial condition
     - Supported schemes
     - Required parameters
     - Notes
   * - Primary mass
     - Any IMF class available in ``posydon.popsyn.IMFs``
     - ``primary_mass_scheme``, ``primary_mass_min``, ``primary_mass_max``
     - The scheme name is resolved dynamically from the IMF module. If it is not
       recognized, a flat distribution is used and a warning is emitted.
   * - Mass ratio
     - ``flat_mass_ratio``, ``power_law_mass_ratio``
     - For ``flat_mass_ratio``: ``secondary_mass_scheme`` plus
       ``secondary_mass_min`` and ``secondary_mass_max`` (or explicit
       ``q_min``/``q_max``). For ``power_law_mass_ratio``: ``mass_ratio_slope``
       and optional ``q_min``/``q_max``.
     - The flat mass-ratio distribution can be defined either from the allowed
       secondary-mass range or from explicit mass-ratio bounds.
   * - Orbital period
     - ``period`` with ``Sana+12_period_extended`` or ``power_law``; or
       ``separation`` with ``log_uniform``
     - For ``period``: ``orbital_period_min`` and ``orbital_period_max``. For
       ``Sana+12_period_extended`` no additional parameters are needed. For
       ``power_law``: ``power_law_slope``. For ``separation``:
       ``orbital_separation_min`` and ``orbital_separation_max``.
     - Periods are interpreted in days. Separation-based inputs are converted
       internally to periods using the orbital-separation relation.
   * - Binary fraction
     - ``const``
     - ``binary_fraction_scheme``, ``binary_fraction_const``
     - This is currently implemented as a constant binary fraction and is used
       to distinguish single-star and binary PDFs.

How the weights are computed
---------------------------

The method evaluates the target PDF and the simulation PDF at the initial
conditions of each system and computes

``weights = (PDF_target / PDF_sim) * (mean_mass_sim / mean_mass_pop) * (1 / M_sim)``

where ``M_sim`` is the total simulated mass, ``mean_mass_sim`` is the mean mass
of the simulated population, and ``mean_mass_pop`` is the mean mass of the
requested target population.

In practice, this means that the reweighted population is only meaningful when
the target and simulation populations share the same support for the sampled
initial-condition distributions. If the target distribution extends beyond the
range covered by the simulation, you might miss systems!
