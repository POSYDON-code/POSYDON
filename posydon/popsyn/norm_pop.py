"""Compute and normalize population synthesis models."""

__authors__ = [
    "Max Briel <max.briel@gmail.com>",
]

import numpy as np
from scipy.integrate import quad

from posydon.popsyn.initial_distributions import (
    get_binary_fraction_distribution,
    get_IMF_distribution,
    get_mass_ratio_distribution,
    get_orbital_distribution,
)
from posydon.utils.posydonwarning import Pwarn


def get_pdf(kwargs, mass_pdf=False):
    """Function that builds a PDF function given the simulation parameters

    The PDF is

    - binaries: f_b(m1) * IMF(m1) * p(q | m1, P) * p(log10 P | m1, q)
    - single stars: (1 - f_b(m1)) * IMF(m1)

    The orbital part is a density per dex in the period, dN/dlog10(P), for
    every orbital scheme; separation distributions are converted to periods
    with `PeriodFromSeparation`. Using the same measure for all schemes lets
    the PDFs of the simulation and the target population be compared.

    For Moe & Di Stefano (2017), the mass ratio depends on the period, and
    the period distribution is normalised to 1 for binaries, so the product
    is `Moe_17_PsandQs.pdf` and the binary fraction is counted once.

    Parameters
    ----------
    kwargs : dict
        Dictionary containing the simulation parameters
    mass_pdf : bool, optional
        If True, the PDF will return the mass distribution only.
        If False, it will return the full PDF including mass ratio, binary fraction,
        and period distributions. Default is False. Not available if the mass
        ratio depends on the period (Moe & Di Stefano 2017).

    Returns
    -------
    pdf_function : function
        Function `pdf_function(m1, q=0, P=nan, binary=False)` that returns
        the probability density. The defaults match single stars in the
        population data: no companion (q = 0) and no period (P = nan).
    """

    IMF_pdf = get_IMF_distribution(kwargs).pdf
    q_dist = get_mass_ratio_distribution(kwargs)
    f_b_pdf = get_binary_fraction_distribution(kwargs).pdf
    period_pdf = get_orbital_distribution(kwargs).pdf_per_dex

    if mass_pdf and q_dist.depends_on_period:
        raise ValueError("The mass ratio distribution depends on the "
                         "period; use mass_pdf=False.")

    def pdf_function(m1, q=0, P=np.nan, binary=False):
        m1, q, P, binary = np.broadcast_arrays(
            np.asarray(m1, dtype=float), np.asarray(q, dtype=float),
            np.asarray(P, dtype=float), np.asarray(binary, dtype=bool))
        # the binary part is also evaluated for single stars (q = 0,
        # P = nan), where it is not used
        with np.errstate(invalid='ignore', divide='ignore'):
            orbit = q_dist.pdf(q, m1, P)
            if not mass_pdf:
                orbit = orbit * period_pdf(P, m1, q)
        return f_b_pdf(binary, m1) * IMF_pdf(m1) * np.where(binary, orbit, 1.0)

    return pdf_function

def get_mean_mass(params):
    '''Calculate the mean mass of the population

    Integrates the mass distribution to calculate the mean (initial) mass of
    a system, single or binary:

        mean_mass = int IMF(m1) m1 (1 + f_b(m1) E[q | m1]) dm1

    where E[q | m1] is the mean of the mass ratio distribution, over its
    support (`MassRatioDistribution.mean`).

    Parameters
    ----------
    params : dict
        Dictionary containing the MODEL parameters

    Returns
    -------
    mean_mass : float
        Mean mass of the population
    '''
    IMF_pdf = get_IMF_distribution(params).pdf
    f_b_dist = get_binary_fraction_distribution(params)
    q_dist = get_mass_ratio_distribution(params)

    m1_min = params['primary_mass_min']
    m1_max = params['primary_mass_max']
    points = [m for m in getattr(f_b_dist, 'breakpoints', ())
              if m1_min < m < m1_max] or None

    # quad requires scalar integrands
    def integrand(m1):
        return np.asarray(IMF_pdf(m1) * m1
                          * (1 + f_b_dist.binary_fraction(m1)
                             * q_dist.mean(m1))).item()

    return quad(integrand, m1_min, m1_max, points=points, limit=200)[0]

def check_binary_fraction_coverage(simulation_parameters,
                                   population_parameters, n_m1=2000):
    """Check that the simulation covers the binary fraction of the population

    The check is done on a grid of primary masses in the overlap of the two
    primary-mass ranges. Where the population requests only binaries but the
    simulation has none (f_b,pop = 1, f_b,sim = 0), or only single stars but
    the simulation has no single stars (f_b,pop = 0, f_b,sim = 1), no system
    can be reweighted and an error is raised. Where the population requests
    binaries or single stars of which the simulation has none, while the
    other kind is covered, those systems are missing from the weights and a
    warning is issued.

    Parameters
    ----------
    simulation_parameters, population_parameters : dict
        Parameters of the simulation and of the requested population.
    n_m1 : int, optional
        Number of grid points in log10(m1).

    Raises
    ------
    ValueError
        If, for some primary masses, the simulation has no system of the kind
        requested by the population.
    """
    f_b_sim = get_binary_fraction_distribution(simulation_parameters)
    f_b_pop = get_binary_fraction_distribution(population_parameters)
    m_lo = max(simulation_parameters['primary_mass_min'],
               population_parameters['primary_mass_min'])
    m_hi = min(simulation_parameters['primary_mass_max'],
               population_parameters['primary_mass_max'])
    if m_lo > m_hi:
        return
    m1 = np.geomspace(m_lo, m_hi, n_m1)
    f_sim = f_b_sim.binary_fraction(m1)
    f_pop = f_b_pop.binary_fraction(m1)
    no_bin = (f_sim == 0) & (f_pop > 0)
    no_single = (f_sim == 1) & (f_pop < 1)

    def mass_range(mask):
        return f"{m1[mask].min():.3g}-{m1[mask].max():.3g} Msun"

    for missing, uncovered, kind in [(no_single, f_pop == 0, "single stars"),
                                     (no_bin, f_pop == 1, "binaries")]:
        if np.any(missing & uncovered):
            raise ValueError(f"No {kind} simulated, but requested for "
                             "primary masses "
                             f"{mass_range(missing & uncovered)}.")
        if np.any(missing):
            Pwarn(f"No {kind} simulated, but requested for primary masses "
                  f"{mass_range(missing)}; they are missing from the "
                  "weights.", "InappropriateValueWarning")

def calculate_model_weights(pop_data,
                            M_sim,
                            simulation_parameters,
                            population_parameters):
    """Reweight each model in the simulation to the requested population

    Uses the PDF of the simulation and the PDF of the requested population to calculate
    the weights for each model in the simulation to match the requested population.

    Parameters
    ----------
    pop_data : dict
        Dictionary containing the population data.
        This needs to contain the following keys:
        - `S1_mass_i`: initial mass of the primary
        - `S2_mass_i`: initial mass of the secondary
        - `orbital_period_i`: initial orbital period
        - `state_i`: initial state of the system (e.g. 'initially_single_star' for single stars)
        These are used to calculate the PDF for each model in the simulation.

    M_sim : float
        Mass of the simulation
    simulation_parameters : dict
        Dictionary containing the simulation parameters.
        This is used to calculate the PDF of the simulation.
        The parameters in this dictionary are the initial conditions of the population.
        The following parameters are required to be present in the dictionary:
        - `primary_mass_scheme`
        - `primary_mass_min`
        - `primary_mass_max`
        - `secondary_mass_scheme`
        - `secondary_mass_min` and `secondary_mass_max`, or `q_min` and
          `q_max`, if `secondary_mass_scheme` is `flat_mass_ratio`
        - `mass_ratio_slope` (and optionally `q_min`, `q_max`) if
          `secondary_mass_scheme` is `power_law_mass_ratio`
        - `binary_fraction_scheme`
        - `binary_fraction_const`
        - `orbital_scheme`
        - `orbital_period_scheme` or `orbital_separation_scheme` depending on the `orbital_scheme`
        - `orbital_period_min` and `orbital_period_max` or `orbital_separation_min` and `orbital_separation_max` depending on the `orbital_scheme`
        - `power_law_slope` if `orbital_period_scheme` is `power_law` or `log_power_law`
        - `log_orbital_separation_mean` and `log_orbital_separation_sigma` if `orbital_separation_scheme` is `log_normal`
        For Moe & Di Stefano (2017) (`secondary_mass_scheme`,
        `orbital_period_scheme` or `eccentricity_scheme` set to
        'Moe+17-PsandQs'), the binary fraction, mass ratio and period come
        from `Moe_17_PsandQs`; only `orbital_period_min` and
        `orbital_period_max` are needed besides the primary mass settings.

    population_parameters : dict
        Dictionary containing the population parameters, which is the requested population to which we want to reweight the simulation. This is used to calculate the PDF of the requested population.
        The parameters in this dictionary are the initial conditions of the population you want to reweight to.
        The following parameters are required to be present in the dictionary:
        - `primary_mass_scheme`
        - `primary_mass_min`
        - `primary_mass_max`
        - `secondary_mass_scheme`
        - `secondary_mass_min` and `secondary_mass_max`, or `q_min` and
          `q_max`, if `secondary_mass_scheme` is `flat_mass_ratio`
        - `mass_ratio_slope` (and optionally `q_min`, `q_max`) if
          `secondary_mass_scheme` is `power_law_mass_ratio`
        - `binary_fraction_scheme`
        - `binary_fraction_const`
        - `orbital_scheme`
        - `orbital_period_scheme` or `orbital_separation_scheme` depending on the `orbital_scheme`
        - `orbital_period_min` and `orbital_period_max` or `orbital_separation_min` and `orbital_separation_max` depending on the `orbital_scheme`
        - `power_law_slope` if `orbital_period_scheme` is `power_law` or `log_power_law`
        - `log_orbital_separation_mean` and `log_orbital_separation_sigma` if `orbital_separation_scheme` is `log_normal`
        For Moe & Di Stefano (2017) (`secondary_mass_scheme`,
        `orbital_period_scheme` or `eccentricity_scheme` set to
        'Moe+17-PsandQs'), the binary fraction, mass ratio and period come
        from `Moe_17_PsandQs`; only `orbital_period_min` and
        `orbital_period_max` are needed besides the primary mass settings.

    Returns
    -------
    output : ndarray of floats
        Weights for each model in the simulation to match the requested population
        This has the units of likelihood of the systems per unit mass (Msun^-1).

    """

    check_binary_fraction_coverage(simulation_parameters,
                                   population_parameters)

    # build the pdf functions
    PDF_sim = get_pdf(simulation_parameters, mass_pdf=False)
    PDF_pop = get_pdf(population_parameters, mass_pdf=False)

    # initial properties
    mean_mass_sim = get_mean_mass(simulation_parameters)
    mean_mass_pop = get_mean_mass(population_parameters)

    factor = (1/M_sim) * (mean_mass_sim / mean_mass_pop)

    # we still need to distinguish between binary and single stars for the PDF
    binary_mask = pop_data['state_i'] != 'initially_single_star'
    weight_pop = PDF_pop(m1=pop_data['S1_mass_i'],
                         q=pop_data['S2_mass_i']/pop_data['S1_mass_i'],
                         P=pop_data['orbital_period_i'],
                         binary=binary_mask)

    weight_sim = PDF_sim(m1=pop_data['S1_mass_i'],
                         q=pop_data['S2_mass_i']/pop_data['S1_mass_i'],
                         P=pop_data['orbital_period_i'],
                         binary=binary_mask)

    # Systems outside the sampled parameter space have zero probability
    # under the simulation parameters (inf/nan); they get weight 0.
    with np.errstate(divide='ignore', invalid='ignore'):
        output = (weight_pop / weight_sim) * factor
    outside = ~np.isfinite(output)
    if np.any(outside):
        Pwarn(f"{np.count_nonzero(outside)} of {len(output)} systems lie "
              "outside the sampled parameter space. This should not happen; "
              "please check that this is expected. Their weights are set "
              "to 0.", "ValueWarning")
        output[outside] = 0.0
    return output
