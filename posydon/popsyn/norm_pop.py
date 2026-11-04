"""Compute and normalize population synthesis models."""

__authors__ = [
    "Max Briel <max.briel@gmail.com>",
]

import numpy as np
from scipy.integrate import dblquad, quad

from posydon.popsyn.initial_distributions import (
    get_binary_fraction_distribution,
    get_IMF_distribution,
    get_mass_ratio_distribution,
    get_orbital_distribution,
)
from posydon.utils.posydonwarning import Pwarn


def get_pdf(kwargs, mass_pdf=False):
    """Function that builds a PDF function given the simulation parameters

    Parameters
    ----------
    kwargs : dict
        Dictionary containing the simulation parameters
    mass_pdf : bool, optional
        If True, the PDF will return the mass distribution only.
        If False, it will return the full PDF including mass ratio, binary fraction,
        and period distributions. Default is False.

    Returns
    -------
    pdf_function : function
        Function `pdf_function(m1, q=0, P=nan, binary=False)` that returns
        the probability density. The defaults match single stars in the
        population data: no companion (q = 0) and no period (P = nan).
    """

    IMF_pdf = get_IMF_distribution(kwargs).pdf
    q_pdf = get_mass_ratio_distribution(kwargs).pdf
    f_b_pdf = get_binary_fraction_distribution(kwargs).pdf
    period_pdf = get_orbital_distribution(kwargs).pdf_per_dex

    if mass_pdf:
        pdf_function = lambda m1, q=0, P=np.nan, binary=False: (
            np.where(# binaries
                     np.asarray(binary),
                     (f_b_pdf(np.asarray(binary))
                      * IMF_pdf(np.asarray(m1))
                      * q_pdf(np.asarray(q), np.asarray(m1))),
                     # single stars
                     (f_b_pdf(np.asarray(binary))
                      * IMF_pdf(np.asarray(m1)))
                    )
        )
    else:
        pdf_function = lambda m1, q=0, P=np.nan, binary=False: (
            np.where(
                np.asarray(binary),
                # binaries
                (f_b_pdf(np.asarray(binary))
                * IMF_pdf(np.asarray(m1))
                * q_pdf(np.asarray(q), np.asarray(m1))
                * period_pdf(np.asarray(P), np.asarray(m1), np.asarray(q))),
                # single stars
                (f_b_pdf(np.asarray(binary))
                * IMF_pdf(np.asarray(m1)))
            )
        )
    return pdf_function

def get_mean_mass(params):
    '''Calculate the mean mass of the population

    Integrates the mass distribution to calculate the mean mass of
    the population

    Parameters
    ----------
    params : dict
        Dictionary containing the MODEL parameters

    Returns
    -------
    mean_mass : float
        Mean mass of the population
    '''

    PDF = get_pdf(params, mass_pdf=True)

    # integration bounds; the q bounds are the support of the
    # mass ratio distribution, which can depend on m1
    m1_min = params['primary_mass_min']
    m1_max = params['primary_mass_max']
    q_dist = get_mass_ratio_distribution(params)

    # binary integration (quad requires scalar integrands)
    I_bin = dblquad(
        lambda q, m: np.asarray((m + m * q)
                                * PDF(m, q, binary=True)).item(),
        m1_min, m1_max,
        lambda m: float(q_dist.support(m)[0]),
        lambda m: float(q_dist.support(m)[1]),
    )[0]

    # single star integration
    I_single = quad(
        lambda m: np.asarray(m * PDF(m, binary=False)).item(),
        m1_min, m1_max,
    )[0]

    mean_mass = I_bin + I_single
    return mean_mass

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

    Returns
    -------
    output : ndarray of floats
        Weights for each model in the simulation to match the requested population
        This has the units of likelihood of the systems per unit mass (Msun^-1).

    """

    f_b_sim = get_binary_fraction_distribution(simulation_parameters).f_b
    f_b_pop = get_binary_fraction_distribution(population_parameters).f_b
    if (f_b_sim == 1) and (f_b_pop == 0):
        raise ValueError("No single stars simulated, but requested")
    if (f_b_sim == 0) and (f_b_pop == 1):
        raise ValueError("No binaries simulated, but requested")

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

    # Systems with zero probability under the simulation parameters give
    # inf/nan; they lie outside the simulated parameter space and get weight 0.
    with np.errstate(divide='ignore', invalid='ignore'):
        output = (weight_pop / weight_sim) * factor
    outside = ~np.isfinite(output)
    if np.any(outside):
        Pwarn(f"{np.count_nonzero(outside)} of {len(output)} systems have "
              "zero probability under the simulation parameters, i.e. they "
              "lie outside the simulated parameter space. Their weights are "
              "set to 0. Check that the simulation parameters match the "
              "sampled population.", "ValueWarning")
        output[outside] = 0.0
    return output
