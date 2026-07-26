"""Compute and normalize population synthesis models."""

__authors__ = [
    "Max Briel <max.briel@gmail.com>",
]

import numpy as np
from scipy.integrate import dblquad, nquad, quad
from scipy.interpolate import RegularGridInterpolator

import posydon.popsyn.IMFs as IMFs
from posydon.popsyn import independent_sample
from posydon.popsyn.distributions import (
    FlatMassRatio,
    LogUniform,
    PowerLawMassRatio,
    PowerLawPeriod,
    Sana12Period,
)
from posydon.popsyn.Moes_distributions import Moe_17_PsandQs
from posydon.utils import constants as const
from posydon.utils.common_functions import orbital_separation_from_period
from posydon.utils.posydonwarning import Pwarn

_gen_Moe_17_PsandQs = None


def _get_moe_generator(kwargs):
    """Return a cached Moe & Di Stefano grid generator."""
    global _gen_Moe_17_PsandQs

    if _gen_Moe_17_PsandQs is None:
        _gen_Moe_17_PsandQs = Moe_17_PsandQs(**kwargs)
    return _gen_Moe_17_PsandQs


def _build_moe_pdf_helpers(kwargs):
    """Build interpolation helpers for the Moe & Di Stefano sampler."""

    moe = _get_moe_generator(kwargs)
    m1_grid = moe.M1v
    logp_grid = moe.logPv
    q_grid = moe.qv

    binary_fraction_grid = np.clip(moe.cumPbindist[-1, :], 0.0, 1.0)
    binary_fraction_pdf = RegularGridInterpolator(
        (m1_grid,),
        binary_fraction_grid,
        bounds_error=False,
        fill_value=None,
    )

    period_pdf_grid = np.zeros_like(moe.cumPbindist)
    q_pdf_grid = np.zeros_like(moe.cumqdist)
    q_mean_grid = np.zeros((len(logp_grid), len(m1_grid)))

    for i in range(len(m1_grid)):
        period_pdf_grid[:, i] = np.clip(
            np.gradient(moe.cumPbindist[:, i], logp_grid, edge_order=2),
            0.0,
            None,
        )
        for j in range(len(logp_grid)):
            q_pdf_grid[:, j, i] = np.clip(
                np.gradient(moe.cumqdist[:, j, i], q_grid, edge_order=2),
                0.0,
                None,
            )
            q_mean_grid[j, i] = np.trapz(q_grid * q_pdf_grid[:, j, i], q_grid)

    period_pdf = RegularGridInterpolator(
        (logp_grid, m1_grid),
        period_pdf_grid,
        bounds_error=False,
        fill_value=0.0,
    )
    q_pdf = RegularGridInterpolator(
        (q_grid, logp_grid, m1_grid),
        q_pdf_grid,
        bounds_error=False,
        fill_value=0.0,
    )
    q_mean = RegularGridInterpolator(
        (logp_grid, m1_grid),
        q_mean_grid,
        bounds_error=False,
        fill_value=0.0,
    )

    return moe, binary_fraction_pdf, period_pdf, q_pdf, q_mean


def _moe_primary_pdf(kwargs):
    IMF_pdf = get_IMF_pdf(kwargs)
    moe, binary_fraction_pdf, _, _, _ = _build_moe_pdf_helpers(kwargs)

    def pdf(m1, binary=False):
        m1 = np.asarray(m1)
        binary = np.asarray(binary)
        m1_clipped = np.clip(m1, moe.M1v[0], moe.M1v[-1])
        f_binary = np.clip(binary_fraction_pdf(m1_clipped), 0.0, 1.0)
        return np.where(binary,
                        IMF_pdf(m1) * f_binary,
                        IMF_pdf(m1) * (1.0 - f_binary))

    return pdf


def _moe_full_pdf(kwargs):
    IMF_pdf = get_IMF_pdf(kwargs)
    moe, binary_fraction_pdf, period_pdf, q_pdf, _ = _build_moe_pdf_helpers(kwargs)

    if kwargs.get('orbital_scheme', 'period') != 'period':
        raise ValueError("Moe+17-PsandQs reweighting currently supports only orbital_scheme='period'.")

    def pdf(m1, q=0, P=0, binary=False):
        m1 = np.asarray(m1)
        q = np.asarray(q)
        P = np.asarray(P)
        binary = np.asarray(binary)

        m1_clipped = np.clip(m1, moe.M1v[0], moe.M1v[-1])
        logP = np.where(np.asarray(P) > 0, np.log10(np.asarray(P)), moe.logPv[0])
        logP_clipped = np.clip(logP, moe.logPv[0], moe.logPv[-1])
        q_clipped = np.clip(q, moe.qv[0], moe.qv[-1])
        f_binary = np.clip(binary_fraction_pdf(m1_clipped), 0.0, 1.0)

        binary_density = IMF_pdf(m1) * period_pdf(np.column_stack([logP_clipped, m1_clipped]))
        q_density = q_pdf(np.column_stack([q_clipped, logP_clipped, m1_clipped]))

        return np.where(binary,
                        binary_density * q_density,
                        IMF_pdf(m1) * (1.0 - f_binary))

    return pdf


def get_moe_pdf(kwargs, mass_pdf=False):
    """Return a Moe & Di Stefano initial-condition PDF."""
    if mass_pdf:
        return _moe_primary_pdf(kwargs)
    return _moe_full_pdf(kwargs)


def get_IMF_pdf(kwargs):
    '''get the IMF pdf function

    Supported schemes based on the IMF module:
    Additional parameters can be passed to the scheme
    by adding them to the kwargs dictionary with the scheme name as the key
    Example:

    ```
    kwargs = {
        'primary_mass_scheme': 'Salpeter',
        'primary_mass_min': 0.1,
        'primary_mass_max': 100,
        'Salpeter': {'alpha': 2.35}
    }
    ```

    Parameters
    ----------
    kwargs : dict
        Dictionary containing the parameters.
        `primary_mass_scheme`, `primary_mass_min`, `primary_mass_max` are required.
        Additional parameters are required depending on the scheme.

    Returns
    -------
    IMF_pdf : function
        Function that returns the IMF PDF
    '''

    primary_mass_scheme = kwargs.get('primary_mass_scheme', '')
    scheme_kwargs = kwargs.get(primary_mass_scheme, {})
    try:
        # dynamically retrieve the IMF class from the IMFs module
        imf_class = getattr(IMFs, primary_mass_scheme)
        imf = imf_class(m_min=kwargs['primary_mass_min'],
                        m_max=kwargs['primary_mass_max'],
                        **scheme_kwargs)
        IMF_pdf = imf.pdf
    except AttributeError:
        # if not found, default to a flat distribution
        IMF_pdf = lambda m1: np.ones_like(m1) / (kwargs['primary_mass_max'] - kwargs['primary_mass_min'])
        Pwarn(f"The primary_mass_scheme '{primary_mass_scheme}' is not recognized. "
              "Using a flat mass distribution instead.", "UnsupportedModelWarning")

    return IMF_pdf

def get_mass_ratio_pdf(kwargs):
    """Function that returns the mass ratio PDF function

    Supported schemes:
    - `flat_mass_ratio` for `secondary_mass_scheme`
        Requires the following parameters:
        - `secondary_mass_min`
        - `secondary_mass_max`
    - `power_law_mass_ratio` for `secondary_mass_scheme`
        Requires the following parameters:
        - `mass_ratio_slope`: exponent alpha in q^alpha
        - `q_min` (optional, default 0.05)
        - `q_max` (optional, default 1.0)

    Parameters
    ----------
    kwargs : dict
        Dictionary containing the simulation parameters

    Returns
    -------
    pdf : function
        Function that returns the mass ratio PDF
    """
    if (kwargs['secondary_mass_scheme'] == 'flat_mass_ratio'
        and ('q_min' not in kwargs and 'q_max' not in kwargs)):
        # flat mass ratio, where bounds are dependent on m1 and min/max m2
        # and q_min = 0.05, q_max = 1
        def get_pdf_for_m1(m1):
            m1 = np.atleast_1d(m1)
            minimum = np.max(
                [kwargs['secondary_mass_min'] / m1, np.zeros(len(m1))],
                axis=0)

            maximum = np.min(
                [kwargs['secondary_mass_max'] / m1, np.ones(len(m1))],
                axis=0)

            # Use FlatMassRatio distribution class
            q_dist = lambda q: np.where((q > minimum) & (q <= maximum),
                                        1/(maximum - minimum),
                                        0)
            return q_dist
        q_pdf = lambda q, m1: get_pdf_for_m1(m1)(q)
    elif kwargs['secondary_mass_scheme'] == 'flat_mass_ratio':
        # flat mass ratio, where bounds are given
        from posydon.popsyn.distributions import FlatMassRatio
        q_dist = FlatMassRatio(q_min=kwargs['q_min'], q_max=kwargs['q_max'])
        q_pdf = lambda q, m1=None: q_dist.pdf(q)

    elif kwargs['secondary_mass_scheme'] == 'power_law_mass_ratio':
        from posydon.popsyn.distributions import PowerLawMassRatio
        q_dist = PowerLawMassRatio(
            alpha=kwargs['mass_ratio_slope'],
            q_min=kwargs.get('q_min', 0.05),
            q_max=kwargs.get('q_max', 1.0),
        )
        q_pdf = lambda q, m1=None: q_dist.pdf(q)

    else:
        # default to a flat distribution
        Pwarn("The secondary_mass_scheme is not defined use a flat mass ratio "
              "distribution in (0,1].", "UnsupportedModelWarning")
        from posydon.popsyn.distributions import FlatMassRatio
        q_dist = FlatMassRatio(q_min=0.0, q_max=1.0)
        q_pdf = lambda q, m1=None: q_dist.pdf(q)
    return q_pdf

def get_binary_fraction_pdf(kwargs):
    '''get the binary fraction pdf function

    Supported schemes:
    - `const` with `binary_fraction_const`
        requires the following parameters to be present:
        - `binary_fraction_const`

    Parameters
    ----------
    kwargs : dict
        Dictionary containing the parameters

    Returns
    -------
    pdf : function
        Function that returns the binary fraction PDF
    '''
    if kwargs['binary_fraction_scheme'] == 'const':
        f_b = kwargs['binary_fraction_const']
        binary_fraction_pdf = lambda binary: np.where(np.asarray(binary),
                                                      f_b,
                                                      1-f_b)
    else:
        raise ValueError("Binary fraction scheme not recognized")

    return binary_fraction_pdf

def get_period_pdf(kwargs):
    '''get the period pdf function

    Parameters
    ----------
    kwargs : dict
        Dictionary containing the simulation parameters

    Returns
    -------
    pdf : function
        Function that returns the period PDF, which expects the following
        parameters; P, m1
    '''
    if (kwargs['orbital_scheme'] == 'period'):
        if kwargs['orbital_period_scheme'] == 'Sana+12_period_extended':
            period = Sana12Period(
                p_min=kwargs['orbital_period_min'],
                p_max=kwargs['orbital_period_max'],
            )
            period_pdf = lambda P, m1, q=None: period.pdf(P, m1)
        elif kwargs['orbital_period_scheme'] == 'power_law':
            period = PowerLawPeriod(
                p_min=kwargs['orbital_period_min'],
                p_max=kwargs['orbital_period_max'],
                slope=kwargs['power_law_slope'],
            )
            period_pdf = lambda P, m1=None, q=None: period.pdf(P)
        else:
            raise ValueError("Orbital period scheme not recognized")
    elif (kwargs['orbital_scheme'] == 'separation'):
        if kwargs['orbital_separation_scheme'] == 'log_uniform':
            separation_log_uniform = LogUniform(
                min=kwargs['orbital_separation_min'],
                max=kwargs['orbital_separation_max'],
            )
            # Since the outer function is set in PDF(P) and inner as P(a),
            # we need to do a change of variables.
            # PDF(P) = PDF(a) * |da/dP|
            # da/dP = (2/3) * (a/P)
            # PDF(P) = (2/3) * log_uniform(log_a)
            period_pdf = lambda P, m1, q: separation_log_uniform.pdf(
                orbital_separation_from_period(P, m1, q*m1)
            ) * 2./3.
        else:
            raise ValueError("Orbital separation scheme not recognized")
    else:
        raise ValueError("Orbital scheme not recognized. "
                         "Please specify either 'period' or 'separation' "
                         "and a valid scheme for that orbital scheme.")

    return period_pdf

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
        Function that returns a probability density function
    """

    if independent_sample.use_Moe_17_PsandQs(**kwargs):
        return get_moe_pdf(kwargs, mass_pdf=mass_pdf)

    IMF_pdf = get_IMF_pdf(kwargs)
    q_pdf = get_mass_ratio_pdf(kwargs)
    f_b_pdf = get_binary_fraction_pdf(kwargs)
    period_pdf = get_period_pdf(kwargs)

    if mass_pdf:
        pdf_function = lambda m1, q=0, P=0, binary=False: (
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
        pdf_function = lambda m1, q=0, P=0, binary=False: (
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

    if independent_sample.use_Moe_17_PsandQs(**params):
        moe, binary_fraction_pdf, _, _, q_mean = _build_moe_pdf_helpers(params)
        IMF_pdf = get_IMF_pdf(params)

        m1_grid = moe.M1v
        logp_grid = moe.logPv

        m1_clipped = np.clip(m1_grid, params['primary_mass_min'], params['primary_mass_max'])
        f_binary = np.clip(binary_fraction_pdf(m1_clipped), 0.0, 1.0)

        binary_mass_at_m1 = np.zeros_like(m1_grid)
        for i, m1 in enumerate(m1_grid):
            p_bin_logp = np.clip(
                np.gradient(moe.cumPbindist[:, i], logp_grid, edge_order=2),
                0.0,
                None,
            )
            q_expect = q_mean(np.column_stack([logp_grid, np.full_like(logp_grid, m1)]))
            binary_mass_at_m1[i] = np.trapz(p_bin_logp * m1 * (1.0 + q_expect), logp_grid)

        single_mass_at_m1 = m1_grid * (1.0 - f_binary)
        mean_mass = np.trapz(
            IMF_pdf(m1_grid) * (single_mass_at_m1 + binary_mass_at_m1),
            m1_grid,
        )
        return mean_mass

    PDF = get_pdf(params, mass_pdf=True)

    # integration bounds
    m1_min = params['primary_mass_min']
    m1_max = params['primary_mass_max']

    if 'q_min' in params:
        q_min = params['q_min']
    else:
        q_min = np.max([params['secondary_mass_min']/params['primary_mass_min'],
                        0])

    if 'q_max' in params:
        q_max = params['q_max']
    else:
        q_max = np.min([params['secondary_mass_max']/params['primary_mass_max'],
                        1])
    if q_min > q_max:
        raise ValueError("q_min must be less than q_max")

    # binary integration
    I_bin = dblquad(lambda q, m: (m + m * q) * PDF(m, q, P=0, binary=True),
                    m1_min, m1_max,
                    q_min, q_max
                    )[0]


    # single star integration
    I_single = quad(lambda m: m * PDF(m, False), m1_min, m1_max)[0]

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
        - `secondary_mass_min`
        - `secondary_mass_max`
        - `binary_fraction_scheme`
        - `binary_fraction_const`
        - `orbital_scheme`
        - `orbital_period_scheme` or `orbital_separation_scheme` depending on the `orbital_scheme`
        - `orbital_period_min` and `orbital_period_max` or `orbital_separation_min` and `orbital_separation_max` depending on the `orbital_scheme`
        - `power_law_slope` if `orbital_period_scheme` is `power_law`
        - `q_min` and `q_max` if `secondary_mass_scheme` is `flat_mass_ratio`

    population_parameters : dict
        Dictionary containing the population parameters, which is the requested population to which we want to reweight the simulation. This is used to calculate the PDF of the requested population.
        The parameters in this dictionary are the initial conditions of the population you want to reweight to.
        The following parameters are required to be present in the dictionary:
        - `primary_mass_scheme`
        - `primary_mass_min`
        - `primary_mass_max`
        - `secondary_mass_scheme`
        - `secondary_mass_min`
        - `secondary_mass_max`
        - `binary_fraction_scheme`
        - `binary_fraction_const`
        - `orbital_scheme`
        - `orbital_period_scheme` or `orbital_separation_scheme` depending on the `orbital_scheme`
        - `orbital_period_min` and `orbital_period_max` or `orbital_separation_min` and `orbital_separation_max` depending on the `orbital_scheme`
        - `power_law_slope` if `orbital_period_scheme` is `power_law`
        - `q_min` and `q_max` if `secondary_mass_scheme` is `flat_mass_ratio`

    Returns
    -------
    output : ndarray of floats
        Weights for each model in the simulation to match the requested population
        This has the units of likelihood of the systems per unit mass (Msun^-1).

    """

    sim_is_moe = independent_sample.use_Moe_17_PsandQs(**simulation_parameters)
    pop_is_moe = independent_sample.use_Moe_17_PsandQs(**population_parameters)

    if not sim_is_moe and not pop_is_moe:
        f_b_sim = simulation_parameters['binary_fraction_const']
        f_b_pop = population_parameters['binary_fraction_const']
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

    output = (weight_pop / weight_sim) * factor
    output[np.isnan(output)] = 0.0
    output[np.isinf(output)] = 0.0
    return output
