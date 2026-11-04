"""Compute and normalize population synthesis models."""

__authors__ = [
    "Max Briel <max.briel@gmail.com>",
]

import inspect

import numpy as np
from scipy.integrate import dblquad, quad

import posydon.popsyn.IMFs as IMFs
from posydon.popsyn.distributions import (
    ConditionalFlatMassRatio,
    FlatMassRatio,
    LogUniform,
    PowerLawMassRatio,
    PowerLawPeriod,
    Sana12Period,
)
from posydon.utils.common_functions import orbital_separation_from_period


def get_IMF_pdf(kwargs):
    '''get the IMF pdf function

    Supported schemes are the (non-abstract) IMF classes in the IMFs module,
    e.g. `Salpeter`, `Kroupa1993`, `Kroupa2001`, `Chabrier2003`.
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

    Raises
    ------
    ValueError
        If `primary_mass_scheme` is not an IMF class of the IMFs module.
    '''

    primary_mass_scheme = kwargs.get('primary_mass_scheme', '')
    # dynamically retrieve the IMF class from the IMFs module
    imf_class = (getattr(IMFs, primary_mass_scheme, None)
                 if isinstance(primary_mass_scheme, str) else None)
    if not _is_IMF_class(imf_class):
        supported = sorted(name for name, obj in vars(IMFs).items()
                           if _is_IMF_class(obj))
        raise ValueError(f"The primary_mass_scheme '{primary_mass_scheme}' "
                         f"is not recognized. Supported schemes: {supported}")

    scheme_kwargs = kwargs.get(primary_mass_scheme, {})
    imf = imf_class(m_min=kwargs['primary_mass_min'],
                    m_max=kwargs['primary_mass_max'],
                    **scheme_kwargs)
    return imf.pdf

def _is_IMF_class(obj):
    '''Check whether obj is a concrete IMF class from the IMFs module.'''
    return (inspect.isclass(obj)
            and issubclass(obj, IMFs.IMFBase)
            and not inspect.isabstract(obj))

def get_mass_ratio_distribution(kwargs):
    """Build the mass ratio distribution from the simulation parameters

    The same distribution classes are used by the population sampler, so
    the PDF, its support and the sampled mass ratios agree.

    Supported schemes:
    - `flat_mass_ratio` for `secondary_mass_scheme`
        Requires either of the following parameter sets:
        - `secondary_mass_min`, `secondary_mass_max` and `primary_mass_min`
          (q bounds depend on m1, as in the sampler)
        - `q_min` and `q_max` (fixed q bounds)
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
    q_dist : MassRatioDistribution
        Distribution with `pdf(q, m1)` and `support(m1)`.

    Raises
    ------
    KeyError
        If only one of `q_min` and `q_max` is given for `flat_mass_ratio`.
    ValueError
        If the `secondary_mass_scheme` is not supported (such as the
        sampler's `q=1`), if `primary_mass_min` is below
        `secondary_mass_min` (as in the sampler), or if the bounds are
        invalid.
    """
    scheme = kwargs.get('secondary_mass_scheme')
    if scheme == 'flat_mass_ratio':
        if 'q_min' in kwargs or 'q_max' in kwargs:
            # fixed bounds; a KeyError is raised if only one is given
            return FlatMassRatio(q_min=kwargs['q_min'], q_max=kwargs['q_max'])
        # bounds depend on m1, as in the sampler
        if kwargs['primary_mass_min'] < kwargs['secondary_mass_min']:
            raise ValueError("`secondary_mass_min` is larger than "
                             "`primary_mass_min`.")
        return ConditionalFlatMassRatio(m2_min=kwargs['secondary_mass_min'],
                                        m2_max=kwargs['secondary_mass_max'])
    elif scheme == 'power_law_mass_ratio':
        return PowerLawMassRatio(alpha=kwargs['mass_ratio_slope'],
                                 q_min=kwargs.get('q_min', 0.05),
                                 q_max=kwargs.get('q_max', 1.0))
    else:
        raise ValueError(f"The secondary_mass_scheme '{scheme}' is not "
                         "supported for normalisation. Supported schemes: "
                         "['flat_mass_ratio', 'power_law_mass_ratio']")

def get_mass_ratio_pdf(kwargs):
    """Function that returns the mass ratio PDF function

    See `get_mass_ratio_distribution` for the supported schemes.

    Parameters
    ----------
    kwargs : dict
        Dictionary containing the simulation parameters

    Returns
    -------
    pdf : function
        Function `pdf(q, m1=None)` that returns the mass ratio PDF
    """
    return get_mass_ratio_distribution(kwargs).pdf

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

    PDF = get_pdf(params, mass_pdf=True)

    # integration bounds; the q bounds are the support of the
    # mass ratio distribution, which can depend on m1
    m1_min = params['primary_mass_min']
    m1_max = params['primary_mass_max']
    q_dist = get_mass_ratio_distribution(params)

    # binary integration (quad requires scalar integrands)
    I_bin = dblquad(
        lambda q, m: np.asarray((m + m * q)
                                * PDF(m, q, P=0, binary=True)).item(),
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
        - `power_law_slope` if `orbital_period_scheme` is `power_law`

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
        - `power_law_slope` if `orbital_period_scheme` is `power_law`

    Returns
    -------
    output : ndarray of floats
        Weights for each model in the simulation to match the requested population
        This has the units of likelihood of the systems per unit mass (Msun^-1).

    """

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
