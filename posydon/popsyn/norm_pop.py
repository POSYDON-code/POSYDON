"""Compute the underlying stellar population weight of each model for a given simulation."""

__authors__ = [
    "Devina Misra <devina.misra@unige.ch>",
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Konstantinos Kovlakas <Konstantinos.Kovlakas@unige.ch>",
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>",
    "Dimitris Souropanis <dsouropanis@ia.forth.gr>",
    "Max Briel <max.briel@gmail.com>",
]

import numpy as np
from posydon.popsyn import independent_sample
from scipy.integrate import quad, nquad
from posydon.utils.posydonwarning import Pwarn
from posydon.popsyn.distributions import (flat_mass_ratio,
                                          Sana12Period,
                                          PowerLawPeriod)
import posydon.popsyn.IMFs as IMFs

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
        IMF_pdf = lambda m1: np.ones_like(m1)
        
    return IMF_pdf

def get_mass_ratio_pdf(kwargs):
    """Function that returns the mass ratio PDF function
    
    Supported schemes:
    - `flat_mass_ratio` for `secondary_mass_scheme`
        Requires the following parameters:
        - `secondary_mass_min`
        - `secondary_mass_max`
    
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
                [kwargs['secondary_mass_min'] / m1, np.ones(len(m1))*0.05],
                axis=0)
            
            maximum = np.min(
                [kwargs['secondary_mass_max'] / m1, np.ones(len(m1))],
                axis=0)
            
            q_dist = lambda q: np.where((q >= minimum) & (q <= maximum),
                                        1/(maximum - minimum),
                                        0)
            return q_dist
        q_pdf = lambda q, m1: get_pdf_for_m1(m1)(q)
    elif kwargs['secondary_mass_scheme'] == 'flat_mass_ratio':
        # flat mass ratio, where bounds are given
        q_pdf = lambda q, m1=None: np.where(
                                (q > kwargs['q_min']) & (q <= kwargs['q_max']),
                                1/(kwargs['q_max'] - kwargs['q_min']),
                                0)
    
    else:
        # default to a flat distribution
        Pwarn("The secondary_mass_scheme is not defined use a flat mass ratio "
              "distribution in (0,1].", "UnsupportedModelWarning")
        q_pdf = lambda q, m1=None: np.where((q > 0.0) & (q<=1.0), 1, 0)
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
    if kwargs['orbital_period_scheme'] == 'Sana+12_period_extended':
        period = Sana12Period(
            p_min=kwargs['orbital_period_min'],
            p_max=kwargs['orbital_period_max'],
        )
        period_pdf = lambda P, m1: period.pdf(P, m1)
    elif kwargs['orbital_period_scheme'] == 'power_law':
        period = PowerLawPeriod(
            p_min=kwargs['orbital_period_min'],
            p_max=kwargs['orbital_period_max'],
            pi=kwargs['power_law_slope'],
        )
        period_pdf = lambda P, m1: period.pdf(P)
    else:
        raise ValueError("Orbital period scheme not recognized")
    
    return period_pdf

def get_pdf(kwargs, mass_pdf=False):
    """Function that build a PDF function given the simulation parameters
    
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
                * period_pdf(np.asarray(P), np.asarray(m1))),
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
    
    # binary integration
    # 
    I_bin = nquad(lambda q, m: (m + m * q) * PDF(m, q, P=0, binary=True),
                  ranges=[(q_min, q_max),
                          (m1_min, m1_max)])[0]

    # single star integration
    I_single = nquad(lambda m: m * PDF(m, False),
                     ranges=[(m1_min, m1_max)])[0]
    
    mean_mass = I_bin + I_single
    return mean_mass

def calculate_model_weights(pop_data,
                            M_sim,
                            simulation_parameters,
                            population_parameters):
    '''reweight each model in the simulation to the requested population'''
    
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
    
    return (weight_pop / weight_sim) * factor