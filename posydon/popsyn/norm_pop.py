"""Compute the underlying stellar population weight of each model for a given simulation."""

__authors__ = [
    "Devina Misra <devina.misra@unige.ch>",
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Konstantinos Kovlakas <Konstantinos.Kovlakas@unige.ch>",
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>",
    "Dimitris Souropanis <dsouropanis@ia.forth.gr>",
    "MMax Briel <max.briel@gmail.com>",
]

import numpy as np
from posydon.popsyn import independent_sample
from scipy.integrate import quad, nquad
from posydon.utils.posydonwarning import Pwarn
from posydon.popsyn.distributions import flat_mass_ratio
from functools import lru_cache
import posydon.popsyn.IMFs as IMFs

def get_IMF_pdf(kwargs):
    '''get the IMF pdf function'''
    
    primary_mass_scheme = kwargs.get('primary_mass_scheme', '')
    try:
        # dynamically retrieve the IMF class from the IMFs module
        imf_class = getattr(IMFs, primary_mass_scheme)
        imf = imf_class(m_min=kwargs['primary_mass_min'],
                        m_max=kwargs['primary_mass_max'])
        IMF_pdf = lambda m1: imf.pdf(m1)
    except AttributeError:
        # if not found, default to a flat distribution
        IMF_pdf = lambda m1: 1
        
    return IMF_pdf

def get_mass_ratio_pdf(kwargs):
    """Function that returns the mass ratio PDF function with caching for optimization."""
    if kwargs['secondary_mass_scheme'] == 'flat_mass_ratio' and ('q_min' not in kwargs and 'q_max' not in kwargs):
        # flat mass ratio, where bounds are dependent on m1 and min/max m2
        # and q_min = 0.05, q_max = 1
        def get_pdf_for_m1(m1):
            m1 = np.atleast_1d(m1)
            minimum = np.max([kwargs['secondary_mass_min'] / m1, np.ones(len(m1))*0.05], axis=0)
            maximum = np.min([kwargs['secondary_mass_max'] / m1, np.ones(len(m1))], axis=0)
            q_dist = lambda q: np.where((q >= minimum) & (q <= maximum), 1/(maximum - minimum), 0)
            return q_dist
        q_pdf = lambda q, m1: get_pdf_for_m1(m1)(q)
    else:
        q_pdf = lambda q, m1=None: np.where((q > 0.0) & (q<=1.0), 1, 0)
    return q_pdf

def get_pdf(kwargs):
    """Function that build a PDF function given the simulation parameters"""
    
    IMF_pdf = get_IMF_pdf(kwargs)
    q_pdf = get_mass_ratio_pdf(kwargs)
    
    f_b = kwargs['binary_fraction_const']
    
    pdf_function = lambda m1, q=0, binary=False: np.where(
        np.asarray(binary),
        f_b * IMF_pdf(np.asarray(m1)) * q_pdf(np.asarray(q), np.asarray(m1)),
        (1-f_b) * IMF_pdf(np.asarray(m1))
    )
    return pdf_function


def get_mean_mass(PDF, params):
    '''Calculate the mean mass of the population'''
    
    # integration bounds
    m1_min = params['primary_mass_min']
    m1_max = params['primary_mass_max']
    
    if 'q_min' in params:
        q_min = params['q_min']
    else:
        q_min = np.max([params['secondary_mass_min']/params['primary_mass_min'], 0])
        
    if 'q_max' in params:
        q_max = params['q_max']
    else:
        q_max = np.min([params['secondary_mass_max']/params['primary_mass_max'], 1])

    # binary integration
    I_bin = nquad(lambda q, m: (m + m * q) * PDF(m, q, True),
                  ranges=[(q_min, q_max),
                          (m1_min, m1_max)])[0]

    # single star integration
    I_single = nquad(lambda m: m * PDF(m, False),
                     ranges=[(m1_min, m1_max)])[0]
    
    mean_mass = I_bin + I_single
    return mean_mass


# not required anymore
def calculate_underlying_mass(population, simulation_parameters, requested_parameters):
    """Calculate the underlying mass of the population"""
    
    # build the pdf functions
    PDF_sim = get_pdf(simulation_parameters, simulation=True)
    PDF_pop = get_pdf(requested_parameters)
    
    # initial properties
    single_mask = population['state_i'] == 'initially_single_star'
    single_stars = population[single_mask].copy()

    binary_stars = population[~single_mask].copy()
    binary_stars['mass_ratio_ZAMS'] = binary_stars['secondary_ZAMS']/binary_stars['primary_ZAMS']
    
    simulated_mass = np.sum(population['primary_ZAMS']) + np.sum(population['secondary_ZAMS'])
    
    f_b_sim = simulation_parameters['binary_fraction_const']
    f_b_pop = requested_parameters['binary_fraction_const']
    
    if (f_b_sim == 1) and (f_b_pop == 0):
        raise ValueError("No single stars simulated, but requested")
    if (f_b_sim == 0) and (f_b_pop == 1):
        raise ValueError("No binaries simulated, but requested")
        
    # counts
    N_sim = len(population)
    
    # binaries
    if f_b_pop == 0:
        binary_weights = np.zeros(len(binary_stars))
    elif f_b_sim == 0:
        binary_weights = np.zeros(len(binary_stars))
    else:
        binary_weights = f_b_sim/f_b_pop * (PDF_sim(m1=binary_stars['primary_ZAMS'],
                           q=binary_stars['mass_ratio_ZAMS'],
                           binary=True)
                   / PDF_pop(m1=binary_stars['primary_ZAMS'],
                             q=binary_stars['mass_ratio_ZAMS'],
                             binary=True))
    # single stars
    if f_b_pop == 1: # no single stars requested
        single_weights = np.zeros(len(single_stars))
    elif f_b_sim == 1: # no single stars simulated
        single_weights = np.zeros(len(single_stars))
    else:
        single_weights = (1-f_b_sim)/(1-f_b_pop) * (PDF_sim(m1=single_stars['primary_ZAMS'],
                                                            binary=False)
                                             / PDF_pop(m1=single_stars['primary_ZAMS'],
                                                         binary=False))
    N_pop = np.nansum(binary_weights) + np.nansum(single_weights)

    
    print('N_sim', N_sim)
    print('N_pop', N_pop)
    
    # mean masses
    mean_mass_sim = get_mean_mass(PDF_sim, simulation_parameters, simulation=True)
    mean_mass_pop = get_mean_mass(PDF_pop, requested_parameters)
    print('sim mass:', mean_mass_sim)
    print('pop mass:', mean_mass_pop)
    
    M_pop = simulated_mass * (N_pop * mean_mass_pop)/ (mean_mass_sim * N_sim)
    
    return M_pop
    

def calculate_model_weights(pop_data, M_sim, simulation_parameters, population_parameters):
    '''reweight each model in the simulation to the requested population'''
    
    f_b_sim = simulation_parameters['binary_fraction_const']
    f_b_pop = population_parameters['binary_fraction_const']
    if (f_b_sim == 1) and (f_b_pop == 0):
        raise ValueError("No single stars simulated, but requested")
    if (f_b_sim == 0) and (f_b_pop == 1):
        raise ValueError("No binaries simulated, but requested")
    
    # build the pdf functions
    PDF_sim = get_pdf(simulation_parameters)
    PDF_pop = get_pdf(population_parameters)
    
    # initial properties
    mean_mass_sim = get_mean_mass(PDF_sim, simulation_parameters)
    mean_mass_pop = get_mean_mass(PDF_pop, population_parameters)
        
    factor = (1/M_sim) * (mean_mass_sim / mean_mass_pop)
    
    # we still need to distinguish between binary and single stars for the PDF
    binary_mask = pop_data['state_i'] != 'initially_single_star'
    weight_pop = PDF_pop(m1=pop_data['S1_mass_i'], q=pop_data['S2_mass_i']/pop_data['S1_mass_i'], binary=binary_mask)
    weight_sim = PDF_sim(m1=pop_data['S1_mass_i'], q=pop_data['S2_mass_i']/pop_data['S1_mass_i'], binary=binary_mask)
    
    return (weight_pop / weight_sim) * factor