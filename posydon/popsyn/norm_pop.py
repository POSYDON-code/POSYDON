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
import posydon.popsyn.IMFs as IMFs

def primary_mass_resample(ZAMS_primary, simulation_parameters, requested_parameters):
    '''Renormalise the population to the requested parameters'''
    
    # TODO: add check if the simulated range is within the requested range!!
    
    # Setup for current simulation weights
    if simulation_parameters['primary_mass_scheme'] == 'Kroupa2001':
        IMF_sim = IMFs.Kroupa2001(m_min=simulation_parameters['primary_mass_min'],
                              m_max=simulation_parameters['primary_mass_max'])
    
    
    # Setup for requested weights
    if requested_parameters['primary_mass_scheme'] == 'Kroupa2001':
        IMF_req = IMFs.Kroupa2001(m_min=requested_parameters['primary_mass_min'],
                                  m_max=requested_parameters['primary_mass_max'])
    
    
    # Calculate the weights    
    P_IMF_sim = IMF_sim.pdf(ZAMS_primary)
    P_IMF_req = IMF_req.pdf(ZAMS_primary)
    
    # New weight of each model in the total population
    f_IMF = P_IMF_req/P_IMF_sim
    
    return f_IMF

def mass_ratio_resample(ZAMS_primary,
                        ZAMS_secondary,
                        simulation_parameters,
                        requested_parameters):
    '''Renormalise the population to the requested parameters'''
    
    q_list = ZAMS_secondary/ZAMS_primary
    
    # Setup for current simulation weights
    if simulation_parameters['secondary_mass_scheme'] == 'flat_mass_ratio':
        m2_min = simulation_parameters['secondary_mass_min']
        m2_max = simulation_parameters['secondary_mass_max']
        mass_ratios_min = np.max([m2_min/ZAMS_primary,
                             np.ones(len(ZAMS_primary))*0.05],
                             axis=0)
        mass_ratios_max = np.min([m2_max/ZAMS_primary,
                              np.ones(len(ZAMS_primary))],
                             axis=0)
        q_dists = [flat_mass_ratio(i,j) for i, j in zip(mass_ratios_min, mass_ratios_max)]
        P_q_sim = np.array([i.pdf(q) for i,q in zip(q_dists, q_list)])
    
    if requested_parameters['secondary_mass_scheme'] == 'flat_mass_ratio':
        # Add check that a q min/max is provided and add that here.
        q_min = 0
        q_max = 1
        q_dist = flat_mass_ratio(q_min, q_max)
        P_q_req = np.array([q_dist.pdf(q) for i,q in zip(q_dists, q_list)])
    
    # Calculate the weights
    f_q = P_q_req/P_q_sim
    
    return f_q


def binary_fraction_resample(ZAMS_primary, simulation_parameters, requested_parameters):
    '''Get the binary fraction correction for the requested parameters'''
    
    # Current the binary fraction is independent of the binary parameters
    # I request ZAMS_primary to get the length of the population
    # this can be used to create a dependent binary fraction
    f_bin_simulated = np.ones(len(ZAMS_primary))*simulation_parameters['binary_fraction_const']
    f_bin_requested = np.ones(len(ZAMS_primary))*requested_parameters['binary_fraction_const']
    
    f_bin = f_bin_requested/f_bin_simulated
    f_sin = (1-f_bin_requested)/(1-f_bin_simulated)
    return f_bin, f_sin
    
def mass_correction(simulation_parameters, requested_parameters):
    
    # Setup for current simulation weights
    # primary mass
    if simulation_parameters['primary_mass_scheme'] == 'Kroupa2001':
        IMF_sim = IMFs.Kroupa2001(m_min=simulation_parameters['primary_mass_min'],
                              m_max=simulation_parameters['primary_mass_max'])
        
    # Setup for requested weights
    if requested_parameters['primary_mass_scheme'] == 'Kroupa2001':
        IMF_req = IMFs.Kroupa2001(m_min=requested_parameters['primary_mass_min'],
                                  m_max=requested_parameters['primary_mass_max'])
    
    if requested_parameters['secondary_mass_scheme'] == 'flat_mass_ratio':
        q_req = flat_mass_ratio(0, 1)
    
    f_b_req = requested_parameters['binary_fraction_const']
    factor = ((1-f_b_req)*quad(lambda m : m*IMF_req.pdf(m), IMF_req.m_min, IMF_req.m_max)[0]
            + f_b_req*nquad(lambda m, q : (m+m*q)*IMF_req.pdf(m)*q_req.pdf(q),
                            ranges=[(IMF_req.m_min, IMF_req.m_max),
                                    (q_req.q_min, q_req.q_max)])[0])
    print('f1:', (1-f_b_req)*quad(lambda m : m*IMF_req.pdf(m), IMF_req.m_min, IMF_req.m_max)[0])
    
    f_b_sim = simulation_parameters['binary_fraction_const']
    factor2 =  ((1-f_b_sim) * quad(lambda m : m*IMF_sim.pdf(m), IMF_sim.m_min, IMF_sim.m_max)[0]
                + f_b_sim * nquad(lambda m, q : (m+m*q)*IMF_sim.pdf(m)*1,
                    ranges=[(IMF_sim.m_min, IMF_sim.m_max),
                            (0.0, 1)])[0])
    print('f2:',(1-f_b_sim) * quad(lambda m : m*IMF_sim.pdf(m), IMF_sim.m_min, IMF_sim.m_max)[0])

    mass_correction = factor/factor2
    return mass_correction
    

def underlying_mass(population, simulation_parameters, requested_parameters):
    
    ZAMS_primary = population['primary_ZAMS']   
    ZAMS_secondary = population['secondary_ZAMS']
    
    mask = population["state_i"] == "initially_single_star"
    single_stars = population[mask]
    binaries = population[~mask]
    
    # correction f_b
    f_bin, f_sin = binary_fraction_resample(ZAMS_primary, simulation_parameters, requested_parameters)

    # correction binaries
    # Only need to input the binary models here (should speed up the calculation!)
    f_IMF_bin = primary_mass_resample(binaries['primary_ZAMS'],
                                  simulation_parameters,
                                  requested_parameters)
    
    # I need to sample the single stars as well
    f_IMF_sin = primary_mass_resample(single_stars['primary_ZAMS'],
                                    simulation_parameters,
                                    requested_parameters)
    
    f_q = mass_ratio_resample(binaries['primary_ZAMS'],
                              binaries['secondary_ZAMS'],
                              simulation_parameters,
                              requested_parameters)
    f_corr_bin = f_IMF_bin
    f_corr_sin = f_IMF_sin
    print(f_corr_sin)
    print(f_sin[mask])
    sample_space_correction = 1/mass_correction(simulation_parameters, requested_parameters)

    underlying_mass_single = single_stars['primary_ZAMS']/(f_corr_sin*(f_sin[mask]))
    underlying_mass_binary = (binaries['primary_ZAMS']+binaries['secondary_ZAMS'])/(f_corr_bin*f_bin[~mask]) 

    underlying_mass = np.sum(underlying_mass_single)/sample_space_correction + np.sum(underlying_mass_binary)/sample_space_correction
    
    return underlying_mass


def prob_q(q_list, ZAMS_primary, m2_min, m2_max):
    
    mass_ratios_min = np.max([m2_min/ZAMS_primary,
                             np.ones(len(ZAMS_primary))*0.05],
                             axis=0)
    mass_ratios_max = np.min([m2_max/ZAMS_primary,
                              np.ones(len(ZAMS_primary))],
                             axis=0)
    q_dists = [flat_mass_ratio(i,j) for i, j in zip(mass_ratios_min, mass_ratios_max)]
    q_prop = [i.pdf(q) for i,q in zip(q_dists, q_list)]
    return np.array(q_prop)