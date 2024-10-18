"""Compute the underlying stellar population mass for a given simulation."""


__authors__ = [
    "Devina Misra <devina.misra@unige.ch>",
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Konstantinos Kovlakas <Konstantinos.Kovlakas@unige.ch>",
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>",
    "Dimitris Souropanis <dsouropanis@ia.forth.gr>",
]

import numpy as np
from posydon.popsyn import independent_sample
from scipy.integrate import quad
from posydon.utils.posydonwarning import Pwarn

def initial_total_underlying_mass(simulated_mass=None, simulated_mass_single=None, simulated_mass_binaries=None, f_bin=0.7,  **kwargs):
    """Compute the initial total mass of the population.


     Parameters
    ----------
    simulated_mass, simulated_mass_single, simulated_mass_binaries : float
        Total simulated mass, simulated mass of binary systems and simulated mass of single stars, respectively.
    f_bin: float
        The binary fraction of your population in "nature". If not provided, the default value is set to 0.7

     Parameters are expected to be provided as keyword arguments (kwargs):
     --------
     primary_mass_min: float
         minimum initial mass of the primary stae
     primary_mass_max: float
         maximum initial mass of the primary star
     binary_fraction_const: float
         Binary fraction used in the simulations
     primary_mass_scheme: string
         Kroupa2001 or Salpeter options
     secondary_mass_scheme: string
         mass ratio distribution 
     
    ----------
        

    Returns
    -------
    underlying_total_mass: float
        The underlying total mass of the population: float
    f_corr_single_stars, f_corr_binaries; float
        Correction factors for singles and binaries, respectively.
    
        

    """
    # ENHANCEMENT: the code should assume default values of the POSYDON sampler
    # if not provided in kwargs.
    
    
    f_bin_nature = f_bin #This parameter represents the preferred fraction of binary systems within your population, allowing users to input any value ranging from 0.1 to 1.                               
    f_bin_simulated = kwargs['binary_fraction_const']

    m_min = 0.01
    m_max = 200.0
    m_a = kwargs['primary_mass_min']
    m_b = kwargs['primary_mass_max']

    if simulated_mass is None:
        initial_ZAMS_mass = independent_sample.generate_independent_samples(
            **kwargs)
        initial_ZAMS_mass_1 = initial_ZAMS_mass[2]
        initial_ZAMS_mass_2 = initial_ZAMS_mass[3]
        n_binaries = int(f_bin_simulated * len(initial_ZAMS_mass_1))
        initial_ZAMS_TOTAL_binaries = (sum(initial_ZAMS_mass_1[:n_binaries]) +  sum(initial_ZAMS_mass_2[:n_binaries]))
        initial_ZAMS_TOTAL_single = sum(initial_ZAMS_mass_1[n_binaries:len(initial_ZAMS_mass_1)])
        initial_ZAMS_TOTAL_mass = initial_ZAMS_TOTAL_binaries + initial_ZAMS_TOTAL_single

    else:
        initial_ZAMS_TOTAL_mass = simulated_mass
        initial_ZAMS_TOTAL_single= simulated_mass_single
        initial_ZAMS_TOTAL_binaries= simulated_mass_binaries
    

    def imf_part_1(m, m_min, alpha1):
        return (m/m_min)**-alpha1

    def imf_part_2(m, m_1, m_min, alpha1, alpha2):
        return ((m_1/m_min)**-alpha1) * ((m/m_1)**-alpha2)

    def imf_part_3(m, m_1, m_2, m_min, alpha1, alpha2, alpha3):
        return ((m_1/m_min)**-alpha1)*((m_2/m_1)**-alpha2)*((m/m_2)**-alpha3)

    def mean_mass_of_binary(f0, f_bin_nature, m_1, m_2, m_min, m_max,
                            alpha1, alpha2, alpha3):
        a1, a2, a3 = alpha1, alpha2, alpha3
        mean_mass = f0 * (1 + (f_bin_nature/2)) * (
            (1/(2-a1)) * ((m_1**(2-a1) - m_min**(2-a1))/(m_min**-a1))
            + ((1/(2-a2))
               * (m_1/m_min)**-a1 * ((m_2**(2-a2) - m_1**(2-a2))/(m_1**-a2)))
            + ((1/(2-a3)) * (m_1/m_min)**-a1 * (m_2/m_1)**-a2
               * ((m_max**(2-a3) - m_2**(2-a3))/(m_2**-a3))))
        return mean_mass

    def mean_mass_simulated(alpha3, m_a, m_b, f_bin_simulated):
        a = alpha3
        mean_mass_sim = (1+(f_bin_simulated/2)) * ((1-a) / (2-a)) * (
            (m_b**(2-a)-m_a**(2-a))/(m_b**(1-a)-m_a**(1-a)))
        return mean_mass_sim

    def fraction_simulated( m_1, m_2, m_min, m_max,
                           alpha1, alpha2, alpha3):
        a1, a2, a3 = alpha1, alpha2, alpha3
        f_model =   f0 * (1/(1-a3)) * (m_1/m_min)**(-a1) * (m_2/m_1)**(
            -a2) * ((m_b**(1-a3)-m_a**(1-a3)) / (m_2**-a3)) 
        return f_model

    # Kroupa P., 2001, MNRAS, 322, 231
    if (kwargs['primary_mass_scheme'] == 'Kroupa2001'
            and kwargs['secondary_mass_scheme'] == 'flat_mass_ratio'):
        m_1 = 0.08
        m_2 = 0.5
        alpha1 = 0.3
        alpha2 = 1.3
        alpha3 = 2.3

    # Salpeter E. E., 1955, ApJ, 121, 161
    elif (kwargs['primary_mass_scheme'] == 'Salpeter'
            and kwargs['secondary_mass_scheme'] == 'flat_mass_ratio'):
        m_1 = 0.08
        m_2 = 0.5
        alpha1 = 2.35
        alpha2 = 2.35
        alpha3 = 2.35
    else:
        Pwarn("Scheme not included yet: primary_mass_scheme="
              f"{kwargs['primary_mass_scheme']}, secondary_mass_scheme"
              f"={kwargs['secondary_mass_scheme']}", "UnsupportedModelWarning")
        return np.nan, np.nan, np.nan

    

    f0 = 1/(quad(imf_part_1, m_min, m_1, args=(m_min, alpha1))[0]
                + quad(imf_part_2, m_1, m_2, args=(m_1, m_min, alpha1, alpha2))[0]
                + quad(imf_part_3, m_2, m_max, args=(m_1, m_2, m_min, alpha1, alpha2, alpha3))[0])

    
    f_corr_single_stars = (fraction_simulated(m_1, m_2, m_min, m_max, alpha1, alpha2, alpha3)
                  * mean_mass_simulated(alpha3, m_a, m_b, f_bin_simulated)
                  / mean_mass_of_binary(f0, f_bin_nature, m_1, m_2, m_min, m_max, alpha1, alpha2, alpha3))

    
    f_corr_binaries = (fraction_simulated(m_1, m_2, m_min, m_max, alpha1, alpha2, alpha3)
                  * mean_mass_simulated(alpha3, m_a, m_b, f_bin_simulated)
                  / mean_mass_of_binary(f0, f_bin_nature, m_1, m_2, m_min, m_max, alpha1, alpha2, alpha3))  


    if (f_bin_simulated == 1): #when you have simulated only binary stars
        if (f_bin_nature != 0): 

            f_corr_binaries = f_corr_binaries*f_bin_nature
            f_corr_single_stars=0


            underlying_total_mass=initial_ZAMS_TOTAL_binaries/f_corr_binaries

            
        else:
            Pwarn("We apologize for not having simulated single stars. Please provide a value greater than 0 and less than or equal to 1.",                        "UnsupportedModelWarning")
            
            return np.nan, np.nan, np.nan

    if (f_bin_simulated != 1): #when you have modeled both binary and single stars
        if (f_bin_nature == 0): #but you want the underlying mass for a population consisting of only single stars

            

            f_corr_single_stars= (f_corr_single_stars*(1-f_bin_nature))/(1-f_bin_simulated)
            f_corr_binaries=0
            underlying_total_mass=initial_ZAMS_TOTAL_single/f_corr_single_stars

            
        if (f_bin_nature == 1): #you want the underlying mass for a population consisting of only binary stars
            
            f_corr_binaries = f_corr_binaries*(f_bin_nature/f_bin_simulated)
            f_corr_single_stars=0

            underlying_total_mass=initial_ZAMS_TOTAL_binaries/f_corr_binaries

       
        else:  #you want the underlying mass for a population consisting of both single and binary stars
            
             f_corr_binaries = f_corr_binaries*(f_bin_nature/f_bin_simulated)
             f_corr_single_stars=  (f_corr_single_stars*(1-f_bin_nature))/(1-f_bin_simulated)
             underlying_total_mass=(initial_ZAMS_TOTAL_single/f_corr_single_stars)+(initial_ZAMS_TOTAL_binaries/f_corr_binaries)
       
    
    

    return underlying_total_mass, f_corr_single_stars, f_corr_binaries
     
