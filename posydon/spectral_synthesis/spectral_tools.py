"""

This is tools and useful funtions that are using in the spectral synthesis and
the creation of CMD diagrams combined.
"""

__authors__ = [
    "Eirini Kasdagli <kasdaglie@ufl.edu>",
    "Jeffrey Andrews <jeffrey.andrews@ufl.edu>"]

from copy import copy
import numpy as np
import pandas as pd
from scipy import integrate
import posydon.utils.constants as constants


#Constants
Zo = 0.0142
key_traslation = {
    'log_L',
    'log_R',
}
stripped_stars = ["stripped_He_Central_C_depletion",
                  'stripped_He_Core_He_burning',
                  'stripped_He_Central_He_depleted',
                  'stripped_He_Central_C_depletion',
                  'stripped_He_non_burning']
keys_to_save= ['state',
               'S1_log_g',
               'S2_log_g',
               'S1_Teff',
               'S2_Teff',
               'S1_state',
               'S2_state',
               'S1_mass',
               'S2_mass',
               'S1_log_L',
               'S2_log_L',
               'S1_log_R',
               'S2_log_R',
               'S1_surface_h1',
               'S2_surface_h1',
               'S1_surface_he4',
               'S2_surface_he4',
               'S1_lg_wind_mdot',
               'S2_lg_wind_mdot'
               #'S1_metallicity'
               #'Z/Zo'
               #'Z/Zo',
               #'log_g',
               #'Teff'
                ]
keys_to_load = ['state',
                'time',
               'S1_state',
               'S2_state',
               'S1_mass',
               'S2_mass',
               'S1_log_L',
               'S2_log_L',
               'S1_log_R',
               'S2_log_R',
               'S1_surface_h1',
               'S2_surface_h1',
               'S1_surface_he4',
               'S2_surface_he4',
               'S1_lg_wind_mdot',
               'S2_lg_wind_mdot']


def find_max_time(history):
    """Find the max time of a population."""
    times = history[history.event == 'maxtime'].time
    return np.max(times)

def load_posydon_population(population_file,metallicity):
    """Loads the data of a POSYDON population h5 file saves 
        and add all the data used for spectral synthesis.

    Args:
        population_file: h5 file 
            POSYDON output population file
        max_number_of_binaries: int. Defaults to None.
            The number of binaries as subset of the total population 
            if value is not the default one (Total population). 

    Returns:
        pop: pd array
            The data of the population that are useful for spectral synthesis. 
            Includes the columns keys_to_save as well as metallicity, Teff and log(g) 
            for each star.
    """

    history = pd.read_hdf(population_file, key='history',usecols = keys_to_load)
    max_time = find_max_time(history)
    final_stars = history[(history.time == max_time) & (history.event == 'END')]
    zams_stars = history[(history.time == 0.0)]
    for col in final_stars:
        if col not in keys_to_save:
            del final_stars[str(col)]
    #Find the stripped stars and put their final mass as M_init (mass at ZAMS).
    """
    for star in ['S1','S2']:
        #stripped_index = copy(final_stars[final_stars[f'{star}_state'].isin(stripped_stars)].index.values)
        index_values = copy(final_stars.index.values)
        M_init = zams_stars.loc[index_values][f'{star}_mass'].values
        #final_stars.loc[stripped_index][f'{star}_mass'] = M_init
    """
    pop = copy(final_stars.reset_index())
    #Get the metallicity values from ZAMS
    pop['Z/Zo'] = np.ones(len(pop))*metallicity#zams_stars['S1_metallicity'].iloc[0]/Zo
    #Add two extra columns in the pop data
    for star in ['S1','S2']:
        M = np.asarray(pop[f'{star}_mass'])*constants.Msun
        R = np.asarray(10**pop[f'{star}_log_R'])*constants.Rsun
        L = np.asarray(10**pop[f'{star}_log_L'])*constants.Lsun
        pop[f'{star}_log_g'] = np.log10(constants.standard_cgrav*M/R**2)
        pop[f'{star}_Teff'] = (L/(4*np.pi*R**2*constants.boltz_sigma))**0.25
    return pop


def calculate_logg(row,x):
    """Calculate log_g of S1 and S2 from the population data"""
    mass = row[f'{x}_mass']*constants.Msun
    R = 10**row[f'{x}_log_R']*constants.Rsun
    return np.log10(constants.standard_cgrav*mass/R**2)
 
def calculate_Teff(row,x):
    """Calculate Teff of S1 and S2 from the population data"""
    L = 10**row[f'{x}_log_L']*constants.Lsun
    R = 10**row[f'{x}_log_R']*constants.Rsun
    return (L/(4*np.pi*R**2*constants.boltz_sigma))**0.25

def grid_global_limits(spectral_grids):
    """Calculates the global limits of 
    total collection of the spectral libraries

    Args:
        spectral_grids: dict of the spectral grids

    Returns:
        T_max,T_min,logg_max,logg_min: floats
    """
    T_max  = 0
    T_min = 100000
    logg_max = 0
    logg_min = 20
    for key in spectral_grids:
        specgrid = spectral_grids[key]
        for label in specgrid.axis_labels:
            if label== 'Teff':
                T_max = max(T_max,specgrid.axis_x_max[label])
                T_min = min(T_min,specgrid.axis_x_min[label])
            elif label == 'log(g)':
                logg_max = max(logg_max,specgrid.axis_x_max[label])
                logg_min = min(logg_min,specgrid.axis_x_min[label])
    return T_max,T_min,logg_max,logg_min

def final_mass(file,CO = True):
    history = pd.read_hdf(file,key = 'history')
    final_time = max(history[history.event == 'maxtime'].time)
    end = history[(history.event == 'END')]
    final = end[end.time == final_time]
    if CO:
        S1_CO = end[(end.time != final_time) & (np.isin(end.S1_state,CO_state))]
        S2_CO = end[(end.time != final_time) & (np.isin(end.S2_state,CO_state))]
        total_mass = sum(final.S1_mass) + sum(final.S2_mass) + sum(S1_CO.S1_mass) + sum(S2_CO.S2_mass)
    else: 
        S1_nCO = final[np.isin(final.S1_state,CO_state) == False]
        S2_nCO = final[np.isin(final.S2_state,CO_state) == False]
        total_mass = sum(S1_nCO.S1_mass) + sum(S2_nCO.S2_mass)
    return total_mass


def smooth_flux_negatives(lam_f, flux):
    flux_c = copy(flux)

    i = 0
    i_start = 0
    i_end = 0

    while i < len(flux_c):

        i_low = None
        i_high = None
        low = None
        high = None

        if flux_c[i] > 0 :
             i += 1
             continue

        i -= 1
        # Find point 1
        while low is None:
            if flux_c[i] > 0:
                low = flux_c[i]
                i_low = i
            i += 1
            if i == len(flux_c)-1:
               break
        # Find point 2
        while high is None:
            if flux_c[i] > 0:
                high = flux_c[i]
                i_high = i
            else:
                i += 1
            if i == len(flux_c) -1:
                break
        
        if i == len(flux_c) -1:
                break
        
        i_start = i_low + 1
        i_end = i_high - 1

        # Interpolate between them in log-space
        slope = (np.log10(flux_c[i_high]) - np.log10(flux_c[i_low])) / (lam_f[i_high] - lam_f[i_low])
        intercept = np.log10(flux_c[i_low]) - slope * lam_f[i_low]
        # Replace bad fluxes
        j = copy(i_start)
        
        while j <= i_end:
            flux_c[j] = 10**(slope * lam_f[j] + intercept)
            #print(j, 10**(slope * lam_f[j] + intercept) ,lam_f[j],i_start ,i_low , i_high , slope, intercept)
            j += 1

    return flux_c

def isochrome_weight(m1,IMF_type='Salpeter'):
    #We will normalise the spectrum for a total mass of 1e6 M_sun
    if IMF_type == 'Salpeter':
        a = 2.35
    else:
        raise ValueError("IMF type not compatible")
    m_min = 0.2
    m_max = 120
    norm = ( (-a + 2 ))/(m_max**(-a + 2) - m_min**(-a + 2))
    if m1 > m_max or m1 < m_min: 
        weight = 0
    elif pd.isna(m1):
        weight = 0
    else:
        weight = norm*m1**(-a)
    return weight



def IMF_WEIGHT(mini):
    imf_lower_limit = 0.099  
    imf_upper_limit = 150
    a = 2.35

    wght = [0.0] * len(mini)

    # mirrors: funcint(imf, imf_lower_limit, imf_upper_limit) in Fortran
    # NOTE: no extra x factor — pure number-weighted IMF normalization
    norm = integrate.quad(
        lambda x: x*x**(-a), imf_lower_limit, imf_upper_limit
    )[0]

    for i in range(len(mini)):
        if mini[i] < imf_lower_limit or mini[i] > imf_upper_limit:
            continue

        if i == 0:
            m1 = imf_lower_limit
        else:
            m1 = mini[i] - 0.5 * (mini[i] - mini[i-1])

        if i == len(mini)-1:
            m2 = imf_upper_limit
        else:
            m2 = mini[i] + 0.5 * (mini[i+1] - mini[i])

        if m2 < m1:
            print('IMF_WEIGHT WARNING: non-monotonic mass!', m1, m2, m2-m1)
            continue

        if m2 == m1:
            continue
        dm=0.1
        m_low = max(mini[i] - dm/2, imf_lower_limit)
        m_high = min(mini[i] + dm/2, imf_upper_limit)
        wght[i] = integrate.quad(lambda x: x**(-a), m_low, m_high)[0] / norm

    return wght