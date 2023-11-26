"""

This is tools and useful funtions that are using in the spectral synthesis and
the creation of CMD diagrams combined.
"""

__authors__ = [
    "Eirini Kasdagli <kasdaglie@ufl.edu>",
    "Jeffrey Andrews <jeffrey.andrews@ufl.edu>"]

from copy import copy
import numpy as np
import astropy.constants as con
import astropy.units as unt
import pandas as pd
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
               #'S1_metallicity'
               #'Z/Zo'
               #'Z/Zo',
               #'log_g',
               #'Teff'
                ]


def find_max_time(history):
    """Find the max time of a population."""
    times = history[history.event == "maxtime"].time

    return np.max(times)

def load_posydon_population(population_file):
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

    history = pd.read_hdf(population_file, key='history')
    max_time = find_max_time(history)
    final_stars = history[(history.time == max_time) & (history.event == "END")]
    zams_stars = history[history.event == 'ZAMS']
    for col in final_stars:
        if col not in keys_to_save:
            del final_stars[str(col)]
    #Find the stripped stars and put their final mass as M_init (mass at ZAMS).
    for star in ['S1','S2']:
        stripped_index = copy(final_stars[final_stars[f'{star}_state'].isin(stripped_stars)].index.values)
        M_init = zams_stars.loc[stripped_index][f'{star}_mass'].values
        final_stars.loc[stripped_index][f'{star}_mass'] = M_init
    pop = copy(final_stars.reset_index())
    #Get the metallicity values from ZAMS
    pop['Z/Zo'] = np.ones(len(pop))*zams_stars['S1_metallicity'].iloc[0]/Zo
    #Add two extra columns in the pop data
    for star in ['S1','S2']:
        M = np.asarray(pop[f'{star}_mass'])*con.M_sun
        R = np.asarray(10**pop[f'{star}_log_R'])*con.R_sun
        L = np.asarray(10**pop[f'{star}_log_L'])*con.L_sun
        pop[f'{star}_log_g'] = np.log10(con.G*M/R**2/(unt.cm/unt.s**2))
        pop[f'{star}_Teff'] = (L/(4*np.pi*R**2*con.sigma_sb))**0.25/unt.K
    return pop


def calculate_logg(row,x):
    """Calculate log_g of S1 and S2 from the population data"""
    mass = row[f'{x}_mass']*con.M_sun
    R = 10**row[f'{x}_log_R']*con.R_sun
    return np.log10(con.G*mass/R**2/(unt.cm/unt.s**2))

def calculate_Teff(row,x):
    """Calculate Teff of S1 and S2 from the population data"""
    L = 10**row[f'{x}_log_L']*con.L_sun
    R = 10**row[f'{x}_log_R']*con.R_sun
    return (L/(4*np.pi*R**2*con.sigma_sb))**0.25/unt.K

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
