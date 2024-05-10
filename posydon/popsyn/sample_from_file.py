"""Get the initial parameters for a binary population."""


__authors__ = [
    "Matthias Kruckow <matthias.kruckow@unige.ch>",
]


import os
import numpy as np
import pandas as pd
import warnings
from posydon.utils.posydonerror import FileError
from posydon.popsyn.independent_sample import (generate_orbital_periods,
                                               generate_orbital_separations,
                                               generate_eccentricities,
                                               generate_primary_masses,
                                               generate_secondary_masses)

PRIMARY_MASS_NAMES = ['primary_mass', 'mass_1', 'm_1', 'm1']
SECONDARY_MASS_NAMES = ['secondary_mass', 'mass_2', 'm_2', 'm2']
PERIOD_NAMES = ['orbital_period', 'period', 'p_orb', 'porb', 'p']
SEPARATION_NAMES = ['orbital_separation', 'separation', 'semi-major_axis',\
                    'semi_major_axis', 'a']
ECCENTRICITY_NAMES = ['eccentricity', 'ecc', 'e']

def infer_key(available_keys=[], allowed_keys=[]):
    """Infer key from list of allowed keys.

    Parameters
    ----------
    available_keys : iterable object of str, e.g. list of str
        Collection of available keys.
    allowed_keys : iterable object of str, e.g. list of str
        Collection of allowed keys.

    Returns
    -------
    key : str
        The first matched key.
    """

    for k in allowed_keys:
        for key in available_keys:
            if key.casefold() == k.casefold():
                return key

    return ''

def get_samples_from_file(orbital_scheme='', **kwargs):
    """Read a population of binaries at ZAMS from a file.

    Parameters
    ----------
    orbital_scheme : str
        Scheme to get the orbit: 'separation', 'period'
    **kwargs : dictionary
        kwargs from BinaryPopulation class, which should contain
        read_samples_from_file

    Returns
    -------
    orbital_scheme_set : ndarray of floats
        orbital separations/periods depending on the scheme
    eccentricity_set : ndarray of floats
        eccentricities
    m1_set : ndarray of floats
        primary masses
    m2_set : ndarray of floats
        secondary masses

    """

    if 'read_samples_from_file' in kwargs:
        filename = kwargs['read_samples_from_file']
    else:
        raise KeyError("In get_samples_from_file no 'read_samples_from_file'"
                       " in kwargs.")
    # Check and read file into dataframe
    if os.path.isfile(filename):
        df = pd.read_csv(filename)
    else:
        raise FileError(f'{filename} not found!')

    # get number of data frame entries
    set_n = len(df)
    
    # Get eccentricities
    key = infer_key(available_keys=df.keys(), allowed_keys=ECCENTRICITY_NAMES)
    if key=='':
        warnings.warn(f'No eccentricity column found in {filename}, hence get'
                      ' independent random ones.')
        eccentricity_set = generate_eccentricities(number_of_binaries=set_n,\
                                                   **kwargs)
    else:
        eccentricity_set = np.array(df[key])

    # Generate primary masses
    key = infer_key(available_keys=df.keys(), allowed_keys=PRIMARY_MASS_NAMES)
    if key=='':
        warnings.warn(f'No primary mass column found in {filename}, hence get'
                      ' independent random ones.')
        m1_set = generate_primary_masses(number_of_binaries=set_n, **kwargs)
    else:
        m1_set = np.array(df[key])

    # Generate secondary masses
    key = infer_key(available_keys=df.keys(),
                    allowed_keys=SECONDARY_MASS_NAMES)
    if key=='':
        warnings.warn(f'No secondary mass column found in {filename}, hence'
                      ' get independent random ones.')
        m2_set = generate_secondary_masses(m1_set, number_of_binaries=set_n,\
                                           **kwargs)
    else:
        m2_set = np.array(df[key])

    if orbital_scheme == 'separation':
        # Generate orbital separations
        key = infer_key(available_keys=df.keys(),
                        allowed_keys=SEPARATION_NAMES)
        if key=='':
            warnings.warn(f'No separation column found in {filename}, hence'
                          ' get independent random ones.')
            orbital_scheme_set = generate_orbital_separations(\
                                  number_of_binaries=set_n, **kwargs)
        else:
            orbital_scheme_set = np.array(df[key])
    elif orbital_scheme == 'period':
        # Generate orbital periods
        key = infer_key(available_keys=df.keys(), allowed_keys=PERIOD_NAMES)
        if key=='':
            warnings.warn(f'No period column found in {filename}, hence get'
                          ' independent random ones.')
            orbital_scheme_set = generate_orbital_periods(m1_set,\
                                  number_of_binaries=set_n, **kwargs)
        else:
            orbital_scheme_set = np.array(df[key])
    else:
        raise ValueError("Allowed orbital schemes are separation or period.")

    if 'number_of_binaries' in kwargs:
        number_of_binaries = kwargs['number_of_binaries']
        if 'index' in kwargs:
            index = kwargs['index']
        else:
            index = 0
        # expand the sets by doubling them as long as it is needed
        while len(orbital_scheme_set)<index+number_of_binaries:
            orbital_scheme_set = np.append(orbital_scheme_set,\
                                           orbital_scheme_set)
            eccentricity_set = np.append(eccentricity_set, eccentricity_set)
            m1_set = np.append(m1_set, m1_set)
            m2_set = np.append(m2_set, m2_set)
        # only return the requested set size
        return orbital_scheme_set[index:index+number_of_binaries],\
               eccentricity_set[index:index+number_of_binaries],\
               m1_set[index:index+number_of_binaries],\
               m2_set[index:index+number_of_binaries]
    else:
        return orbital_scheme_set, eccentricity_set, m1_set, m2_set

