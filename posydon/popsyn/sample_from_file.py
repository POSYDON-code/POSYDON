"""Get the initial parameters for a binary population."""


__authors__ = [
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>",
]


import os
import numpy as np
import pandas as pd
import warnings
from posydon.popsyn.independent_sample import (generate_orbital_periods,
                                               generate_orbital_separations,
                                               generate_eccentricities,
                                               generate_primary_masses,
                                               generate_secondary_masses)

PRIMARY_MASS_NAMES = ['s1_mass', 'primary_mass', 'mass_1', 'm_1', 'm1']
SECONDARY_MASS_NAMES = ['s2_mass', 'secondary_mass', 'mass_2', 'm_2', 'm2']
PERIOD_NAMES = ['orbital_period', 'period', 'p_orb', 'porb', 'p']
SEPARATION_NAMES = ['orbital_separation', 'separation', 'semi-major_axis',\
                    'semi_major_axis', 'a']
ECCENTRICITY_NAMES = ['eccentricity', 'ecc', 'e']
PRIMARY_KICK_VELOCITY_NAMES = ['s1_natal_kick_array_0', 'w_1', 'w1']
SECONDARY_KICK_VELOCITY_NAMES = ['s2_natal_kick_array_0', 'w_2', 'w2']
PRIMARY_KICK_AZIMUTHAL_ANGLE_NAMES = ['s1_natal_kick_array_1', 'phi_1', 'phi1']
SECONDARY_KICK_AZIMUTHAL_ANGLE_NAMES = ['s2_natal_kick_array_1', 'phi_2',\
                                        'phi2']
PRIMARY_KICK_POLAR_ANGLE_NAMES = ['s1_natal_kick_array_2', 'theta_1', 'theta1']
SECONDARY_KICK_POLAR_ANGLE_NAMES = ['s2_natal_kick_array_2', 'theta_2',\
                                    'theta2']
PRIMARY_KICK_MEAN_ANOMALY_NAMES = ['s1_natal_kick_array_3', 'mean_anomaly_1',\
                                   'mean_anomaly1']
SECONDARY_KICK_MEAN_ANOMALY_NAMES = ['s2_natal_kick_array_3',\
                                     'mean_anomaly_2', 'mean_anomaly2']

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
        raise FileNotFoundError(f'{filename} not found!')

    # Get number of data frame entries
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

    # Get primary masses
    key = infer_key(available_keys=df.keys(), allowed_keys=PRIMARY_MASS_NAMES)
    if key=='':
        warnings.warn(f'No primary mass column found in {filename}, hence get'
                      ' independent random ones.')
        m1_set = generate_primary_masses(number_of_binaries=set_n, **kwargs)
    else:
        m1_set = np.array(df[key])

    # Get secondary masses
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
        # Get orbital separations
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
        # Get orbital periods
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
        # Expand the sets by doubling them as long as it is needed
        while len(orbital_scheme_set)<index+number_of_binaries:
            orbital_scheme_set = np.append(orbital_scheme_set,\
                                           orbital_scheme_set, axis=0)
            eccentricity_set = np.append(eccentricity_set, eccentricity_set,\
                                         axis=0)
            m1_set = np.append(m1_set, m1_set, axis=0)
            m2_set = np.append(m2_set, m2_set, axis=0)
        # Only return the requested set size
        return orbital_scheme_set[index:index+number_of_binaries],\
               eccentricity_set[index:index+number_of_binaries],\
               m1_set[index:index+number_of_binaries],\
               m2_set[index:index+number_of_binaries]
    else:
        return orbital_scheme_set, eccentricity_set, m1_set, m2_set

def get_kick_samples_from_file(**kwargs):
    """Read a kicks for population of binaries from a file.

    Parameters
    ----------
    **kwargs : dictionary
        kwargs from BinaryPopulation class, which should contain
        read_samples_from_file

    Returns
    -------
    s1_natal_kick_array_set : ndarray of floats
        natal kick array for the primary star
        containing: kick velocity, azimuthal angle, polar angle, mean anomaly
    s2_natal_kick_array_set : ndarray of floats
        natal kick array for the secondary star
        containing: kick velocity, azimuthal angle, polar angle, mean anomaly

    """

    if 'read_samples_from_file' in kwargs:
        filename = kwargs['read_samples_from_file']
    else:
        raise KeyError("In get_kick_samples_from_file no "
                       "'read_samples_from_file' in kwargs.")
    # Check and read file into dataframe
    if os.path.isfile(filename):
        df = pd.read_csv(filename)
    else:
        raise FileNotFoundError(f'{filename} not found!')

    # Get number of data frame entries
    set_n = len(df)
    
    ## Get primary kick
    # Velocity
    key = infer_key(available_keys=df.keys(),\
                    allowed_keys=PRIMARY_KICK_VELOCITY_NAMES)
    if key=='':
        warnings.warn('No kick velocity column of primary found in '
                      f'{filename}, hence do not set them.')
        w_set = np.array(set_n*[None])
    else:
        w_set = np.array(df[key])
    # Azimuthal angle
    key = infer_key(available_keys=df.keys(),\
                    allowed_keys=PRIMARY_KICK_AZIMUTHAL_ANGLE_NAMES)
    if key=='':
        warnings.warn('No azimuthal angle column of primary found in '
                      f'{filename}, hence do not set them.')
        phi_set = np.array(set_n*[None])
    else:
        phi_set = np.array(df[key])
    # Polar angle
    key = infer_key(available_keys=df.keys(),\
                    allowed_keys=PRIMARY_KICK_POLAR_ANGLE_NAMES)
    if key=='':
        warnings.warn('No polar angle column of primary found in '
                      f'{filename}, hence do not set them.')
        theta_set = np.array(set_n*[None])
    else:
        theta_set = np.array(df[key])
    # Mean anomaly
    key = infer_key(available_keys=df.keys(),\
                    allowed_keys=PRIMARY_KICK_MEAN_ANOMALY_NAMES)
    if key=='':
        warnings.warn('No mean anomaly column of primary found in '
                      f'{filename}, hence do not set them.')
        mean_anomaly_set = np.array(set_n*[None])
    else:
        mean_anomaly_set = np.array(df[key])
    kick_1_set = np.transpose(np.array([w_set, phi_set, theta_set,\
                                        mean_anomaly_set]))

    ## Get secondary kick
    # Velocity
    key = infer_key(available_keys=df.keys(),\
                    allowed_keys=SECONDARY_KICK_VELOCITY_NAMES)
    if key=='':
        warnings.warn('No kick velocity column of secondary found in '
                      f'{filename}, hence do not set them.')
        w_set = np.array(set_n*[None])
    else:
        w_set = np.array(df[key])
    # Azimuthal angle
    key = infer_key(available_keys=df.keys(),\
                    allowed_keys=SECONDARY_KICK_AZIMUTHAL_ANGLE_NAMES)
    if key=='':
        warnings.warn('No azimuthal angle column of secondary found in '
                      f'{filename}, hence do not set them.')
        phi_set = np.array(set_n*[None])
    else:
        phi_set = np.array(df[key])
    # Polar angle
    key = infer_key(available_keys=df.keys(),\
                    allowed_keys=SECONDARY_KICK_POLAR_ANGLE_NAMES)
    if key=='':
        warnings.warn('No polar angle column of secondary found in '
                      f'{filename}, hence do not set them.')
        theta_set = np.array(set_n*[None])
    else:
        theta_set = np.array(df[key])
    # Mean anomaly
    key = infer_key(available_keys=df.keys(),\
                    allowed_keys=SECONDARY_KICK_MEAN_ANOMALY_NAMES)
    if key=='':
        warnings.warn('No mean anomaly column of secondary found in '
                      f'{filename}, hence do not set them.')
        mean_anomaly_set = np.array(set_n*[None])
    else:
        mean_anomaly_set = np.array(df[key])
    kick_2_set = np.transpose(np.array([w_set, phi_set, theta_set,\
                                        mean_anomaly_set]))

    if 'number_of_binaries' in kwargs:
        number_of_binaries = kwargs['number_of_binaries']
        if 'index' in kwargs:
            index = kwargs['index']
        else:
            index = 0
        # expand the sets by doubling them as long as it is needed
        while len(kick_1_set)<index+number_of_binaries:
            kick_1_set = np.append(kick_1_set, kick_1_set, axis=0)
            kick_2_set = np.append(kick_2_set, kick_2_set, axis=0)
        # only return the requested set size
        return kick_1_set[index:index+number_of_binaries],\
               kick_2_set[index:index+number_of_binaries]
    else:
        return kick_1_set, kick_2_set

