"""Input/output functions for spectral synthesis parameters.
Mirrors posydon/popsyn/io.py [3].
"""

__authors__ = [
    "Eirini Kasdagli <kasdaglie@ufl.edu>",
]

import configparser
import os


# mirrors: saved_ini_parameters in binarypopulation.py [2]
SAVED_SPECTRAL_PARAMETERS = [
    'metallicity',
    'save_data',
    'population_file',
    'spectral_type',
    'isochromes',
    'main_grid',
    'secondary_grid',
    'lam_min',
    'lam_max',
    'lam_res',
    'ostar_temp_cut_off',
    'bstar_temp_cut_off',
]

# mirrors: type conversion in binarypop_kwargs_from_ini [3]
SPECTRAL_PARAM_TYPES = {
    'metallicity'        : float,
    'save_data'          : bool,
    'population_file'    : str,
    'spectral_type'      : bool,
    'isochromes'         : bool,
    'main_grid'          : str,
    'secondary_grid'     : str,
    'lam_min'            : float,
    'lam_max'            : float,
    'lam_res'            : int,
    'ostar_temp_cut_off' : float,
    'bstar_temp_cut_off' : float,
}

# mirrors: default_kwargs in posydon/popsyn/defaults.py
DEFAULT_SPECTRAL_KWARGS = {
    'metallicity'        : 1.0,
    'save_data'          : True,
    'population_file'    : None,
    'spectral_type'      : False,
    'isochromes'         : False,
    'main_grid'          : None,
    'secondary_grid'     : None,
    'lam_min'            : 110.0,
    'lam_max'            : 100000.0,
    'lam_res'            : 300000,
    'ostar_temp_cut_off' : 28000.0,
    'bstar_temp_cut_off' : 15000.0,
}


def spectral_kwargs_from_ini(ini_file):
    """Read spectral synthesis parameters from an ini file.
    Mirrors binarypop_kwargs_from_ini in posydon/popsyn/io.py [3].

    Parameters
    ----------
    ini_file : str
        Path to the ini file.

    Returns
    -------
    dict
        Dictionary of spectral synthesis parameters.

    Raises
    ------
    FileNotFoundError
        If the ini file does not exist.
    ValueError
        If a required parameter is missing or has an invalid type.
    """
    if not os.path.exists(ini_file):
        raise FileNotFoundError(f'ini file not found: {ini_file}')

    config = configparser.ConfigParser()
    config.read(ini_file)

    # mirrors: default_kwargs.copy() in BinaryPopulation.__init__ [2]
    kwargs = DEFAULT_SPECTRAL_KWARGS.copy()

    if 'spectral_synthesis' not in config:
        raise ValueError(
            f'{ini_file} does not contain a [spectral_synthesis] section!'
        )

    section = config['spectral_synthesis']

    # mirrors: type conversion loop in binarypop_kwargs_from_ini [3]
    for key, val in section.items():
        if key not in SPECTRAL_PARAM_TYPES:
            continue
        param_type = SPECTRAL_PARAM_TYPES[key]
        try:
            # mirrors: bool handling in binarypop_kwargs_from_ini [3]
            if param_type == bool:
                kwargs[key] = section.getboolean(key)
            else:
                kwargs[key] = param_type(val)
        except (ValueError, TypeError) as e:
            raise ValueError(
                f'Invalid value for {key}={val} in {ini_file}: {e}'
            )

    return kwargs


def validate_spectral_ini(ini_file):
    """Validate the spectral synthesis ini file.
    Mirrors validate_ini_file in posydon/CLI/popsyn/setup.py [4].

    Parameters
    ----------
    ini_file : str
        Path to the ini file.

    Raises
    ------
    FileNotFoundError
        If the ini file does not exist.
    ValueError
        If required parameters are missing.
    """
    if not os.path.exists(ini_file):
        raise FileNotFoundError(f'ini file not found: {ini_file}')

    kwargs = spectral_kwargs_from_ini(ini_file)

    # mirrors: required parameter check in validate_ini_file [4]
    if kwargs.get('population_file') is None:
        raise ValueError(
            'population_file is required in the ini file!'
        )