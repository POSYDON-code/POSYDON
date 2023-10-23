"""Handle I/O operations for the population synthesis code."""

from configparser import ConfigParser
import ast
import importlib
import os
import errno
import pprint

import numpy as np
import pandas as pd

from posydon.binary_evol.simulationproperties import SimulationProperties


"""
POSYDON Data Types are enforced when converting BinaryStar
and SingleStar instances to Pandas DataFrames. This is done to
ensure memory efficient data storage and to solve problems
combining temp batch files.

Check Pandas docs for allowed data types in DataFrames.
"""

BINARYPROPERTIES_DTYPES = {
    # The state and event of the system. For more information, see
    # `posydon.utils.common_functions.get_binary_state_and_event_and_mt_case()
    'state' : 'string',                    #
    'event': 'string',
    'time': 'float64',                     # age of the system (yr)
    'separation': 'float64',               # binary orbital separation (solar radii)
    'orbital_period': 'float64',           # binary orbital period (days)
    'eccentricity': 'float64',             # binary eccentricity
    'V_sys': 'object',                    # list of the 3 systemic velocity coordinates
    # (R_{star} - R_{Roche_lobe}) / R_{Roche_lobe}...
    'rl_relative_overflow_1': 'float64',   # ...for star 1
    'rl_relative_overflow_2': 'float64',   # ...for star 2
    'lg_mtransfer_rate': 'float64',        # log10 of mass lost from the donor (Msun/yr)
    'mass_transfer_case': 'string',       # current mass transfer case of the system.
                                    # See `get_binary_state_and_event_and_mt_case`
                                    # in `posydon.utils.common_functions`.
    'trap_radius': 'float64',
    'acc_radius': 'float64',
    't_sync_rad_1': 'float64',
    't_sync_conv_1': 'float64',
    't_sync_rad_2': 'float64',
    't_sync_conv_2': 'float64',
    'nearest_neighbour_distance': 'object',   # the distance of system from its nearest
                                            # neighbour of MESA binary system  in case
                                            # of interpolation during the the end of
                                            # the previous step including MESA psygrid.
                                            # The distance is normalized in the
                                            # parameter space and limits at which it
                                            # was calculated. See `mesa_step` for more.
}


STARPROPERTIES_DTYPES = {
    'state': 'string',            # the evolutionary state of the star. For more info see
                                # `posydon.utils.common_functions.check_state_of_star`
    'metallicity': 'float64',      # initial mass fraction of metals
    'mass': 'float64',             # mass (solar units)
    'log_R': 'float64',            # log10 of radius (solar units)
    'log_L': 'float64',            # log10 luminosity (solar units)
    'lg_mdot': 'float64',          # log10 of the absolulte mass change of the star due to
                        # winds and RLO coming from binary history (Msun/yr)
    'lg_system_mdot': 'float64',   # log10 of the system's absolute mass loss from around
                        # the star due to inneficient mass transfer (Msun/yr)
    'lg_wind_mdot': 'float64',     # log10 of the absolute wind mass loss of the star from
                        # `lg_wind_mdot_1/2` in binary_history (Msun/yr)
    'he_core_mass': 'float64',     # in solar units
    'he_core_radius': 'float64',   # in solar units
    'c_core_mass': 'float64',      # in solar units
    'c_core_radius': 'float64',    # in solar units
    'o_core_mass': 'float64',      # in solar units
    'o_core_radius': 'float64',    # in solar units
    'co_core_mass': 'float64',     # in solar units
    'co_core_radius': 'float64',   # in solar units
    # Central mass fractions
    'center_h1': 'float64',
    'center_he4': 'float64',
    'center_c12': 'float64',
    'center_n14': 'float64',
    'center_o16': 'float64',
    'surface_h1': 'float64',
    # Mass fractions at the surface
    'surface_he4': 'float64',
    'surface_c12': 'float64',
    'surface_n14': 'float64',
    'surface_o16': 'float64',
    # Energy production from nuclear burning (solar units)
    'log_LH': 'float64',           # H burning
    'log_LHe': 'float64',          # He burning
    'log_LZ': 'float64',           # non-H and non-He burning
    'log_Lnuc': 'float64',         # total nuclear burning energy production
    'c12_c12': 'float64',          # Carbon burning through c12-c12
    'center_gamma': 'float64',
    'avg_c_in_c_core': 'float64',  # average carbon mass fraction at the carbon core
    'surf_avg_omega': 'float64',   # surface average rotational angular velocity (rad/sec)
    'surf_avg_omega_div_omega_crit': 'float64',    # ratio of `surf_avg_omega` to the
                                        # surface critical rotation limit,
                                        # taking into account centrifugal force
                                        # & radiation pressure (Eddington lim.)
    'total_moment_of_inertia': 'float64',      # total moment of inertia (gr*cm^2)
    'log_total_angular_momentum': 'float64',   # log10 of total ang. momentum (gr*cm^2/s)
    'spin': 'float64',                         # the dimesionless spin of the star, if it
                                    # is a compact object, which is equal to
                                    # c*J/(GM^2).
    'conv_env_top_mass': 'float64',
    'conv_env_bot_mass': 'float64',
    'conv_env_top_radius': 'float64',
    'conv_env_bot_radius': 'float64',
    'conv_env_turnover_time_g': 'float64',
    'conv_env_turnover_time_l_b': 'float64',
    'conv_env_turnover_time_l_t': 'float64',
    'envelope_binding_energy': 'float64',
    'mass_conv_reg_fortides': 'float64',
    'thickness_conv_reg_fortides': 'float64',
    'radius_conv_reg_fortides': 'float64',
    'lambda_CE_1cent': 'float64',
    'lambda_CE_10cent': 'float64',
    'lambda_CE_30cent': 'float64',
    'lambda_CE_pure_He_star_10cent': 'float64',
    'profile': 'object',    # the profile of the star, including extended information of
                            # its internal structure, for a specific timestep, usually for
                            # the end of the previous step including MESA psygrid.
}

EXTRA_BINARY_COLUMNS_DTYPES = {
    'step_names' : 'string',
    'step_times' : 'float64', 
}

# no default extras for history attributes
EXTRA_STAR_COLUMNS_DTYPES = {}

SCALAR_NAMES_DTYPES = {
    'natal_kick_array_0': 'float64',
    'natal_kick_array_1': 'float64',
    'natal_kick_array_2': 'float64',
    'natal_kick_array_3': 'float64',
    'SN_type':'string',
    'f_fb': 'float64',
    'sping_orbit_tilt': 'float64',
}


def clean_binary_history_df(binary_df, extra_binary_dtypes_user=None, 
                            extra_S1_dtypes_user=None, extra_S2_dtypes_user=None):
    """Take a posydon binary history DataFrame from the
    BinaryStar.to_df method and clean the data for saving
    by setting Data Types of the columns explicitly.
    
    Parameters
    ---------
    binary_df : DataFrame
        A pandas Dataframe containing binary history
    extra_binary_dtypes_user : dict, optional
        A dictionary with extra column names as keys, and their
        associated data types as values.
    extra_S1_dtypes_user : dict, optional
        Same as above, but only for star 1.
    extra_S2_dtypes_user : dict, optional
        Same as above, but only for star 2.
        
    Returns
    -------
    binary_df : DataFrame
        A cleaned binary history ready for saving to HDF.
    """
    
    assert isinstance( binary_df, pd.DataFrame )
    
    # User specified extra binary and star columns
    if extra_binary_dtypes_user is None:
        extra_binary_dtypes_user = {}
    if extra_S1_dtypes_user is None:
        extra_S1_dtypes_user = {}
    if extra_S2_dtypes_user is None:
        extra_S2_dtypes_user = {}

    # try to coerce data types automatically first
    binary_df = binary_df.infer_objects()

    # Set data types for all columns explicitly
    hist_columns =  set(binary_df.columns)
    
    # combine columns for binary history with extras
    BP_comb_extras_dict = {**BINARYPROPERTIES_DTYPES, **EXTRA_BINARY_COLUMNS_DTYPES, 
                           **extra_binary_dtypes_user}
    # combine columns for star 1 and star 2 history with extras
    SP_comb_S1_dict = {**STARPROPERTIES_DTYPES, **EXTRA_STAR_COLUMNS_DTYPES, 
                       **extra_S1_dtypes_user}
    SP_comb_S2_dict = {**STARPROPERTIES_DTYPES, **EXTRA_STAR_COLUMNS_DTYPES, 
                       **extra_S2_dtypes_user}
    
    # All default & user passed keys which we have mappings for
    binary_keys = set( BP_comb_extras_dict.keys() )
    S1_keys = set( ['S1_' + key for key in SP_comb_S1_dict.keys()] )
    S2_keys = set( ['S2_' + key for key in SP_comb_S2_dict.keys()] )
    
    # Find common keys between the binary_df and our column-dtype mapping
    common_keys =  hist_columns & ( binary_keys | S1_keys | S2_keys )

    # Create a dict with column-dtype mapping only for columns in binary_df
    common_dtype_dict = {}
    for key in common_keys:
        if key in binary_keys:
            common_dtype_dict[key] = BP_comb_extras_dict.get( key )
        elif key in S1_keys:
            common_dtype_dict[key] = SP_comb_S1_dict.get( key.replace('S1_', '') )
        elif key in S2_keys:
            common_dtype_dict[key] = SP_comb_S2_dict.get( key.replace('S2_', '') )
        else:
            raise ValueError(f'No data type found for {key}. Dtypes must be explicity declared.')
    # set dtypes
    binary_df = binary_df.astype( common_dtype_dict )
    # unset clean str data because pandas strings are broken for hdf saving
    convert_to_obj = {}
    for key, val in common_dtype_dict.items():
        if val == 'string':
            convert_to_obj[key] = 'object'
    binary_df = binary_df.astype( convert_to_obj )
    return binary_df


def clean_binary_oneline_df(oneline_df, extra_binary_dtypes_user=None, 
                            extra_S1_dtypes_user=None, extra_S2_dtypes_user=None):
    """Take a posydon binary oneline DataFrame from the
    BinaryStar.to_oneline_df method and clean the data for saving
    by setting Data Types of the columns explicitly.
    
    This method is similar to clean_binary_history_df since they
    have many overalapping columns, with a few extras and different naming.
    
    Note: there may be edge cases not handed if new scalar_names are added.
    
    Parameters
    ---------
    binary_df : DataFrame
        A pandas Dataframe containing binary history
    extra_binary_dtypes_user : dict, optional
        A dictionary with extra column names as keys, and their
        associated data types as values.
    extra_S1_dtypes_user : dict, optional
        Same as above, but only for star 1.
    extra_S2_dtypes_user : dict, optional
        Same as above, but only for star 2.
        
    Returns
    -------
    binary_df : DataFrame
        A cleaned binary history ready for saving to HDF.
    """
    assert isinstance( oneline_df, pd.DataFrame )
    
    # User specified extra binary and star columns
    if extra_binary_dtypes_user is None:
        extra_binary_dtypes_user = {}
    if extra_S1_dtypes_user is None:
        extra_S1_dtypes_user = {}
    if extra_S2_dtypes_user is None:
        extra_S2_dtypes_user = {}

    # try to coerce data types automatically first
    oneline_df = oneline_df.infer_objects()

    # Set data types for all columns explicitly
    oneline_columns =  set(oneline_df.columns)
    
    # combine columns for binary history with extras
    BP_comb_extras_dict = {**BINARYPROPERTIES_DTYPES, **EXTRA_BINARY_COLUMNS_DTYPES, 
                           **extra_binary_dtypes_user}
    # combine columns for star 1 and star 2 history with extras
    SP_comb_S1_dict = {**STARPROPERTIES_DTYPES, **EXTRA_STAR_COLUMNS_DTYPES, 
                       **SCALAR_NAMES_DTYPES, **extra_S1_dtypes_user}
    SP_comb_S2_dict = {**STARPROPERTIES_DTYPES, **EXTRA_STAR_COLUMNS_DTYPES, 
                       **SCALAR_NAMES_DTYPES, **extra_S2_dtypes_user}
    
    # All default & user passed keys which we have mappings for
    binary_keys = set( [key + '_i' for key in BP_comb_extras_dict.keys()] 
                     + [key + '_f' for key in BP_comb_extras_dict.keys()] )
    S1_keys = set( ['S1_' + key + '_i' for key in SP_comb_S1_dict.keys()] 
                 + ['S1_' + key + '_f' for key in SP_comb_S1_dict.keys()] )
    S2_keys = set( ['S2_' + key + '_i' for key in SP_comb_S2_dict.keys()]
                 + ['S2_' + key + '_f' for key in SP_comb_S2_dict.keys()] )
    
    # Find common keys between the binary_df and our column-dtype mapping
    common_keys =  oneline_columns & ( binary_keys | S1_keys | S2_keys )
    

    # Create a dict with column-dtype mapping only for columns in binary_df
    common_dtype_dict = {}
    # helper function to remove the '_i' and '_f' 
    strip_prefix_and_suffix = lambda key : key.strip('_i').strip('_f').strip('S1_').strip('S2_')
    for key in common_keys:
        if key in binary_keys:
            common_dtype_dict[key] = BP_comb_extras_dict.get( strip_prefix_and_suffix(key) )
        elif key in S1_keys:
            common_dtype_dict[key] = SP_comb_S1_dict.get( strip_prefix_and_suffix(key) )
        elif key in S2_keys:
            common_dtype_dict[key] = SP_comb_S2_dict.get( strip_prefix_and_suffix(key) )
        else:
            raise ValueError(f'No data type found for {key}. Dtypes must be explicity declared.')
    # set dtypes
    oneline_df = oneline_df.astype( common_dtype_dict )
    # unset clean str data because pandas strings are broken for hdf saving
    convert_to_obj = {}
    for key, val in common_dtype_dict.items():
        if val == 'string':
            convert_to_obj[key] = 'object'
    oneline_df = oneline_df.astype( convert_to_obj )
    
    return oneline_df



def parse_inifile(path, verbose=False):
    """Parse an inifile for evolving binary populations.

    Parameters
    ---------
    path : str or list like
        Path to inifile. If multiple files are given,
        duplicate args are overwritten (stacked) first
        to last.
    verbose : bool
        Print helpful info.

    Returns
    -------
    parser : <class, ConfigParser>
        An instance of ConfigParser.

    """
    parser = ConfigParser()
    # Keys do not become lowercase by default
    parser.optionxform = lambda option: option

    if isinstance(path, str):
        path = os.path.abspath(path)
        if verbose:
            print('Reading inifile: \n\t{}'.format(path))
        if not os.path.exists(path):
            raise FileNotFoundError(
                  errno.ENOENT, os.strerror(errno.ENOENT), path)
    elif isinstance(path, (list, np.ndarray)):
        path = [os.path.abspath(f) for f in path]

        if verbose:
            print('Reading inifiles: \n{}'.format(pprint.pformat(path)))
        bad_files = []
        for f in path:
            if not os.path.exists(f):
                bad_files.append(f)
        if bool(bad_files):
            raise FileNotFoundError(
                  errno.ENOENT, os.strerror(errno.ENOENT), bad_files)

    else:
        raise ValueError("Path must be a string or list of strings. Given {}".
                         format(path))

    files_read = parser.read(path)
    # Catch silent errors from configparser.read
    if len(files_read) == 0:
        raise ValueError("No files were read successfully. Given {}.".
                         format(path))
    return parser


def simprop_kwargs_from_ini(path, verbose=False):
    """Convert an inifile into kwargs for the SimulationProperties class.

    Parameters
    ---------
    path : str or list like
        Path to inifile. If multiple files are given,
        duplicate args are overwritten (stacked) first
        to last.

    verbose : bool
        Print helpful info.

    Returns
    -------
    parser_dict : <class, dict>
        The inifile converted to the kwargs.
    """
    parser = parse_inifile(path, verbose=verbose)
    parser_dict = {}
    for section in parser:
        # skip default section
        if section == 'DEFAULT':
            continue

        # evaluate str values as literal python and put
        # into dict because parser only handles strings
        sect_dict = {key: ast.literal_eval(val)
                     for key, val in parser[section].items()}

        # only try imports from sections with 'import'
        if 'import' in list(sect_dict.keys()):
            # from posydon.... import ....
            #      ^ import      as   ^ name
            import_and_name = sect_dict.pop('import')
            package = sect_dict.pop('absolute_import', None)

            # import module
            module = importlib.import_module(import_and_name[0],
                                             package=package)
            # extract class or function
            cls = getattr(module, import_and_name[1])
            # match the form SimulationProperties expects
            parser_dict[section] = (cls, sect_dict)

        # Try to find user defined hooks using absolute or
        # relative imports. Requires slighltly different syntax.
        if section == 'extra_hooks':
            hooks_list = []
            names = [name for name in list(sect_dict.keys())
                     if 'extra' not in name]
            # do not take extra_pre_step options etc.

            # loop to try and find all hooks classes
            for h in range(1, len(names)):
                valid_names = [name for name in names
                               if int(name.split('_')[-1]) == h]

                # if there are no options with the index, continue
                if not valid_names:
                    continue

                import_and_name = sect_dict.pop('import_{}'.format(h))
                package = sect_dict.pop('absolute_import_{}'.format(h), None)
                class_kwargs = sect_dict.pop('kwargs_{}'.format(h), dict())

                module = importlib.import_module(import_and_name[0],
                                                 package=package)
                cls = getattr(module, import_and_name[1])

                hooks_list.append((cls, class_kwargs))

            parser_dict[section] = hooks_list

    return parser_dict


def binarypop_kwargs_from_ini(path, verbose=False):
    """Convert an inifile into kwargs for the BinaryPopulation class.

    Parameters
    ---------
    path : str or list like
        Path to inifile. If multiple files are given,
        duplicate args are overwritten (stacked) first
        to last.

    verbose : bool
        Print helpful info.

    Returns
    -------
    parser_dict : <class, dict>
        The inifile converted to the kwargs.
    """
    parser = parse_inifile(path, verbose=verbose)

    # make sure the inifile has all of these sections
    for sect_name in ['BinaryPopulation_options', 'BinaryStar_output',
                      'SingleStar_1_output', 'SingleStar_2_output']:
        assert sect_name in list(parser.sections())

    pop_kwargs = dict()
    for section in parser.sections():
        if section == 'BinaryPopulation_options':
            for key, val in parser[section].items():
                pop_kwargs[key] = ast.literal_eval(val)

            # hard code for running with MPI, only try to import if use_MPI
            if pop_kwargs['use_MPI']:
                from mpi4py import MPI
                pop_kwargs['comm'] = MPI.COMM_WORLD
            else:
                pop_kwargs['comm'] = None

        # right now binary, S1, and S2 output kwargs are all passed
        # into the BinaryPopulation during init
        elif section == 'BinaryStar_output':
            binary_kwargs = dict()
            for key, val in parser[section].items():
                binary_kwargs[key] = ast.literal_eval(val)
            pop_kwargs = {**pop_kwargs, **binary_kwargs}

        elif section == 'SingleStar_1_output':
            S1_kwargs = dict()
            for key, val in parser[section].items():
                S1_kwargs[key] = ast.literal_eval(val)
            pop_kwargs['include_S1'] = S1_kwargs.pop('include_S1')
            if pop_kwargs['include_S1']:
                pop_kwargs['S1_kwargs'] = S1_kwargs

        elif section == 'SingleStar_2_output':
            S2_kwargs = dict()
            for key, val in parser[section].items():
                S2_kwargs[key] = ast.literal_eval(val)
            pop_kwargs['include_S2'] = S2_kwargs.pop('include_S2')
            if pop_kwargs['include_S2']:
                pop_kwargs['S2_kwargs'] = S2_kwargs

    # finally get the population properties
    sim_prop_kwargs = simprop_kwargs_from_ini(path)
    pop_kwargs['population_properties'] = SimulationProperties(
        **sim_prop_kwargs)

    return pop_kwargs
