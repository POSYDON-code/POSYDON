"""Handle I/O operations for the population synthesis code."""

from configparser import ConfigParser
import ast
import importlib
import os
import errno
import pprint

import numpy as np

from posydon.binary_evol.simulationproperties import SimulationProperties


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
