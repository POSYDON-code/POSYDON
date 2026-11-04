"""Read default values from the packaged population inifile.

The default population inifile, ``posydon/popsyn/population_params_default.ini``,
is the single source of default values for the evolution steps and the binary
population. The classes take their ``DEFAULT_KWARGS`` from the sections of
this file, so building a class directly gives the same defaults as building
it from an inifile.

This module only depends on the standard library, so it can be imported by
any POSYDON module without creating circular imports.
"""


__authors__ = [
    "Max Briel <max.briel@gmail.com>",
]


import ast
import copy
import os
from configparser import ConfigParser
from functools import lru_cache

DEFAULT_INI_PATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    'popsyn', 'population_params_default.ini')

# keys of a section, which select the class and are not class parameters
IMPORT_KEYS = ('import', 'absolute_import')


@lru_cache(maxsize=None)
def _read_default_sections():
    """Read the default inifile once and evaluate all its values.

    Returns
    -------
    dict
        Maps each section name to a dict of its evaluated values.
    """
    parser = ConfigParser()
    # keys keep their case
    parser.optionxform = lambda option: option
    if not parser.read(DEFAULT_INI_PATH):
        raise FileNotFoundError(f"Default inifile not found: {DEFAULT_INI_PATH}")
    return {section: {key: ast.literal_eval(value)
                      for key, value in parser[section].items()}
            for section in parser.sections()}


def section_defaults(section, exclude=()):
    """Get the default values of one section of the default inifile.

    Parameters
    ----------
    section : str
        Name of the section, e.g. 'step_SN' or 'BinaryPopulation_options'.
    exclude : iterable of str
        Keys to leave out, e.g. the ones handled by another class. The keys
        'import' and 'absolute_import' are always left out.

    Returns
    -------
    dict
        A new dict with the default values; mutable values are copies, so it
        can be changed without affecting other users of the defaults.
    """
    sections = _read_default_sections()
    if section not in sections:
        raise KeyError(f"Section [{section}] is not in the default inifile "
                       f"{DEFAULT_INI_PATH}.")
    skip = set(IMPORT_KEYS) | set(exclude)
    return copy.deepcopy({key: value for key, value in sections[section].items()
                          if key not in skip})
