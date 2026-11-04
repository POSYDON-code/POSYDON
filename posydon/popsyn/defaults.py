"""Default binary population parameters.

The defaults are the values of the [BinaryPopulation_options] section of the
default population inifile, ``posydon/popsyn/population_params_default.ini``,
which is the single source of POSYDON's default parameters.
"""


__authors__ = [
    "Kyle Akira Rocha <kylerocha2024@u.northwestern.edu>",
    "Jeffrey Andrews <jeffrey.andrews@northwestern.edu>",
    "Max Briel <max.briel@gmail.com>",
]


from posydon.utils.ini_defaults import section_defaults

default_kwargs = section_defaults('BinaryPopulation_options')
