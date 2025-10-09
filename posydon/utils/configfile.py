"""Module for providing support for configuration files.

The ConfigClass loads, edits and saves configuration files in JSON format. The
data are encapsulated in a single dictionary mapping the variable names (keys)
to their values.

ConfigFile is an alternative to the `configparse` module in Python's standard
libray. While `configparser` can handle only string values and requires using
sections (as in .ini files), `ConfigFile` is simpler and faster to use, while
permits keys and values of any type that a Python dictionary permits.

As it relies on the Python `json` module, an internal function is called when
a values is not known to `json`, e.g. numpy arrays.

Finally, the function `parse_inifile` is defined to help reading .ini files.

Examples
--------
(1) Saving a dictionary:

    D = {"username": "konstantinos",
         "nodes": [1, 5, 3],
         "memory": 1024,
         "code": 42,
         "std.out": "/data/output/std.out"}

    config_file = ConfigFile("tmp.cfg")
    config_file.update(D)
    config_file.save()

    OR

    my_config_file = ConfigFile()
    my_config_file.update(D)
    my_config_file.save(filename)


(2) Loading and printing configuration:

    config = ConfigFile(filename)
    print("Loaded a config file with {} entries.".format(len(config)))
    print("The code is", config["code"])
    print(config)

(4) Altering entries (and creating if non-existing):

    config["machine"] = "supercomputer"
    config["code"] = 24
    config["H0 and error"] = (67.8, 0.9)

(3) Loading configuration from multiple files:

    config = CongfigFile("config1.json")
    config.load("config2.json")

    OR

    config.load("config1.json")
    config.load("config2.json")

    If the two loaded configuration files have common keys, then an Exception
    will occur. To allow updates, e.g. in the case of default configuration
    and user configuration overriding the former, then:

    config = ConfigFile("default.cfg")
    config.load("user.cfg", can_update=True)

(5) Iterating entries:

    config = ConfigFile("my_conf.json")

    print("All configuration items:")
    for key in config:
        print("    {}: {}".format(key, config[key]))

    print("Ok, I'll repeat that...:")
    for key, value in config.items():
        print("    {}: {}".format(key, value))

"""


__authors__ = [
    "Konstantinos Kovlakas <Konstantinos.Kovlakas@unige.ch>",
    "Scott Coughlin <scottcoughlin2014@u.northwestern.edu>",
    "Jeffrey Andrews <jeffrey.andrews@northwestern.edu>",
]


import ast
import configparser
import copy
import json
import operator
import os

import numpy as np


class ConfigFile:
    """Class handling input, process and output of configurations."""

    def __init__(self, path=None):
        """Initialize a ConfigFile with or without a path."""
        self.entries = {}
        self.path = path
        if self.path is not None:
            if os.path.exists(self.path):
                self.load(path)

    def deepcopy(self):
        """Make a deep copy of the object."""
        newobj = ConfigFile()
        newobj.path = copy.deepcopy(self.path)
        newobj.entries = copy.deepcopy(self.entries)
        return newobj

    @staticmethod
    def _serialize(data):
        """Serialize data of types unknown to Python's `json` module."""
        if isinstance(data, np.ndarray):
            return data.tolist()

    def save(self, path=None, overwrite=True):
        """Save the configuration entries into a JSON file.

        Parameters
        ----------
        path : str or None
            Where to save. If None, save to the path from the initialization.
        overwrite : bool
            If True, it will overwrite if the path exists.

        """
        if path is None:
            if self.path is None:
                raise ValueError("No path passed.")
            path = self.path

        if os.path.exists(path) and not overwrite:
            raise PermissionError("JSON file not saved: overwrite not "
                                  "permitted.")

        with open(path, "wt") as f:
            json.dump(self.entries, f, sort_keys=True, indent=4,
                      ensure_ascii=True, default=self._serialize)

    def load(self, path=None, can_update=False):
        """Load the entries from a JSON file containing a dictionary.

        Parameters
        ----------
        path : str or None
            The path of the JSON file. If `None` it will use the path which
            which the ConfigFile instance was initialized.
        can_update : bool
            If True, if a key already exists, it will get the new value.
            If False, an Exception is thrown.

        """
        if path is None:
            if self.path is None:
                raise ValueError("No path passed.")
            path = self.path

        with open(path, "rt") as f:
            new_entries = json.load(f)

            if not can_update:
                current_keys = set(self.entries.keys())
                new_keys = set(new_entries.keys())
                common = list(current_keys & new_keys)
                if len(common) != 0:
                    raise PermissionError("Not allowed to update the entries"
                                          " {}".format(common))

            self.entries.update(new_entries)

    def __getattr__(self, key):
        """Return the value of an entry."""
        return self.entries[key]

    def __getitem__(self, key):
        """Return the value of an entry."""
        return self.entries[key]

    def __setitem__(self, key, value):
        """Create new or updates an entry."""
        self.entries[key] = value

    def __delitem__(self, key):
        """Delete an entry by it's key."""
        del self.entries[key]

    def __iter__(self):
        """Allow the iteration of entries."""
        return iter(self.entries)

    def update(self, dictionary):
        """Create new or update entries from an external dictionary."""
        self.entries.update(dictionary)

    def keys(self):
        """Return the keys of the configuration dictionary."""
        return self.entries.keys()

    def values(self):
        """Return the values of the configuration dictionary."""
        return self.entries.values()

    def items(self):
        """Return (key, value) tuples of the configuration dictionary."""
        return self.entries.items()

    def __repr__(self):
        """Represent the configuration dictionary as a string."""
        output_str = ''
        for key, value in self.entries.items():
            output_str = output_str + key + ": " + str(value) + "\n"
        return output_str

    def __contains__(self, item):
        """Search if a specific entry exists in the configuration."""
        return item in self.entries

    def __len__(self):
        """Return the number of configuration entries."""
        return len(self.entries)


def parse_inifile(inifile):
    """Parse an inifile and return dicts of each section."""
    binOps = {
        ast.Add: operator.add,
        ast.Sub: operator.sub,
        ast.Mult: operator.mul,
        ast.Div: operator.truediv,
        ast.Mod: operator.mod
    }

    def arithmetic_eval(s):
        """Control how the strings from the inifile get parsed."""
        node = ast.parse(s, mode='eval')

        def _eval(node):
            """Different strings receive different evaluation."""
            if isinstance(node, ast.Expression):
                return _eval(node.body)
            elif isinstance(node, ast.Str):
                if ',' in node.s:
                    return node.s.replace(' ', '').split(',')
                else:
                    return node.s
            elif isinstance(node, ast.Num):
                return node.n
            elif isinstance(node, ast.BinOp):
                return binOps[type(node.op)](_eval(node.left),
                                             _eval(node.right))
            elif isinstance(node, ast.List):
                return [_eval(x) for x in node.elts]
            elif isinstance(node, ast.Name): # pragma: no cover
                result = VariableKey(item=node)
                constants_lookup = {
                    'True': True,
                    'False': False,
                    'None': None,
                }
                value = constants_lookup.get(result.name, result,)
                if type(value) == VariableKey:
                    # return regular string
                    return value.name
                else:
                    # return special string like True or False
                    return value
            elif isinstance(node, ast.NameConstant):
                # None, True, False are nameconstants in python3 but names in 2
                return node.value
            else:
                raise Exception('Unsupported type {}'.format(node))

        return _eval(node.body)

    # ---- Create configuration-file-parser object and read parameters file.
    cp = configparser.ConfigParser(
        {'MESA_DIR': os.environ['MESA_DIR']},
        interpolation=configparser.ExtendedInterpolation()
    )
    cp.read(inifile)

    # ---- Read needed variables from the inifile
    dictionary = {}
    for section in cp.sections():
        dictionary[section] = {}
        for option in cp.options(section):
            opt = cp.get(section, option)
            try:
                try:
                    dictionary[section][option] = arithmetic_eval(opt)
                except Exception:
                    dictionary[section][option] = json.loads(opt)
            except Exception:
                if ',' in opt:
                    dictionary[section][option] = opt.replace(
                        ' ', '').split(',')
                else:
                    dictionary[section][option] = opt

    run_parameters = dictionary['run_parameters']
    mesa_inlists = dictionary['mesa_inlists']
    mesa_extras = dictionary['mesa_extras']
    slurm = dictionary['slurm']

    return run_parameters, slurm, mesa_inlists, mesa_extras


class VariableKey(object): # pragma: no cover
    """A dictionary key which is a variable.

    @ivar item: The variable AST object.
    """

    def __init__(self, item):
        """Construct the object by giving a `name` to it."""
        self.name = item.id

    def __eq__(self, compare):
        """Equality if the names are the same."""
        return (
            compare.__class__ == self.__class__
            and compare.name == self.name
        )

    def __hash__(self):
        """Allow hashing using the name of the variable."""
        return hash(self.name)
