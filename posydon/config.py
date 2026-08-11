"""Configuration helpers for environment paths in POSYDON."""

import os

from dotenv import load_dotenv

load_dotenv()


def ensure_path(name, default=None):
    """Get the 'PATH_*' variable set in the environment.

    Parameters
    ----------
    name : str
        Name of the environment variable.
    default : str or None
        Fallback path to use when the environment variable is missing or
        invalid. Only used for the local POSYDON repository defaults.

    Returns
    -------
    str
        Defined path or raises an error if the variable is not defined and
        pointing to a valid path.

    """
    if not isinstance(name, str):
        raise TypeError("'name' has to be a string.")

    value = os.getenv(name)
    if value is None:
        if default is not None:
            value = default
        else:
            raise NameError(f"{name} is not defined in the environment.")

    if not os.path.isdir(value):
        if default is not None and os.path.isdir(default):
            value = default
        else:
            raise NotADirectoryError(f"{value} given in {name} is an invalid "
                                     "path.")
    return value


# POSYDON environment variables
repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
PATH_TO_POSYDON = ensure_path("PATH_TO_POSYDON", default=repo_root)
path_to_posydon_data = ensure_path("PATH_TO_POSYDON_DATA",
                                   default=os.path.join(repo_root, "POSYDON_data"))
if path_to_posydon_data.endswith("POSYDON_data"):
    PATH_TO_POSYDON_DATA = path_to_posydon_data
else:
    PATH_TO_POSYDON_DATA = os.path.join(path_to_posydon_data, "POSYDON_data")
