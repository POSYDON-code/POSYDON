import os
from dotenv import load_dotenv

load_dotenv()

def ensure_path(name):
    """Get the PATH set in the environment
    
    Parameters
    ----------
    name : str
        Name of the environment variable.
    
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
        raise NameError(f"{name} is not defined in the environment.")
    if not os.path.isdir(value):
        raise NotADirectoryError(f"{value} given in {name} is an invalid "
                                 "path.")
    return value

# POSYDON environment variables
PATH_TO_POSYDON = ensure_path("PATH_TO_POSYDON")
PATH_TO_POSYDON_DATA = os.path.join(ensure_path("PATH_TO_POSYDON_DATA"),
                                    "POSYDON_data")
