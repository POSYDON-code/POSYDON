"""The POSYDON exception class and subclasses for more specific errors."""

import copy
from posydon.binary_evol.binarystar import BinaryStar
from posydon.binary_evol.singlestar import SingleStar


class POSYDONError(Exception):
    """General POSYDON exception class."""

    def __init__(self, message="", objects=None):
        """General POSYDON exception.

        Parameters
        ----------
        message : str
            POSYDONError message.
        objects : None or list of objects.
            A list of accompanied objects (or None if not set), that will be
            handled differently when `str` method is used. This error can
            only accept BinaryStar, SingleStar, or a list of each.

        """
        self.message = message
        # copy the objects: we must know their state at the moment of the error
        self.objects = copy.copy(objects)
        super().__init__(self.message)

    def __str__(self):
        """Create the text that accompanies this exception."""
        result = ""
        if self.objects is not None:
            if isinstance(self.objects,list):
                for i, obj in enumerate(self.objects):
                    if isinstance(obj, (BinaryStar, SingleStar)):
                        result += f"\n\nOBJECT #{i+1} ({type(obj)}):\n{str(obj)}"
            elif isinstance(obj, (BinaryStar, SingleStar)):
                result += f"\n\nOBJECT #({type(self.objects)}):\n{str(self.objects)}"
            else:
                pass
        return result + '\n'+ super().__str__()

class FileError(POSYDONError):
    """POSYDON error specific for reading/writing files."""

class FlowError(POSYDONError):
    """POSYDON error specific for binary evolution flow errors."""

class GridError(POSYDONError):
    """POSYDON error specific for PSyGrid operations."""

class MatchingError(POSYDONError):
    """POSYDON error specific for matching stars to single grid during the detached."""

class ModelError(POSYDONError):
    """POSYDON error specific for when a binary FAILS due to limitations of physics modeling assumptions."""

class NumericalError(POSYDONError):
    """POSYDON error specific for when a binary FAILS due to limitations of numerical methods."""

def initial_condition_message(binary,ini_params = None ):
    if ini_params is None:
        ini_params = ["\nFailed Binary Initial Conditions:\n",
                    f"S1 mass: {binary.star_1.mass_history[0]} \n",
                    f"S2 mass: {binary.star_2.mass_history[0]} \n",
                    f"S1 state: {binary.star_1.state_history[0]} \n",
                    f"S2 state: { binary.star_2.state_history[0]}\n",
                    f"orbital period: { binary.orbital_period_history[0] } \n",
                    f"eccentricity: { binary.eccentricity_history[0]} \n",
                    f"binary state: { binary.state_history[0] }\n",
                    f"binary event: { binary.state_history[0] }\n",
                    f"S1 natal kick array: { binary.star_1.natal_kick_array }\n",
                    f"S2 natal kick array: { binary.star_2.natal_kick_array}\n"]
    message = ""
    for i in ini_params:
        message +=  i 
    return message
