"""The POSYDON exception class and subclasses for more specific errors."""


__authors__ = [
    "Camille Liotine <cliotine@u.northwestern.edu>",
    "Eirini Kasdagli <kasdaglie@ufl.edu>",
    "Konstantinos Kovlakas <Konstantinos.Kovlakas@unige.ch>",
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>",
]


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
        objects : None or list of objects
            A list of accompanied objects (or None if not set), that will be
            handled differently when `str` method is used. This error can
            only accept BinaryStar, SingleStar, or a list of each.

        """
        if not isinstance(message, str):
            raise TypeError("The error message must be a string.")
        if ((objects is not None) and
            (not isinstance(objects, (list, SingleStar, BinaryStar)))):
            raise TypeError("The error objects must be None, a list, a "
                            "SingleStar object, or a BinaryStar object.")
        self.message = message
        # copy the objects: we must know their state at the moment of the error
        self.objects = copy.copy(objects)
        super().__init__(self.message)

    def __str__(self):
        """Create the text that accompanies this exception."""
        result = ""
        if self.objects is not None:
            if isinstance(self.objects, list):
                for i, obj in enumerate(self.objects):
                    if isinstance(obj, (BinaryStar, SingleStar)):
                        result += f"\n\nOBJECT #{i+1} ({type(obj)}):\n{str(obj)}"
            elif isinstance(self.objects, (BinaryStar, SingleStar)):
                result += f"\n\nOBJECT #({type(self.objects)}):\n{str(self.objects)}"
            else: # pragma: no cover
                pass
        return result + '\n'+ super().__str__()

# Subclasses of POSYODONError in alphabetic order
# Before you add a new subclass check the list of python error classes at the
# end of this file

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

def initial_condition_message(binary, ini_params=None):
    """Generate a message with the initial conditions.

        Parameters
        ----------
        binary : BinaryStar
            BinaryStar object to take the initial conditions from.
        ini_params : None or iterable of str
            If None take the initial conditions from the binary, otherwise add
            each item of it to the message.

        Returns
        -------
        string
            The message with the initial conditions.

    """
    if not isinstance(binary, BinaryStar):
        raise TypeError("The binary must be a BinaryStar object.")
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
        message += i 
    return message


# There are 5 base exception classes in python: BaseException, Exception,
# ArithmeticError, BufferError, LookupError. Those are not aimed to be raised
# directly but their subclasses.
#
# List of Python's default exception (sub)classes in alphabetic order
# AssertionError: Raised when an assert statement fails.
# AttributeError: Raised when an attribute reference or assignment fails.
# EOFError: Raised when the input() function hits an end-of-file condition
#     (EOF) without reading any data.
# GeneratorExit: Raised when a generator or coroutine is closed.
# ImportError: Raised when the import statement has troubles trying to load a
#     module or a requested object does not exist in the module.
#     "ModuleNotFoundError" is a subclass of this error.
# IndexError: Raised when a sequence subscript is out of range.
# KeyError: Raised when a mapping key is not found in the set of existing keys.
# KeyboardInterrupt: Raised when the user hits the interrupt key (normally
#     Control-C or Delete).
# MemoryError: Raised when an operation runs out of memory.
# NameError: Raised when a local or global name is not found.
#     "UnboundLocalError" is a subclass of this error.
# RuntimeError: Raised when an error is detected that doesn’t fall in any of
#     the other categories. "NotImplementedError"(indicates a planned but not
#     yet implemented class/function) and "RecursionError"(reached maximum
#     recursion depth) are subclasses of this error.
# OSError: This exception is raised when a system function returns a
#     system-related error. This error has different names in old or platform
#     specfic versions like "EnvironmentError", "IOError", "WindowsError".
#     There is a bunch of subclasses: "BlockingIOError", "ChildProcessError",
#     "ConnectionError", "BrokenPipeError", "ConnectionAbortedError",
#     "ConnectionRefusedError", "ConnectionResetError", "FileExistsError",
#     "FileNotFoundError", "InterruptedError", "IsADirectoryError",
#     "NotADirectoryError", "PermissionError", "ProcessLookupError",
#     "TimeoutError".
# OverflowError: Raised when the result of an arithmetic operation is too large
#     to be represented.
# ReferenceError: Raised when a reference is no longer valid because the
#     referred object is garbage collected.
# StopIteration: Raised if there is no next item, while requested.
# StopAsyncIteration: Raised to stop asynchronous interation.
# SyntaxError: Raised when the parser encounters a syntax error.
#     "IndentationError"(has an own subclass TabError) is a subclasses of this
#     error.
# SystemError: Raised when the interpreter finds an internal error.
# SystemExit: This exception is subclass of BaseException instead of Exception
#     to be not caugth with other exceptions.
# TypeError: Raised when an operation or function is applied to an object of
#     inappropriate type.
# ValueError: Raised when an operation or function receives an argument that
#     has the right type but an inappropriate value. "UnicodeError"(has own
#     subclasses: "UnicodeEncodeError", "UnicodeDecodeError",
#     "UnicodeTranslateError") is a subclass of this error.
# ZeroDivisionError: Raised when the second argument of a division or modulo
#     operation is zero.
#
# for more details see https://docs.python.org/3/library/exceptions.html
