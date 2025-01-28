"""The POSYDON exception class and subclasses for more specific errors."""


__authors__ = [
    "Camille Liotine <cliotine@u.northwestern.edu>",
    "Eirini Kasdagli <kasdaglie@ufl.edu>",
    "Jeff Andrews <jeffrey.andrews@ufl.edu>",
    "Konstantinos Kovlakas <Konstantinos.Kovlakas@unige.ch>",
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>",
]

class POSYDONError(Exception):
    """General POSYDON exception class."""

    def __init__(self, message=""):
        """General POSYDON exception.

        Parameters
        ----------
        message : str
            POSYDONError message.
        """
        if not isinstance(message, str):
            raise TypeError("The error message must be a string.")
        
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        """Create the text that accompanies this exception."""
        return super().__str__()

# Subclasses of POSYODONError in alphabetic order
# Before you add a new subclass check the list of python error classes at the
# end of this file
    
class ClassificationError(POSYDONError):
    """POSYDON error specific for binary classification errors."""

class DataError(POSYDONError):
    """POSYDON error specific for data errors."""

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
