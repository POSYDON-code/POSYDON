""""The POSYDON warnings class and subclasses for more detailed warning handling."""


__authors__ = [
    "Monica Gallegos-Garcia <monicagallegosgarcia@u.northwestern.edu>",
    "Camille Liotine <cliotine@u.northwestern.edu>",
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>",   
]


import warnings


class POSYDONWarning(Warning):
    """General POSYDON warning class."""
    def __init__(self, message=''):
        self.message = message
    
    def __str__(self):
        return repr(self.message)

class ApproximationWarning(POSYDONWarning):
    """Warning that a physical approximation was used during binary evolution."""
    def __init__(self, message=''):
        super().__init__(message)

class BinaryParsingWarning(POSYDONWarning):
    """Warnings related to parsing of binaries."""
    def __init__(self, message=''):
        super().__init__(message)

class ClassificationWarning(POSYDONWarning):
    """Warnings related to classification during binary evolution."""
    def __init__(self, message=''):
        super().__init__(message)

class EvolutionWarning(POSYDONWarning):
    """Warning that something unexpeted occurred during the binary evolution, but
    the evolution is able to continue (binary did not fail)."""
    def __init__(self, message=''):
        super().__init__(message)
        
class InappropriateValueWarning(POSYDONWarning):
    """Warnings that a strange value is used."""
    def __init__(self, message=''):
        super().__init__(message)

class IncompletenessWarning(POSYDONWarning):
    """Warnings when not all tasks could be done."""
    def __init__(self, message=''):
        super().__init__(message)

class InterpolationWarning(POSYDONWarning):
    """Warnings related to interpolation during binary evolution."""
    def __init__(self, message=''):
        super().__init__(message)

class MissingFilesWarning(POSYDONWarning):
    """Warnings related to missing files."""
    def __init__(self, message=''):
        super().__init__(message)

class OverwriteWarning(POSYDONWarning):
    """Warning that a data will get overwritten."""
    def __init__(self, message=''):
        super().__init__(message)

class ReplaceValueWarning(POSYDONWarning):
    """Warning that a value got replaced."""
    def __init__(self, message=''):
        super().__init__(message)

class UnsupportedModelWarning(POSYDONWarning):
    """Warnings related to selecting a model that is not supported."""
    def __init__(self, message=''):
        super().__init__(message)


# all POSYDON warnings subclasses should be defined beforehand
POSYDONWarning_subclasses = {cls.__name__: cls for cls in\
                             POSYDONWarning.__subclasses__()}

def _get_POSYDONWarning_class(category):
    """Inferring the POSYDONWarning class.

        Parameters
        ----------
        category : str or a warnings category class
            A string, which will be converted into a POSYDONWarning class.

        Returns
        -------
        POSYDONWarning or None
            A POSYDONWarning class or None if the input is neither a
            POSYDONWarning class nor convertable to one.
    """
    if category is None:
        return category
    if isinstance(category,str):
        if category in POSYDONWarning_subclasses.keys():
            category = POSYDONWarning_subclasses[category]
        else:
            category = POSYDONWarning
    if isinstance(category(), POSYDONWarning):
        return category
    else:
        return None

def Pwarn(message, category=None, stacklevel=2, **kwargs):
    """Issueing a warning via warnings.warn.

        Parameters
        ----------
        message : str
            The message printed in the warning.
        category : str or a warnings category class (default: None)
            A string, which will be converted into a POSYDONWarning class.
        stacklevel : int (default: 2)
            The stack level passed to warnings.warn, defaults to 2.
        **kwargs : dict (optional)
            Dictionary containing extra options passed to warnings.warn.
    """
    category = _get_POSYDONWarning_class(category)
    warnings.warn(message=message, category=category, stacklevel=stacklevel,
                  **kwargs)

def SetPOSYDONWarnings(action="default", category=POSYDONWarning, **kwargs):
    """Add the warnings filter for POSYDON warnings.

        Parameters
        ----------
        action : str (default: "default")
            The behaviour for those warnings. Should be one out of
            "default", "error", "ignore", "always", "module", or "once"
        category : str or a warnings category class (default: POSYDONWarning)
            A string, which will be converted into a POSYDONWarning class.
        **kwargs : dict (optional)
            Dictionary containing extra options passed to
            warnings.filterwarnings.
    """
    category = _get_POSYDONWarning_class(category)
    if isinstance(category(), POSYDONWarning):
        warnings.filterwarnings(action=action, category=category, **kwargs)

def NoPOSYDONWarnings(category=POSYDONWarning):
    """Switch the warnings filter to ignore for POSYDON warnings.

        Parameters
        ----------
        category : str or a warnings category class (default: POSYDONWarning)
            A string, which will be converted into a POSYDONWarning class.
    """
    SetPOSYDONWarnings(action="ignore", category=category)

def AllPOSYDONWarnings(category=POSYDONWarning):
    """Switch the warnings filter to always for POSYDON warnings.

        Parameters
        ----------
        category : str or a warnings category class (default: POSYDONWarning)
            A string, which will be converted into a POSYDONWarning class.
    """
    SetPOSYDONWarnings(action="always", category=category)


# The base class of all warnings is "Warning", which is derived from Exception
#
# List of Python's default warnings (sub)classes in alphabetic order
# BytesWarning: Base category for warnings related to bytes and bytearray.
# DeprecationWarning: Base category for warnings about deprecated features when
#     those warnings are intended for other Python developers (ignored by
#     default, unless triggered by code in __main__).
# EncodingWarning: Base class for warnings related to encodings.
# FutureWarning: Base category for warnings about deprecated features when
#     those warnings are intended for end users of applications that are
#     written in Python.
# ImportWarning: Base category for warnings triggered during the process of
#     importing a module (ignored by default).
# PendingDeprecationWarning: Base category for warnings about features that
#     will be deprecated in the future (ignored by default).
# ResourceWarning: Base category for warnings related to resource usage
#     (ignored by default).
# RuntimeWarning: Base category for warnings about dubious runtime features.
# SyntaxWarning: Base category for warnings about dubious syntactic features.
# UnicodeWarning: Base category for warnings related to Unicode.
# UserWarning: The default category for warn().
# 
# for more details see https://docs.python.org/3/library/exceptions.html#warnings
# and https://docs.python.org/3/library/warnings.html
