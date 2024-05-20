""""The POSYDON warnings class and subclasses for more detailed warning handling."""


__authors__ = [
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>",
    "Monica Gallegos-Garcia <monicagallegosgarcia@u.northwestern.edu>",
    "Camille Liotine <cliotine@u.northwestern.edu>"
]


import warnings


class POSYDONWarning(Warning):
    """General POSYDON warning class."""

class ReplaceValueWarning(POSYDONWarning):
    """Warning that a value got replaced."""
    def __init__(self):
        warnings.filterwarnings("default", category=ReplaceValueWarning)

class InterpolationWarning(POSYDONWarning):
    """Warnings related to interpolation during binary evolution."""
    def __init__(self):
        warnings.filterwarnings("default", category=InterpolationWarning)

class ClassificationWarning(POSYDONWarning):
    """Warnings related to interpolation during binary evolution."""
    def __init__(self):
        warnings.filterwarnings("default", category=ClassificationWarning)

class ApproximationWarning(POSYDONWarning):
    """Warnings related to interpolation during binary evolution."""
    def __init__(self):
        warnings.filterwarnings("default", category=ApproximationWarning)

class UnsupportedModelWarning(POSYDONWarning):
    """Warnings related to selecting a model that is not supported."""
    def __init__(self):
        warnings.filterwarnings("default", category=UnsupportedModelWarning)

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
