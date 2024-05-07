""""The POSYDON warnings class and subclasses for more detailed warning handling."""

import warnings

__authors__ = [
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>",
]

class POSYDONWarning(Warning):
    """General POSYDON warning class."""

class ReplaceValueWarning(POSYDONWarning):
    """Warning that a value got replaced."""
    def __init__(self):
        warning.filterwarnings("default", category=ReplaceValueWarning)
