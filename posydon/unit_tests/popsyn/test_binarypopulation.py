"""Unit tests of posydon/popsyn/binarypopulation.py
"""

__authors__ = [
    "Elizabeth Teng <elizabethteng@u.northwestern.edu>"
]

# import the module which will be tested
import posydon.popsyn.binarypopulation as totest

# aliases
pd = totest.pd
np = totest.np

from inspect import isclass, isroutine

# import other needed code for the tests, which is not already imported in the
# module you like to test
from pytest import approx, fixture, raises, warns
