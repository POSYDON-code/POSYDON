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
traceback = totest.traceback
atexit = totest.atexit
os = totest.os
tqdm = totest.tqdm
psutil = totest.psutil
sys = totest.sys

BinaryStar = totest.BinaryStar
SingleStar = totest.SingleStar
properties_massless_remnant = totest.properties_massless_remnant
SimulationProperties = totest.SimulationProperties
get_formation_times = totest.get_formation_times
generate_independent_samples = totest.generate_independent_samples
binary_fraction_value = totest.binary_fraction_value
get_samples_from_file = totest.get_samples_from_file
get_kick_samples_from_file = totest.get_kick_samples_from_file
initial_total_underlying_mass = totest.initial_total_underlying_mass
default_kwargs = totest.default_kwargs

binarypop_kwargs_from_ini = totest.binarypop_kwargs_from_ini
Zsun = totest.Zsun
POSYDONError = totest.POSYDONError
Pwarn = totest.Pwarn
Catch_POSYDON_Warnings = totest.Catch_POSYDON_Warnings
orbital_period_from_separation = totest.orbital_period_from_separation
orbital_separation_from_period = totest.orbital_separation_from_period
set_binary_to_failed = totest.set_binary_to_failed


# import other needed code for the tests, which is not already imported in the
# module you like to test
from pytest import fixture, raises, warns, approx
from inspect import isroutine, isclass