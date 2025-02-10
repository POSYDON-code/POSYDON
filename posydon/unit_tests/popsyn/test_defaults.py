"""Unit tests of posydon/popsyn/defaults.py
"""

__authors__ = [
    "Elizabeth Teng <elizabethteng@u.northwestern.edu>"
]

# import the module which will be tested
import posydon.popsyn.analysis as totest
from posydon.utils.constants import age_of_universe

# import other needed code for the tests, which is not already imported in the
# module you like to test
import pytest

# define test classes collecting several test functions
class TestElements:
    # check for objects, which should be an element of the tested module

    def test_dir(self):
        elements = ['default_kwargs', '__authors__',\
                    '__builtins__', '__cached__', '__doc__', '__file__',\
                    '__loader__', '__name__', '__package__', '__spec__']
        assert dir(totest) == elements, "There might be added or removed "\
                                        + "objects without an update on the "\
                                        + "unit test."
    
    def test_kwargs(self):
        elements = [
            'entropy',
            'number_of_binaries',
            'metallicity',
            'star_formation',
            'max_simulation_time',
            'orbital_scheme',
            'orbital_separation_scheme',
            'orbital_separation_min',
            'orbital_separation_max',
            'log_orbital_seperation_mean',
            'log_orbital_seperation_sigma',
            'orbital_period_scheme',
            'orbital_period_min',
            'orbital_period_max',
            'eccentricity_scheme',
            'primary_mass_scheme',
            'primary_mass_min',
            'primary_mass_max',
            'secondary_mass_scheme',
            'secondary_mass_min',
            'secondary_mass_max',
            'binary_fraction_const',
            'binary_fraction_scheme'
        ]
        assert set(totest.default_kwargs.keys()) == elements, \
            "The default_kwargs dictionary keys have changed. Please update the test."
    
    def test_instance_default_kwargs(self):
        assert isinstance(totest.default_kwargs, (dir)),\
                "default_kwargs is of type: "\
                + str(type(totest.default_kwargs))

    def test_instance_entropy(self):
        assert isinstance(totest.default_kwargs['entropy'], (type(None), float)), \
            "entropy should be None or a float"
        
    def test_instance_number_of_binaries(self):
        assert isinstance(totest.default_kwargs['number_of_binaries'], int), \
            "number_of_binaries should be an integer"

    def test_instance_metallicity(self):
        assert isinstance(totest.default_kwargs['metallicity'], float), \
            "metallicity should be a float"
        
    def test_instance_star_formation(self):
        assert isinstance(totest.default_kwargs['star_formation'], str), \
            "star_formation should be a string"

    def test_instance_max_simulation_time(self):
        assert isinstance(totest.default_kwargs['max_simulation_time'], (float, int)), \
            "max_simulation_time should be a float or int"
        assert totest.default_kwargs['max_simulation_time'] == age_of_universe, \
            "max_simulation_time should be equal to age_of_universe"
        
    def test_instance_orbital_scheme(self):
        assert isinstance(totest.default_kwargs['orbital_scheme'], str), \
            "orbital_scheme should be a string"
        
    def test_instance_orbital_separation_scheme(self):
        assert isinstance(totest.default_kwargs['orbital_separation_scheme'], str), \
            "orbital_scheme should be a string"
        
    def test_instance_orbital_separation_min(self):
        assert isinstance(totest.default_kwargs['orbital_separation_min'], float), \
            "orbital_separation_min should be a float"
        
    def test_instance_orbital_separation_max(self):
        assert isinstance(totest.default_kwargs['orbital_separation_max'], float), \
            "orbital_separation_max should be a float"
        
    def test_instance_log_orbital_seperation_mean(self):
        assert isinstance(totest.default_kwargs['log_orbital_seperation_mean'], (type(None), float)), \
            "log_orbital_seperation_mean should be None or a float"
        
    def test_instance_log_orbital_seperation_sigma(self):
        assert isinstance(totest.default_kwargs['log_orbital_seperation_sigma'], (type(None), float)), \
            "log_orbital_seperation_sigma should be None or a float"
        
    def test_instance_orbital_period_min(self):
        assert isinstance(totest.default_kwargs['orbital_period_min'], float), \
            "orbital_period_min should be a float"
        
    def test_instance_orbital_period_max(self):
        assert isinstance(totest.default_kwargs['orbital_period_max'], float), \
            "orbital_period_max should be a float"
        
    def test_instance_eccentricity_scheme(self):
        assert isinstance(totest.default_kwargs['eccentricity_scheme'], str), \
            "eccentricity_scheme should be a string"
        
    def test_instance_primary_mass_min(self):
        assert isinstance(totest.default_kwargs['primary_mass_min'], float), \
            "primary_mass_min should be a float"
        
    def test_instance_primary_mass_max(self):
        assert isinstance(totest.default_kwargs['primary_mass_max'], float), \
            "primary_mass_max should be a float"
        
    def test_instance_secondary_mass_min(self):
        assert isinstance(totest.default_kwargs['secondary_mass_min'], float), \
            "secondary_mass_min should be a float"
        
    def test_instance_secondary_mass_max(self):
        assert isinstance(totest.default_kwargs['secondary_mass_max'], float), \
            "secondary_mass_max should be a float"
        
    def test_instance_binary_fraction_const(self):
        assert isinstance(totest.default_kwargs['binary_fraction_const'], int), \
            "binary_fraction_const should be an integer"
        
    def test_instance_binary_fraction_scheme(self):
        assert isinstance(totest.default_kwargs['binary_fraction_scheme'], str), \
            "binary_fraction_scheme should be a string"

class TestValues:
    # check that the values fit

    def test_value_(self):
        assert totest.default_kwargs['number_of_binaries'] == 100, \
            "number_of_binaries default value has changed"

    def test_value_(self):
        assert totest.default_kwargs['metallicity'] == 1.0, \
            "metallicity default value has changed"

    def test_value_(self):
        assert totest.default_kwargs['star_formation'] == 'constant', \
            "star_formation default value has changed"

    def test_value_(self):
        assert totest.default_kwargs['max_simulation_time'] == age_of_universe, \
            "max_simulation_time default value should be age_of_universe"

    def test_value_(self):
        assert totest.default_kwargs['orbital_scheme'] == 'period', \
            "orbital_scheme default value has changed"

    def test_value_(self):
        assert totest.default_kwargs['orbital_separation_min'] == 5.0, \
            "orbital_separation_min default value has changed"

    def test_value_(self):
        assert totest.default_kwargs['orbital_separation_max'] == 1.0e5, \
            "orbital_separation_max default value has changed"

    def test_value_(self):
        assert totest.default_kwargs['orbital_period_min'] == 0.75, \
            "orbital_period_min default value has changed"

    def test_value_(self):
        assert totest.default_kwargs['orbital_period_max'] == 6000, \
            "orbital_period_max default value has changed"

    def test_value_(self):
        assert totest.default_kwargs['eccentricity_scheme'] == 'zero', \
            "eccentricity_scheme default value has changed"

    def test_value_(self):
        assert totest.default_kwargs['primary_mass_min'] == 7.0, \
            "primary_mass_min default value has changed"

    def test_value_(self):
        assert totest.default_kwargs['primary_mass_max'] == 120.0, \
            "primary_mass_max default value has changed"

    def test_value_(self):
        assert totest.default_kwargs['secondary_mass_min'] == 0.35, \
            "secondary_mass_min default value has changed"

    def test_value_(self):
        assert totest.default_kwargs['secondary_mass_max'] == 120.0, \
            "secondary_mass_max default value has changed"

    def test_value_(self):
        assert totest.default_kwargs['binary_fraction_const'] == 1, \
            "binary_fraction_const default value has changed"

    def test_value_(self):
        assert totest.default_kwargs['binary_fraction_scheme'] == 'const', \
            "binary_fraction_scheme default value has changed"

