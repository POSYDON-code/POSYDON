"""Unit tests of posydon/popsyn/defaults.py
"""

__authors__ = [
    "Elizabeth Teng <elizabethteng@u.northwestern.edu>"
]

# import the module which will be tested
import posydon.popsyn.defaults as totest

# import other needed code for the tests, which is not already imported in the
# module you like to test
import pytest

# define test classes collecting several test functions
class TestElements:
    # check for objects, which should be an element of the tested module

    def test_dir(self):
        elements = ['default_kwargs', '__authors__',\
                    '__builtins__', '__cached__', '__doc__', '__file__',\
                    '__loader__', '__name__', '__package__', '__spec__','age_of_universe']
        totest_elements = set(dir(totest))
        missing_in_test = set(elements) - totest_elements
        assert len(missing_in_test) == 0, "There are missing objects in "\
                                          +f"{totest.__name__}: "\
                                          +f"{missing_in_test}. Please "\
                                          +"check, whether they have been "\
                                          +"removed on purpose and update "\
                                          +"this unit test."
        new_in_test = totest_elements - set(elements)
        assert len(new_in_test) == 0, "There are new objects in "\
                                      +f"{totest.__name__}: {new_in_test}. "\
                                      +"Please check, whether they have been "\
                                      +"added on purpose and update this "\
                                      +"unit test."
    
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
        assert set(totest.default_kwargs.keys()) == set(elements), \
            "The default_kwargs dictionary keys have changed. Please update the test."
    
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
        assert isinstance(totest.default_kwargs['binary_fraction_const'], (float, int)), \
            "binary_fraction_const should be a float or int"
        
    def test_instance_binary_fraction_scheme(self):
        assert isinstance(totest.default_kwargs['binary_fraction_scheme'], str), \
            "binary_fraction_scheme should be a string"