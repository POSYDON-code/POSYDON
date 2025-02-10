"""Unit tests of posydon/popsyn/rate_calculation.py
"""

__authors__ = [
    "Elizabeth Teng <elizabethteng@u.northwestern.edu>"
]

# import the module which will be tested
import posydon.popsyn.rate_calculation as totest
# aliases
np = totest.np
sp = totest.sp

# import other needed code for the tests, which is not already imported in the
# module you like to test
from pytest import fixture, raises, warns, approx
from inspect import isroutine, isclass
from posydon.utils.constants import Zsun
from astropy.cosmology import Planck15 as cosmology
from astropy import constants as const
from astropy.cosmology import z_at_value
from scipy.interpolate import CubicSpline
from astropy import units as u


# define test classes collecting several test functions
class TestElements:
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = ['get_shell_comoving_volume', 'get_comoving_distance_from_redshift', 'get_cosmic_time_from_redshift', 'redshift_from_cosmic_time_interpolator',\
                    'get_redshift_from_cosmic_time','get_redshift_bin_edges',\
                    'get_redshift_bin_centers','__authors__',\
                    '__builtins__', '__cached__', '__doc__', '__file__',\
                    '__loader__', '__name__', '__package__', '__spec__']
        assert dir(totest) == elements, "There might be added or removed "\
                                        + "objects without an update on the "\
                                        + "unit test."

    def test_instance_get_shell_comoving_volume(self):
        assert isroutine(totest.GRB_selection)

    def test_instance_get_comoving_distance_from_redshift(self):
        assert isroutine(totest.chi_eff)

    def test_instance_get_cosmic_time_from_redshift(self):
        assert isroutine(totest.get_cosmic_time_from_redshift)

    def test_instance_redshift_from_cosmic_time_interpolator(self):
        assert isroutine(totest.redshift_from_cosmic_time_interpolator)

    def test_instance_get_redshift_from_cosmic_time(self):
        assert isroutine(totest.get_redshift_from_cosmic_time)

    def test_instance_get_redshift_bin_edges(self):
        assert isroutine(totest.get_redshift_bin_edges)
        
    def test_instance_get_redshift_bin_centers(self):
        assert isroutine(totest.get_redshift_bin_centers)

class TestFunctions:
    
    def test_get_shell_comoving_volume():
        
        
    def test_get_comoving_distance_from_redshift():
        
        
    def test_get_cosmic_time_from_redshift():
        
        
    def test_redshift_from_cosmic_time_interpolator():
        
    def test_get_redshift_from_cosmic_time():
        
    def test_get_redshift_bin_edges():
        
    def test_get_redshift_bin_centers():
        
        