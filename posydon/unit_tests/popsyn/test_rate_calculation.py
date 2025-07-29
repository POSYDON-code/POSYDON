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

# define test classes collecting several test functions
class TestElements:
    
    def test_defaults(self):
        elements = ["delta_t",
                    "SFR",
                    "sigma_SFR",
                    "Z_max",
                    "select_one_met",
                    "dlogZ",
                    "Zsun"]
        assert set(totest.DEFAULT_MODEL.keys()) == elements, \
            "The DEFAULT_MODEL dictionary keys have changed. Please update the test."    
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = ['get_shell_comoving_volume', 'get_comoving_distance_from_redshift', \
                    'get_cosmic_time_from_redshift', 'redshift_from_cosmic_time_interpolator',\
                    'get_redshift_from_cosmic_time','get_redshift_bin_edges',\
                    'get_redshift_bin_centers','__authors__',\
                    '__builtins__', '__cached__', '__doc__', '__file__',\
                    '__loader__', '__name__', '__package__', '__spec__']
        assert dir(totest) == elements, "There might be added or removed "\
                                        + "objects without an update on the "\
                                        + "unit test."

    def test_instance_get_shell_comoving_volume(self):
        assert isroutine(totest.get_shell_comoving_volume)

    def test_instance_get_comoving_distance_from_redshift(self):
        assert isroutine(totest.get_comoving_distance_from_redshift)

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
    
    # test functions
    def test_get_shell_comoving_volume(self):
        # 2 missing arguments
        with raises(TypeError, match="missing 2 required positional arguments: 'z_hor_i' and 'z_hor_f'"):
            totest.test_get_shell_comoving_volume()
        # 1 missing argument
        with raises(TypeError, match="missing 1 required positional argument: 'z_hor_f'"):
            totest.test_get_shell_comoving_volume(0.1)
        # examples
        tests = [(0.1, 1.0, approx(97.7972132977263, abs=6e-12)),\
                 (0.3, 2.0, approx(277.8780499884267, abs=6e-12))]
        for (z1, z2, v) in tests:
            assert totest.test_get_shell_comoving_volume(z1, z2) == v
            
    def test_get_comoving_distance_from_redshift(self):
        # missing argument
        with raises(TypeError, match="missing 1 required positional argument: 'z'"):
            totest.get_comoving_distance_from_redshift()
        # examples
        tests = [(0.1, approx(432.1244883487781, abs=6e-12)),\
                 (1.0, approx(3395.905311975348, abs=6e-12))]
        for (z, d) in tests:
            assert totest.get_comoving_distance_from_redshift(z) == d
        
    def test_get_cosmic_time_from_redshift(self):
        # missing argument
        with raises(TypeError, match="missing 1 required positional argument: 'z'"):
            totest.test_get_cosmic_time_from_redshift()
        # examples
        tests = [(0.1, approx(12.453793290949799, abs=6e-12)),\
                 (1.0, approx(5.862549255024051, abs=6e-12))]
        for (z, t) in tests:
            assert totest.test_get_cosmic_time_from_redshift(z) == t
        
    def test_redshift_from_cosmic_time_interpolator(self):
        interp = totest.redshift_from_cosmic_time_interpolator()
        assert isinstance(interp, CubicSpline) 
        
    def test_get_redshift_from_cosmic_time(self):
        # missing argument
        with raises(TypeError, match="missing 1 required positional argument: 't_cosm'"):
            totest.test_get_redshift_from_cosmic_time()
        # examples
        tests = [(0.1, approx(29.8325299, abs=6e-12)),\
                 (1.0, approx(5.67584779, abs=6e-12))]
        for (t, z) in tests:
            assert totest.test_get_cosmic_time_from_redshift(t) == z        
        
    def test_get_redshift_bin_edges(self):
        # missing argument
        with raises(TypeError, match="missing 1 required positional argument: 'delta_t'"):
            totest.test_get_redshift_bin_edges()
        # examples
        tests = [(100., approx(0.006963184181145605, abs=6e-12)),\
                 (1000., approx(0.07301543666184201, abs=6e-12))]
        for (t,arr) in tests:
            assert totest.test_get_redshift_bin_edges(t)[1] == arr
        
    def test_get_redshift_bin_centers(self):
        # missing argument
        with raises(TypeError, match="missing 1 required positional argument: 'delta_t'"):
            totest.test_get_redshift_bin_centers()
        # examples
        tests = [(100., approx(49.33542627789386, abs=6e-12)),\
                 (1000., approx(13.957133275502315, abs=6e-12))]
        for (t,arr) in tests:
            assert totest.test_get_redshift_bin_centers(t)[-1] == arr
        