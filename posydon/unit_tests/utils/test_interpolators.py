"""Unit tests of posydon/utils/interpolators.py
"""

__authors__ = [
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>"
]

# import the module which will be tested
import posydon.utils.interpolators as totest
# aliases
np = totest.np

# import other needed code for the tests, which is not already imported in the
# module you like to test
from pytest import fixture, raises, warns, approx
from inspect import isroutine, isclass

# define test classes collecting several test functions
class TestElements:
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = ['PchipInterpolator', 'PchipInterpolator2', '__authors__',\
                    '__builtins__', '__cached__', '__doc__', '__file__',\
                    '__loader__', '__name__', '__package__', '__spec__',\
                    'interp1d', 'np']
        assert dir(totest) == elements, "There might be added or removed "\
                                        + "objects without an update on the "\
                                        + "unit test."

    def test_instance_interp1d(self):
        assert isclass(totest.interp1d)

    def test_instance_PchipInterpolator2(self):
        assert isclass(totest.PchipInterpolator2)


class Testinterp1d:
    @fixture
    def data(self):
        # test data
        return ([0.0, 1.0], [1.0, 0.0])

    @fixture
    def interp1d(self, data):
        # initialize an instance of the class with defaults
        return totest.interp1d(data[0], data[1])

    # test the interp1d class
    def test_init(self, data, interp1d):
        assert isroutine(interp1d.__init__)
        # check that the instance is of correct type and all code in the
        # __init__ got executed: the elements are created and initialized
        assert isinstance(interp1d, totest.interp1d)
        assert np.array_equal(interp1d.x, np.array(data[0]))
        assert np.array_equal(interp1d.y, np.array(data[1]))
        assert interp1d.kind == 'linear'
        assert interp1d.below is None
        assert interp1d.above is None
        assert interp1d.extrapolate == False
        # bad input
        with raises(TypeError, match="missing 2 required positional "\
                                     +"arguments: 'x' and 'y'"):
            totest.interp1d()
        with raises(NotImplementedError, match="kind = test is not supported"):
            totest.interp1d(x=data[0], y=data[1], kind='test')
        with raises(ValueError):
            totest.interp1d(None, None)
        with raises(ValueError):
            totest.interp1d('test', 'test')
        with raises(NotImplementedError, match="fill_value has to be "\
                                                +"'extrapolate' or a tuple "\
                                                +"with 1 or 2 elements"):
            totest.interp1d(x=data[0], y=data[1], fill_value='test')
        with raises(TypeError):
            totest.interp1d(x=data[0], y=data[1], fill_value=('test', 'test'))
        with raises(TypeError):
            totest.interp1d(x=data[0], y=data[1], left='test')
        with raises(TypeError):
            totest.interp1d(x=data[0], y=data[1], right='test')
        
        # incorrect input arrays
        x_data = [0.0, 1.0, 2.0, 1.5]
        y_data = [0.0, 1.0, 0.0, 0.5]
        with raises(ValueError, match="x values must be strictly increasing or strictly decreasing."):
            interpolator = totest.interp1d(x_data, y_data)
        
        # invertible data
        interpolator = totest.interp1d(data[0][::-1], data[1][::-1])
        assert np.allclose(interpolator.x, data[0])
        assert np.allclose(interpolator.y, data[1])
        
        # examples
        kinds = ['linear']
        for k in kinds:
            test_interp1d = totest.interp1d(x=data[0], y=data[1], kind=k)
            assert test_interp1d.kind == k
        # examples
        test_interp1d = totest.interp1d(x=data[0], y=data[1],\
                                        fill_value='extrapolate')
        assert test_interp1d.extrapolate == True
        # examples
        tests = [({'fill_value': (2.0,)}, 2.0, 2.0),\
                 ({'fill_value': (2.0, 3.0)}, 2.0, 3.0),\
                 ({'left': 2.0}, 2.0, None),\
                 ({'right': 3.0}, None, 3.0),\
                 ({'left': 2.0, 'right': 3.0}, 2.0, 3.0),\
                 ({'fill_value': (1.2, 1.3), 'left': 2.0, 'right': 3.0},\
                  2.0, 3.0)]
        for (kwargs, b, a) in tests:
            test_interp1d = totest.interp1d(x=data[0], y=data[1], **kwargs)
            if b is None:
                assert test_interp1d.below is None
            else:
                assert test_interp1d.below == b
            if a is None:
                assert test_interp1d.above is None
            else:
                assert test_interp1d.above == a

    def test_call(self, data, interp1d):
        assert isroutine(interp1d.__call__)
        # examples
        for x,y in zip(data[0], data[1]):
            assert interp1d(x) == y
        tests = [(0.0, 1.0), (0.2, approx(0.8, abs=6e-12)),\
                 (0.5, approx(0.5, abs=6e-12)), (0.9, approx(0.1, abs=6e-12)),\
                 (1.4, 0.0), (-0.1, 1.0)]
        for (x,y) in tests:
            assert interp1d(x) == y
        # examples: extrapolate
        test_interp1d = totest.interp1d(x=data[0], y=data[1],\
                                        fill_value='extrapolate')
        tests[-1] = (-0.1, approx(1.1, abs=6e-12))
        tests[-2] = (1.4, approx(-0.4, abs=6e-12))
        for (x,y) in tests:
            assert test_interp1d(x) == y
        # bad value
        interp1d.kind = 'test'
        with raises(NotImplementedError, match="kind = test is not supported"):
            interp1d(0.2)



class TestPchipInterpolator2:
    @fixture
    def PchipInterpolator2(self):
        # initialize an instance of the class with defaults
        return totest.PchipInterpolator2([0.0, 1.0], [1.0, 0.0])

    @fixture
    def PchipInterpolator2_True(self):
        # initialize an instance of the class with defaults
        return totest.PchipInterpolator2([0.0, 1.0], [-0.5, 0.5],\
                                         positive=True)

    @fixture
    def PchipInterpolator2_False(self):
        # initialize an instance of the class with defaults
        return totest.PchipInterpolator2([0.0, 1.0], [-0.5, 0.5],\
                                         positive=False)

    # test the PchipInterpolator2 class
    def test_init(self, PchipInterpolator2, PchipInterpolator2_True,\
                  PchipInterpolator2_False):
        assert isroutine(PchipInterpolator2.__init__)
        # check that the instance is of correct type and all code in the
        # __init__ got executed: the elements are created and initialized
        assert isinstance(PchipInterpolator2, totest.PchipInterpolator2)
        assert isinstance(PchipInterpolator2.interpolator,\
                          totest.PchipInterpolator)
        assert PchipInterpolator2.positive == False
        assert PchipInterpolator2_True.positive == True
        assert PchipInterpolator2_False.positive == False

    def test_call(self, PchipInterpolator2, PchipInterpolator2_True,\
                  PchipInterpolator2_False):
        assert isroutine(PchipInterpolator2.__call__)
        assert PchipInterpolator2(0.1) == 0.9
        assert PchipInterpolator2_True(0.1) == 0.0
        assert PchipInterpolator2_False(0.1) == -0.4
        assert np.allclose(PchipInterpolator2([0.1, 0.8]),\
                           np.array([0.9, 0.2]))
        assert np.allclose(PchipInterpolator2_True([0.1, 0.8]),\
                           np.array([0.0, 0.3]))
        assert np.allclose(PchipInterpolator2_False([0.1, 0.8]),\
                           np.array([-0.4, 0.3]))

