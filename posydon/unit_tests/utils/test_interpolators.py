"""Unit tests of posydon/utils/interpolators.py
"""

__authors__ = [
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>"
]

# import the module which will be tested
import posydon.utils.interpolators as totest

# aliases
np = totest.np

from inspect import isclass, isroutine

# import other needed code for the tests, which is not already imported in the
# module you like to test
from pytest import approx, fixture, raises, warns


# define test classes collecting several test functions
class TestElements:
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = {'interp1d', 'SingleStarInterpolator', '__authors__',\
                    '__builtins__', '__cached__', '__doc__', '__file__',\
                    '__loader__', '__name__', '__package__', '__spec__',\
                    'np', 'PchipInterpolator', 'copy'}
        totest_elements = set(dir(totest))
        missing_in_test = elements - totest_elements
        assert len(missing_in_test) == 0, "There are missing objects in "\
                                          +f"{totest.__name__}: "\
                                          +f"{missing_in_test}. Please "\
                                          +"check, whether they have been "\
                                          +"removed on purpose and update "\
                                          +"this unit test."
        new_in_test = totest_elements - elements
        assert len(new_in_test) == 0, "There are new objects in "\
                                      +f"{totest.__name__}: {new_in_test}. "\
                                      +"Please check, whether they have been "\
                                      +"added on purpose and update this "\
                                      +"unit test."

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
        assert np.allclose(interp1d.x, np.array(data[0]))
        assert np.allclose(interp1d.y, np.array(data[1]))
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
        # x not strictly monotonic
        example_data = ([0.0, 1.0, 0.5], [0.0, 0.5, 1.0])
        test_interp1d = totest.interp1d(example_data[0], example_data[1])
        assert np.allclose(test_interp1d.x, [0.0, 0.5, 1.0])
        assert np.allclose(test_interp1d.y, [0.0, 1.0, 0.5])

        # examples: reversible data
        test_interp1d = totest.interp1d(data[0][::-1], data[1][::-1])
        assert np.allclose(test_interp1d.x, np.array(data[0]))
        assert np.allclose(test_interp1d.y, np.array(data[1]))
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


class TestSingleStarInterpolator:
    @fixture
    def SSI_simple(self):
        # Single y_data, no positives, no derivatives
        t = np.array([0.0, 1.0])
        y_data = [np.array([1.0, 0.0])]
        y_keys = ["y"]
        return totest.SingleStarInterpolator(t, y_data, y_keys)

    @fixture
    def SSI_positive(self):
        t = np.array([0.0, 1.0])
        y_data = [np.array([-0.5, 0.5])]
        y_keys = ["y"]
        positives = [True]
        return totest.SingleStarInterpolator(t, y_data, y_keys, positives=positives)

    @fixture
    def SSI_derivative(self):
        t = np.array([0.0, 1.0])
        y_data = [np.array([-0.5, 0.5])]
        y_keys = ["y"]
        derivatives = [True]
        return totest.SingleStarInterpolator(t, y_data, y_keys, derivatives=derivatives)

    @fixture
    def SSI_multiple(self):
        t = np.array([0.0, 1.0])
        y_data = [np.array([1.0, 0.0]), np.array([-1.0, 2.0])]
        y_keys = ["a", "b"]
        positives = [False, True]
        derivatives = [False, False]
        return totest.SingleStarInterpolator(t, y_data, y_keys, positives=positives, derivatives=derivatives)

    @fixture
    def SSI_full_combinations(self):
        t = np.array([0.0, 1.0])
        y_data = [
            np.array([1.0, 2.0]),   # positive=False, derivative=False
            np.array([3.0, 1.0]),   # positive=True, derivative=False
            np.array([0.5, 0.0]),   # positive=False, derivative=True
            np.array([-1.0, 2.0])   # positive=True, derivative=True
        ]
        y_keys = ["a", "b", "c", "d"]
        positives = [False, True, False, True]
        derivatives = [False, False, True, True]
        return totest.SingleStarInterpolator(t, y_data, y_keys, positives=positives, derivatives=derivatives)

    @fixture
    def SSI_2D_coverage(self):
        t = np.array([0.0, 1.0])
        y_data = [
            [1.0, 2.0],  # positive=False, derivative=False
            [3.0, 4.0]   # positive=False, derivative=False
        ]
        y_keys = ["a", "b"]
        positives = [False, False]
        derivatives = [False, False]  # both same combination
        return SingleStarInterpolator(t, y_data, y_keys, positives=positives, derivatives=derivatives)

    # test the SingleStarInterpolator class
    def test_init(self, SSI_simple, SSI_positive, SSI_derivative, SSI_multiple):
        # Basic type checks
        for instance in [SSI_simple, SSI_positive, SSI_derivative, SSI_multiple]:
            assert isinstance(instance, totest.SingleStarInterpolator)
            assert hasattr(instance, "_interpolators")
            assert hasattr(instance, "_keys")
            assert hasattr(instance, "offset")

        # Check that keys match
        assert SSI_simple.keys == ["y"]
        assert SSI_positive.keys == ["y"]
        assert SSI_derivative.keys == ["y"]
        assert set(SSI_multiple.keys) == {"a", "b"}

    def test_call(self, SSI_simple, SSI_positive, SSI_derivative, SSI_multiple):
        # Single output, no positive, no derivative
        res = SSI_simple(0.1)
        assert isinstance(res, dict)
        assert "y" in res
        assert np.allclose(res["y"], 0.9)

        # Positive output
        res_pos = SSI_positive(0.1)
        assert np.allclose(res_pos["y"], 0.0)  # clipped to 0

        # Derivative output
        res_deriv = SSI_derivative(0.1)
        # derivative between -0.5 -> 0.5 is 1.0 slope
        assert np.allclose(res_deriv["y"], 1.0)

        # Multiple outputs
        res_multi = SSI_multiple(0.5)
        assert np.allclose(res_multi["a"], 0.5)  # linear interp 1 -> 0
        # 'b' is positive, should clip -0.5 to 0
        expected_b = max(-1.0 + (2.0 + 1.0) * 0.5, 0.0)  # linear interp
        assert np.allclose(res_multi["b"], expected_b)

    def test_offset(self, SSI_simple):
        SSI_simple.offset = 0.5
        res = SSI_simple(0.5)
        # effectively querying at t=0.0
        assert np.allclose(res["y"], 1.0)

    def test_full_combinations_call(self,SSI_full_combinations):
        res = SSI_full_combinations(0.5)
        assert set(res.keys()) == {"a", "b", "c", "d"}
        assert np.allclose(res["a"], 1.5)
        assert np.allclose(res["b"], 2.0)
        assert np.allclose(res["c"], -0.5)
        assert np.allclose(res["d"], 3.0)

    def test_values_2D_branch(self,SSI_2D_coverage):
        res = SSI_2D_coverage(0.5)
        assert set(res.keys()) == {"a", "b"}
