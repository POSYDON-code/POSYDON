"""Unit tests of posydon/grids/downsampling.py
"""

__authors__ = [
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>"
]

# import the module which will be tested
import posydon.grids.downsampling as totest
# aliases
np = totest.np

# import other needed code for the tests, which is not already imported in the
# module you like to test
from pytest import fixture, raises, warns
from inspect import isclass, isroutine
from posydon.utils.posydonwarning import InappropriateValueWarning

# define test classes collecting several test functions
class TestElements:
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = {'Pwarn', 'TrackDownsampler', '__authors__',\
                    '__builtins__', '__cached__', '__doc__', '__file__',\
                    '__loader__', '__name__', '__package__', '__spec__', 'np',\
                    'rescale_from_0_to_1', 'sys'}
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

    def test_instance_rescale_from_0_to_1(self):
        assert isroutine(totest.rescale_from_0_to_1)

    def test_instance_TrackDownsampler(self):
        assert isclass(totest.TrackDownsampler)


class TestFunctions:
    # test functions
    def test_rescale_from_0_to_1(self):
        # missing argument
        with raises(TypeError, match="missing 1 required positional "\
                                     +"argument: 'array'"):
            totest.rescale_from_0_to_1()
        # bad input
        with raises(ValueError, match="Scaler works for matrices only."):
            totest.rescale_from_0_to_1(None)
        # examples
        test_values = [[0, 1, 3, -1], [0.5, 2, 3, -2], [0.4, 1.5, 3, -1]]
        norm_values = [[0.0, 0.0, 0.0, 1.0], [1.0, 1.0, 0.0, 0.0],\
                       [0.8, 0.5, 0.0, 1.0]]
        min_values = [0.0, 1.0, 3.0, -2.0]
        max_values = [0.5, 2.0, 3.0, -1.0]
        assert np.array_equal(totest.rescale_from_0_to_1(test_values),\
                              np.array(norm_values))
        arr, minima, maxima = totest.rescale_from_0_to_1(test_values,\
                                                         return_minmax=True)
        assert np.array_equal(arr, np.array(norm_values))
        assert np.array_equal(minima, np.array(min_values))
        assert np.array_equal(maxima, np.array(max_values))


class TestTrackDownsampler:
    @fixture
    def independent_sample(self):
        return np.array([0.0, 0.1, 1.0, 1.2])

    @fixture
    def dependent_sample(self):
        return np.array([[1.1, 0.0], [1.0, 0.1], [0.7, 0.5], [0.1, 0.1]])

    @fixture
    def norm_dependent_sample(self, dependent_sample):
        return totest.rescale_from_0_to_1(dependent_sample)

    @fixture
    def TrackDownsampler(self, independent_sample, dependent_sample):
        # initialize an instance of the class with defaults
        return totest.TrackDownsampler(independent_sample, dependent_sample)

    # test the TrackDownsampler class
    def test_init(self, independent_sample, dependent_sample,\
                  norm_dependent_sample, TrackDownsampler):
        assert isroutine(TrackDownsampler.__init__)
        # check that the instance is of correct type and all code in the
        # __init__ got executed: the elements are created and initialized
        assert TrackDownsampler.verbose == False
        assert TrackDownsampler.keep == None
        assert np.array_equal(TrackDownsampler.t, independent_sample)
        assert np.array_equal(TrackDownsampler.X, norm_dependent_sample)
        # examples: dependend is 1D
        test_TrackDownsampler = totest.TrackDownsampler([0.0, 0.1], [1.0, 2.0])
        assert np.array_equal(test_TrackDownsampler.t, np.array([0.0, 0.1]))
        assert np.array_equal(test_TrackDownsampler.X, np.array([[0.0],\
                                                                 [1.0]]))
        # bad input
        with raises(ValueError, match="`independent` is not iterable."):
            test_TrackDownsampler = totest.TrackDownsampler(None, None)
        with raises(ValueError, match="`independent` should be "\
                                      +"one-dimensional."):
            test_TrackDownsampler = totest.TrackDownsampler([[], []], None)
        with raises(IndexError, match="list index out of range"):
            test_TrackDownsampler = totest.TrackDownsampler([], None)
        with raises(ValueError, match="`dependent` is not iterable."):
            test_TrackDownsampler = totest.TrackDownsampler([0.0], None)
        with raises(ValueError, match="Number of dimensions in `dependent` "\
                                      +"> 2."):
            test_TrackDownsampler = totest.TrackDownsampler([0.0], [[[], []],\
                                                                    [[], []]])
        with raises(ValueError, match="`independent` should be strictly "\
                                      +"increasing."):
            test_TrackDownsampler = totest.TrackDownsampler([0.0, 0.0],\
                                                            [1.0, 2.0])

    def test_say(self, TrackDownsampler, capsys):
        assert isroutine(TrackDownsampler.say)
        # missing argument
        with raises(TypeError, match="missing 1 required positional "\
                                     +"argument: 'message'"):
            TrackDownsampler.say()
        # examples
        TrackDownsampler.say("Test")
        assert capsys.readouterr().out == ""
        TrackDownsampler.verbose = True
        TrackDownsampler.say("Test")
        assert capsys.readouterr().out == "Test\n"

    def test_rescale(self, dependent_sample, TrackDownsampler):
        assert isroutine(TrackDownsampler.rescale)
        # examples
        TrackDownsampler.X = dependent_sample
        TrackDownsampler.minima = None
        TrackDownsampler.maxima = None
        arr, minima, maxima = totest.rescale_from_0_to_1(dependent_sample,\
                                                         return_minmax=True)
        TrackDownsampler.rescale()
        assert np.array_equal(arr, TrackDownsampler.X)
        assert np.array_equal(minima, TrackDownsampler.minima)
        assert np.array_equal(maxima, TrackDownsampler.maxima)

    def test_extract_downsample(self, independent_sample, dependent_sample,\
                                norm_dependent_sample, TrackDownsampler):
        assert isroutine(TrackDownsampler.extract_downsample)
        # examples
        for k in [None, independent_sample<1.0]:
            TrackDownsampler.keep = k
            t, X = TrackDownsampler.extract_downsample()
            assert np.array_equal(t, independent_sample[k])
            assert np.array_equal(X, dependent_sample[k, :])
            t, X = TrackDownsampler.extract_downsample(scale_back=False)
            assert np.array_equal(t, independent_sample[k])
            assert np.array_equal(X, norm_dependent_sample[k, :])

    def test_find_downsample(self, TrackDownsampler, capsys):
        assert isroutine(TrackDownsampler.find_downsample)
        # examples: keep all
        TrackDownsampler.find_downsample()
        for k in TrackDownsampler.keep:
            assert k
        # examples: keep only edges
        with warns(InappropriateValueWarning, match="`max_err` >= 1.0 is "\
                                                    +"like disabling "\
                                                    +"downsampling."):
            TrackDownsampler.find_downsample(max_err=1.0)
            for (i, k) in enumerate(TrackDownsampler.keep):
                if ((i==0) or (i+1==len(TrackDownsampler.keep))):
                    assert k == True
                else:
                    assert k == False
        # examples: max_err can skip element 1; max_interval may keeps it
        for mi in [None, 1.0, 0.99, -0.85, -0.8]:
            TrackDownsampler.find_downsample(max_err=0.12, max_interval=mi)
            for (i, k) in enumerate(TrackDownsampler.keep):
                if ((i==1) and ((mi is None) or (mi>=1.0) or (mi<=-1.0/1.2))):
                    assert k == False
                else:
                    assert k == True
        # bad input
        with raises(ValueError, match="`max_err` must be a non-negative "\
                                      +"finite number."):
            TrackDownsampler.find_downsample(max_err=-1.0)
        with raises(ValueError, match="`max_err` must be a non-negative "\
                                      +"finite number."):
            TrackDownsampler.find_downsample(max_err=np.inf)

    def test_downsample(self, independent_sample, dependent_sample, TrackDownsampler):
        assert isroutine(TrackDownsampler.downsample)
        # examples
        t, X = TrackDownsampler.downsample()
        assert np.array_equal(t, independent_sample)
        assert np.array_equal(X, dependent_sample)
