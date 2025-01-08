"""Unit tests of posydon/grids/scrubbing.py
"""

__authors__ = [
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>"
]

# import the module which will be tested
import posydon.grids.scrubbing as totest
# aliases
np = totest.np

# import other needed code for the tests, which is not already imported in the
# module you like to test
from pytest import fixture, raises, warns
from inspect import isroutine
#from posydon.grids.psygrid import PSyGrid
from posydon.utils.posydonwarning import InappropriateValueWarning

# define test classes collecting several test functions
class TestElements:
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = ['LG_MTRANSFER_RATE_THRESHOLD', 'Pwarn',\
                    'RL_RELATIVE_OVERFLOW_THRESHOLD',\
                    'THRESHOLD_CENTRAL_ABUNDANCE',\
                    'THRESHOLD_CENTRAL_ABUNDANCE_LOOSE_C', '__authors__',\
                    '__builtins__', '__cached__', '__doc__', '__file__',\
                    '__loader__', '__name__', '__package__', '__spec__',\
                    'keep_after_RLO', 'keep_till_central_abundance_He_C',\
                    'np', 'scrub']
        assert dir(totest) == elements, "There might be added or removed "\
                                        +"objects without an update on the "\
                                        +"unit test."

    def test_instance_scrub(self):
        assert isroutine(totest.scrub)

    def test_instance_keep_after_RLO(self):
        assert isroutine(totest.keep_after_RLO)

    def test_instance_keep_till_central_abundance_He_C(self):
        assert isroutine(totest.keep_till_central_abundance_He_C)


class TestFunctions:
    @fixture
    def tables(self):
        return np.array([1.0, 0.9, 0.8, 0.7, 0.6], dtype=([('mass', '<f8')]))

    @fixture
    def models(self):
        return np.array(range(5))

    @fixture
    def ages(self, models):
        return 10**models

    @fixture
    def star_history(self):
        # a temporary star history for testing
        return np.array([(1.0, 0.2, 0.0), (1.0e+2, 0.9, 0.0),\
                         (1.0e+3, 0.2, 0.8)],\
                        dtype=[('star_age', '<f8'), ('center_he4', '<f8'),\
                               ('center_c12', '<f8')])

    @fixture
    def star_history2(self):
        # a temporary star history for testing
        return np.array([(1.0, 1.0, 0.0), (1.0e+2, 0.0, 0.4),\
                         (1.0e+3, 0.0, 0.0)],\
                        dtype=[('star_age', '<f8'), ('center_he4', '<f8'),\
                               ('center_c12', '<f8')])

    @fixture
    def binary_history(self):
        # a temporary binary history for testing
        return np.array([(totest.RL_RELATIVE_OVERFLOW_THRESHOLD-1.0,\
                          totest.RL_RELATIVE_OVERFLOW_THRESHOLD-1.0,\
                          totest.LG_MTRANSFER_RATE_THRESHOLD-1.0, 1.0),\
                         (totest.RL_RELATIVE_OVERFLOW_THRESHOLD,\
                          totest.RL_RELATIVE_OVERFLOW_THRESHOLD-1.0,\
                          totest.LG_MTRANSFER_RATE_THRESHOLD, 1.0e+2),\
                         (totest.RL_RELATIVE_OVERFLOW_THRESHOLD-1.0,\
                          totest.RL_RELATIVE_OVERFLOW_THRESHOLD,\
                          totest.LG_MTRANSFER_RATE_THRESHOLD, 1.0e+3)],\
                        dtype=[('rl_relative_overflow_1', '<f8'),\
                               ('rl_relative_overflow_2', '<f8'),\
                               ('lg_mtransfer_rate', '<f8'), ('age', '<f8')])

    # test functions
    def test_scrub(self, tables, models, ages):
        # missing argument
        with raises(TypeError, match="missing 3 required positional "\
                                     +"arguments: 'tables', 'models', and "\
                                     +"'ages'"):
            totest.scrub()
        # bad input
        with raises(TypeError, match="'NoneType' object is not iterable"):
            totest.scrub(None, None, None)
        # examples: nothing
        assert totest.scrub([None], [None], [None]) == [None]
        assert totest.scrub([None], [models], [ages]) == [None]
        assert totest.scrub([tables], [None], [ages]) == [None]
        assert totest.scrub([tables], [models], [None]) == [None]
        # examples: no scrubbing for two tables and a None type object
        for (t, r) in zip(totest.scrub([tables, tables, None],\
                                       [models, models, None],\
                                       [ages, ages, None]),\
                                       [tables, tables, None]):
            if (isinstance(t, np.ndarray) and isinstance(r, np.ndarray)):
                assert np.array_equal(t, r)
            else:
                assert t == r
        # examples: scrubb element 1 on age
        ages[2] = ages[1]
        assert np.array_equal(totest.scrub([tables], [models], [ages])[0],\
               tables[np.array([True, False, True, True, True])])
        # examples: additionally scrubb element 3 on model
        models[4] = models[3]
        assert np.array_equal(totest.scrub([tables], [models], [ages])[0],\
               tables[np.array([True, False, True, False, True])])

    def test_keep_after_RLO(self, star_history, binary_history):
        # missing argument
        with raises(TypeError, match="missing 3 required positional "\
                                     +"arguments: 'bh', 'h1', and 'h2'"):
            totest.keep_after_RLO()
        # bad input
        with raises(ValueError, match="No `lg_mtransfer_rate` in binary "\
                                      +"history."):
            totest.keep_after_RLO(binary_history[['rl_relative_overflow_1']],\
                                  None, None)
        with raises(ValueError, match="No `rl_relative_overflow` of any star "\
                                      +"in binary history."):
            totest.keep_after_RLO(binary_history[['lg_mtransfer_rate']], None,\
                                  None)
        # examples: nothing
        for h1 in [None, np.array([])]:
            for h2 in [None, np.array([])]:
                assert totest.keep_after_RLO(None, h1, h2) == (None, h1, h2)
        # examples: cut out element 0 and history of star 1; it should be
        # noted, that the function although modifies the late part of input
        # objects
        bh, h1, h2 = totest.keep_after_RLO(binary_history[[\
         'rl_relative_overflow_1', 'lg_mtransfer_rate', 'age']], star_history,\
         None)
        assert np.array_equal(binary_history[['rl_relative_overflow_1',\
                                              'lg_mtransfer_rate',\
                                              'age']][1:], bh)
        assert np.array_equal(star_history[1:], h1)
        assert h2 is None
        # examples: cut out element 0 and history of star 2; test numerical
        # correction, incl. fail
        binary_history['rl_relative_overflow_1'][0] =\
         totest.RL_RELATIVE_OVERFLOW_THRESHOLD
        with raises(Exception, match="Numerical precision fix failed."):
            binary_history['age'] = np.array([0.3, 1.0, 1.0])
            totest.keep_after_RLO(binary_history[['rl_relative_overflow_1',\
                                                  'lg_mtransfer_rate',\
                                                  'age']], None, star_history)
        binary_history['age'] = np.array([0.3, 0.99999999999999994, 1.0])
        bh, h1, h2 = totest.keep_after_RLO(binary_history[[\
         'rl_relative_overflow_1', 'lg_mtransfer_rate', 'age']], None,\
         star_history)
        assert np.array_equal(binary_history[['rl_relative_overflow_1',\
                                              'lg_mtransfer_rate', 'age']], bh)
        assert h1 is None
        assert np.array_equal(star_history, h2)
        # examples: no RLO
        for i in range(len(binary_history)):
            binary_history['rl_relative_overflow_1'][i] =\
             totest.RL_RELATIVE_OVERFLOW_THRESHOLD - 1.0
            binary_history['rl_relative_overflow_2'][i] =\
             totest.RL_RELATIVE_OVERFLOW_THRESHOLD - 1.0
            binary_history['lg_mtransfer_rate'][i] =\
             totest.LG_MTRANSFER_RATE_THRESHOLD - 1.0
        assert totest.keep_after_RLO(binary_history, None, None) == None

    def test_keep_till_central_abundance_He_C(self, star_history,\
                                              star_history2, binary_history):
        # missing argument
        with raises(TypeError, match="missing 3 required positional "\
                                     +"arguments: 'bh', 'h1', and 'h2'"):
            totest.keep_till_central_abundance_He_C()
        # bad input
        with raises(AttributeError, match="'list' object has no attribute "\
                                          +"'dtype'"):
            totest.keep_till_central_abundance_He_C([], [], [])
        # examples: nothing
        for args in [(None, star_history, star_history),\
                     (binary_history, None, star_history),\
                     (binary_history, star_history, None),\
                     (binary_history[['lg_mtransfer_rate']],\
                      star_history[['center_he4']], star_history)]:
            bh, h1, h2, tf = totest.keep_till_central_abundance_He_C(*args)
            if args[0] is None:
                assert bh is None
            else:
                assert np.array_equal(bh, args[0])
            if args[1] is None:
                assert h1 is None
            else:
                assert np.array_equal(h1, args[1])
            if args[2] is None:
                assert h2 is None
            else:
                assert np.array_equal(h2, args[2])
            assert tf == ""
        for cols in [['star_age'], ['star_age', 'center_he4', 'center_c12']]:
            bh, h1, h2, tf = totest.keep_till_central_abundance_He_C(\
                              binary_history[['age']], star_history[cols],\
                              star_history[cols])
            assert np.array_equal(binary_history[['age']], bh)
            assert np.array_equal(star_history[cols], h1)
            assert np.array_equal(star_history[cols], h2)
            assert tf == ""
        # examples: one star is depleted
        bh, h1, h2, tf = totest.keep_till_central_abundance_He_C(\
                          binary_history[['age']], star_history2,\
                          star_history)
        assert np.array_equal(binary_history[['age']][:2], bh)
        assert np.array_equal(star_history2[:2], h1)
        assert np.array_equal(star_history[:2], h2)
        assert tf == "Primary got stopped before central carbon depletion"
        bh, h1, h2, tf = totest.keep_till_central_abundance_He_C(\
                          binary_history[['age']], star_history,\
                          star_history2)
        assert np.array_equal(binary_history[['age']][:2], bh)
        assert np.array_equal(star_history[:2], h1)
        assert np.array_equal(star_history2[:2], h2)
        assert tf == "Secondary got stopped before central carbon depletion"
        # examples: both stars are depleted
        bh, h1, h2, tf = totest.keep_till_central_abundance_He_C(\
                          binary_history[['age']], star_history2,\
                          star_history2)
        assert np.array_equal(binary_history[['age']][:2], bh)
        assert np.array_equal(star_history2[:2], h1)
        assert np.array_equal(star_history2[:2], h2)
        assert tf == "Primary got stopped before central carbon depletion"
        star_history['center_he4'][-1] = 0.0
        star_history['center_c12'][-1] = 0.0
        cols = ['center_he4', 'center_c12']
        bh, h1, h2, tf = totest.keep_till_central_abundance_He_C(\
                          binary_history[['age']], star_history[cols],\
                          star_history2[cols], XCstop=0.5)
        assert np.array_equal(binary_history[['age']][:2], bh)
        assert np.array_equal(star_history[cols][:2], h1)
        assert np.array_equal(star_history2[cols][:2], h2)
        assert tf == "Secondary got stopped before central carbon depletion"
        # examples: stars have history of length 0
        bh, h1, h2, tf = totest.keep_till_central_abundance_He_C(\
                          binary_history[['age']], star_history[0:0],\
                          star_history[0:0])
        assert np.array_equal(binary_history[['age']], bh)
        assert np.array_equal(star_history[0:0], h1)
        assert np.array_equal(star_history[0:0], h2)
        assert tf == ""
