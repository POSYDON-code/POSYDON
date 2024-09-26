"""Unit tests of posydon/utils/ignorereason.py
"""

__authors__ = [
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>"
]

# import the module which will be tested
import posydon.utils.ignorereason as totest

# import other needed code for the tests, which is not already imported in the
# module you like to test
from pytest import fixture, raises
from inspect import isclass, isroutine

# define test classes collecting several test functions
class TestElements:
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = ['IGNORE_REASONS_PRIORITY', 'IgnoreReason', '__authors__',
                    '__builtins__', '__cached__', '__doc__', '__file__',
                    '__loader__', '__name__', '__package__', '__spec__']
        assert dir(totest) == elements, "There might be added or removed "+\
               "objects without an update on the unit test."

    def test_instance_IGNORE_REASONS_PRIORITY(self):
        assert isinstance(totest.IGNORE_REASONS_PRIORITY, (list)),\
               "IGNORE_REASONS_PRIORITY is of type: "+\
               str(type(totest.IGNORE_REASONS_PRIORITY))

    def test_instance_IgnoreReason(self):
        assert isclass(totest.IgnoreReason)


class TestValues:
    # check that the values fit
    def test_value_IGNORE_REASONS_PRIORITY(self):
        v_last = None
        for v in ['ignored_no_history1', 'ignored_no_binary_history',
                  'corrupted_history1', 'corrupted_binary_history',
                  'corrupted_history2', 'ignored_scrubbed_history',
                  'ignored_no_final_profile', 'ignored_no_RLO']:
            # check required values
            assert v in totest.IGNORE_REASONS_PRIORITY, "missing entry"
            # check required order
            if v_last is not None:
                assert totest.IGNORE_REASONS_PRIORITY.index(v_last) <\
                       totest.IGNORE_REASONS_PRIORITY.index(v),\
                       f"the priority order has changed for: {v_last} or {v}"
            v_last = v


class TestIgnoreReason:
    # test the IgnoreReason class
    @fixture
    def NewInstance(self):
        # initialize an instance of the class for each test
        return totest.IgnoreReason()

    def test_init(self, NewInstance):
        assert isroutine(NewInstance.__init__)
        # check that the instance is of correct type and all code in the
        # __init__ got executed: the elements are created and initialized
        assert isinstance(NewInstance, totest.IgnoreReason)
        # strictly speaking, the following test does not test the variable
        # asignments in the __init__, because they call the __setattr__
        assert NewInstance.reason is None
        assert NewInstance.order is None

    def test_bool(self, NewInstance):
        assert isroutine(NewInstance.__bool__)
        assert bool(NewInstance) == False
        NewInstance.reason = totest.IGNORE_REASONS_PRIORITY[0]
        assert bool(NewInstance)

    def test_setattr(self, NewInstance):
        assert isroutine(NewInstance.__setattr__)
        # try to set order: shouldn't change anything
        NewInstance.order = 0
        assert NewInstance.reason is None
        assert NewInstance.order is None
        # set all reasons in decreasing order and compare it
        for r in reversed(totest.IGNORE_REASONS_PRIORITY):
            o = totest.IGNORE_REASONS_PRIORITY.index(r)
            NewInstance.reason = r
            assert NewInstance.order == o
            assert NewInstance.reason == r
        # try to set reason of lowest priority: reason of higher pririty will
        # be kept
        NewInstance.reason = totest.IGNORE_REASONS_PRIORITY[-1]
        assert NewInstance.order == o
        assert NewInstance.reason == r
        # unset the reason: set back to None
        NewInstance.reason = None
        assert NewInstance.reason is None
        assert NewInstance.order is None
        # try error on non existing reason
        assert '' not in totest.IGNORE_REASONS_PRIORITY
        with raises(ValueError):
            NewInstance.reason = ''
