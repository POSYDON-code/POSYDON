"""Unit tests of posydon/utils/ignorereason.py
"""

__authors__ = [
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>"
]

# import the unittest module and the module which will be tested
import unittest
import posydon.utils.ignorereason as totest

# import other needed code for the tests
from inspect import isclass, isroutine

# define test classes
class TestElements(unittest.TestCase):
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = ['IGNORE_REASONS_PRIORITY', 'IgnoreReason', '__authors__',
                    '__builtins__', '__cached__', '__doc__', '__file__',
                    '__loader__', '__name__', '__package__', '__spec__']
        self.assertListEqual(dir(totest), elements,
                             msg="There might be added or removed objects "
                                 "without an update on the unit test.")

    def test_instance_IGNORE_REASONS_PRIORITY(self):
        self.assertIsInstance(totest.IGNORE_REASONS_PRIORITY, (list))

    def test_instance_IgnoreReason(self):
        self.assertTrue(isclass(totest.IgnoreReason))


class TestValues(unittest.TestCase):
    # check that the values fit
    def test_value_IGNORE_REASONS_PRIORITY(self):
        self.assertIn('ignored_no_history1', totest.IGNORE_REASONS_PRIORITY)
        self.assertIn('ignored_no_binary_history',
                      totest.IGNORE_REASONS_PRIORITY)
        self.assertIn('corrupted_history1', totest.IGNORE_REASONS_PRIORITY)
        self.assertIn('corrupted_binary_history',
                      totest.IGNORE_REASONS_PRIORITY)
        self.assertIn('corrupted_history2', totest.IGNORE_REASONS_PRIORITY)
        self.assertIn('ignored_scrubbed_history',
                      totest.IGNORE_REASONS_PRIORITY)
        self.assertIn('ignored_no_final_profile',
                      totest.IGNORE_REASONS_PRIORITY)
        self.assertIn('ignored_no_RLO', totest.IGNORE_REASONS_PRIORITY)


class TestIgnoreReason(unittest.TestCase):
    # test the IgnoreReason class
    def setUp(self):
        # initialize an instance of the class for each test
        self.IgnoreReason = totest.IgnoreReason()

    def test_init(self):
        self.assertTrue(isroutine(self.IgnoreReason.__init__))
        # check that the instance is of correct type and all code in the
        # __init__ got executed: the elements are created and initialized
        self.assertIsInstance(self.IgnoreReason, totest.IgnoreReason)
        self.assertIsNone(self.IgnoreReason.reason)
        self.assertIsNone(self.IgnoreReason.order)

    def test_bool(self):
        self.assertTrue(isroutine(self.IgnoreReason.__bool__))
        self.assertFalse(self.IgnoreReason)
        self.IgnoreReason.reason = totest.IGNORE_REASONS_PRIORITY[0]
        self.assertTrue(self.IgnoreReason)

    def test_setattr(self):
        self.assertTrue(isroutine(self.IgnoreReason.__setattr__))
        # try to set order: shouldn't change anything
        self.IgnoreReason.order = 0
        self.assertIsNone(self.IgnoreReason.reason)
        self.assertIsNone(self.IgnoreReason.order)
        # set all reasons in decreasing order and compare it
        for r in reversed(totest.IGNORE_REASONS_PRIORITY):
            o = totest.IGNORE_REASONS_PRIORITY.index(r)
            self.IgnoreReason.reason = r
            self.assertEqual(self.IgnoreReason.order, o)
            self.assertEqual(self.IgnoreReason.reason, r)
        # try to set reason of lowest priority: reason of higher pririty will
        # be kept
        self.IgnoreReason.reason = totest.IGNORE_REASONS_PRIORITY[-1]
        self.assertEqual(self.IgnoreReason.order, o)
        self.assertEqual(self.IgnoreReason.reason, r)
        # unset the reason: set back to None
        self.IgnoreReason.reason = None
        self.assertIsNone(self.IgnoreReason.reason)
        self.assertIsNone(self.IgnoreReason.order)
        # try error on non existing reason
        self.assertNotIn('', totest.IGNORE_REASONS_PRIORITY)
        with self.assertRaises(ValueError):
            self.IgnoreReason.reason = ''


if __name__ == "__main__":
    unittest.main()
