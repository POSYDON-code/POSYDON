"""Unit tests of posydon/utils/posydonerror.py
"""

__authors__ = [
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>"
]

# import the unittest module and the module which will be tested
import unittest
import posydon.utils.posydonerror as totest

# import other needed code for the tests
from inspect import isclass, isroutine

# define test classes
class TestElements(unittest.TestCase):
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = ['BinaryStar', 'FlowError', 'GridError', 'MatchingError',
                    'ModelError', 'NumericalError', 'POSYDONError',
                    'SingleStar', '__authors__', '__builtins__', '__cached__',
                    '__doc__', '__file__', '__loader__', '__name__',
                    '__package__', '__spec__', 'copy',
                    'initial_condition_message']
        self.assertListEqual(dir(totest), elements,
                             msg="There might be added or removed objects "
                                 "without an update on the unit test.")

    def test_instance_POSYDONError(self):
        self.assertTrue(isclass(totest.POSYDONError))
        self.assertTrue(issubclass(totest.POSYDONError, Exception))

    def test_instance_FlowError(self):
        self.assertTrue(isclass(totest.FlowError))
        self.assertTrue(issubclass(totest.FlowError, totest.POSYDONError))

    def test_instance_GridError(self):
        self.assertTrue(isclass(totest.GridError))
        self.assertTrue(issubclass(totest.GridError, totest.POSYDONError))

    def test_instance_MatchingError(self):
        self.assertTrue(isclass(totest.MatchingError))
        self.assertTrue(issubclass(totest.MatchingError, totest.POSYDONError))

    def test_instance_ModelError(self):
        self.assertTrue(isclass(totest.ModelError))
        self.assertTrue(issubclass(totest.ModelError, totest.POSYDONError))

    def test_instance_NumericalError(self):
        self.assertTrue(isclass(totest.NumericalError))
        self.assertTrue(issubclass(totest.NumericalError, totest.POSYDONError))

    def test_instance_initial_condition_message(self):
        self.assertTrue(isroutine(totest.initial_condition_message))


class TestFunctions(unittest.TestCase):
    # test functions
    def setUp(self):
        # initialize a BinaryStar instance, which is a required argument
        self.BinaryStar = totest.BinaryStar()

    def test_initial_condition_message(self):
        test_object = {'Test': 'object'}
        with self.assertRaises(TypeError):
            message = totest.initial_condition_message(binary=test_object)
        self.assertIn("Failed Binary Initial Conditions",
            totest.initial_condition_message(binary=self.BinaryStar))
        with self.assertRaises(TypeError):
            message = totest.initial_condition_message(binary=test_object,
                                                       ini_params=1)
        with self.assertRaises(TypeError):
            message = totest.initial_condition_message(binary=test_object,
                                                       ini_params=[1,2])
        self.assertEqual("a: 1\nb: 2\n",
            totest.initial_condition_message(binary=self.BinaryStar,
                                             ini_params=["a: 1\n", "b: 2\n"]))


class TestPOSYDONError(unittest.TestCase):
    # test the POSYDONError class
    def setUp(self):
        # initialize an instance of the class for each test
        self.POSYDONError = totest.POSYDONError("test message on posittion")

    def test_init(self):
        self.assertTrue(isroutine(self.POSYDONError.__init__))
        # check that the instance is of correct type and all code in the
        # __init__ got executed: the elements are created and initialized
        self.assertIsInstance(self.POSYDONError, totest.POSYDONError)
        self.assertEqual(self.POSYDONError.message,
                         "test message on posittion")
        self.assertIsNone(self.POSYDONError.objects)
        # check defaults
        error_object = totest.POSYDONError()
        self.assertEqual(error_object.message, "")
        self.assertIsNone(error_object.objects)
        # test a passed message
        error_object = totest.POSYDONError(message="test message with key")
        self.assertEqual(error_object.message, "test message with key")
        self.assertIsNone(error_object.objects)
        # test a passed object
        test_object = {'Test': 'object'}
        error_object = totest.POSYDONError(objects=[test_object])
        self.assertEqual(error_object.message, "")
        self.assertListEqual(error_object.objects, [test_object])
        # test requests on input parameters
        with self.assertRaises(TypeError):
            error_object = totest.POSYDONError(message=test_object)
        with self.assertRaises(TypeError):
            error_object = totest.POSYDONError(objects=test_object)


    def test_str(self):
        self.assertTrue(isroutine(self.POSYDONError.__str__))
        self.assertEqual(str(self.POSYDONError), "\ntest message on posittion")
        # test passed objects
        test_object1 = {'Test': 'object'}
        error_object = totest.POSYDONError(objects=[test_object1])
        self.assertEqual(str(error_object), "\n")
        test_object2 = totest.SingleStar()
        error_object = totest.POSYDONError(objects=test_object2)
        self.assertIn("OBJECT #(<class 'posydon.binary_evol.singlestar.SingleStar'>)", str(error_object))
        test_object3 = totest.BinaryStar()
        error_object = totest.POSYDONError(objects=test_object3)
        self.assertIn("OBJECT #(<class 'posydon.binary_evol.binarystar.BinaryStar'>)", str(error_object))
        error_object = totest.POSYDONError(objects=[test_object1, test_object2, test_object3])
        self.assertNotIn("OBJECT #1", str(error_object))
        self.assertIn("OBJECT #2 (<class 'posydon.binary_evol.singlestar.SingleStar'>)", str(error_object))
        self.assertIn("OBJECT #3 (<class 'posydon.binary_evol.binarystar.BinaryStar'>)", str(error_object))


if __name__ == "__main__":
    unittest.main()
