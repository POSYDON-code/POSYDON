"""Unit tests of posydon/utils/posydonwarning.py
"""

__authors__ = [
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>"
]

# import the unittest module and the module which will be tested
import unittest
import posydon.utils.posydonwarning as totest

# import other needed code for the tests
from inspect import isclass, isroutine
from io import StringIO
from contextlib import redirect_stdout, redirect_stderr

# define test classes
class TestElements(unittest.TestCase):
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = ['AllPOSYDONWarnings', 'ApproximationWarning',
                    'BinaryParsingWarning', 'Catch_POSYDON_Warnings',
                    'ClassificationWarning', 'EvolutionWarning',
                    'InappropriateValueWarning', 'IncompletenessWarning',
                    'InterpolationWarning', 'MissingFilesWarning',
                    'NoPOSYDONWarnings', 'OverwriteWarning', 'POSYDONWarning',
                    'Pwarn', 'ReplaceValueWarning', 'SetPOSYDONWarnings',
                    'UnsupportedModelWarning', '_CAUGHT_POSYDON_WARNINGS',
                    '_Caught_POSYDON_Warnings', '_POSYDONWarning_subclasses',
                    '_POSYDON_WARNINGS_REGISTRY', '__authors__',
                    '__builtins__', '__cached__', '__doc__', '__file__',
                    '__loader__', '__name__', '__package__', '__spec__',
                    '_apply_POSYDON_filter', '_get_POSYDONWarning_class',
                    '_issue_warn', 'copy', 'get_stats', 'print_stats', 'sys',
                    'warnings']
        self.assertListEqual(dir(totest), elements,
                             msg="There might be added or removed objects "
                                 "without an update on the unit test.")

    def test_instance_POSYDONWarning(self):
        self.assertTrue(isclass(totest.POSYDONWarning))
        self.assertTrue(issubclass(totest.POSYDONWarning, Warning))

    def test_instance_ApproximationWarning(self):
        self.assertTrue(isclass(totest.ApproximationWarning))
        self.assertTrue(issubclass(totest.ApproximationWarning,
                                   totest.POSYDONWarning))

    def test_instance_BinaryParsingWarning(self):
        self.assertTrue(isclass(totest.BinaryParsingWarning))
        self.assertTrue(issubclass(totest.BinaryParsingWarning,
                                   totest.POSYDONWarning))

    def test_instance_ClassificationWarning(self):
        self.assertTrue(isclass(totest.ClassificationWarning))
        self.assertTrue(issubclass(totest.ClassificationWarning,
                                   totest.POSYDONWarning))

    def test_instance_EvolutionWarning(self):
        self.assertTrue(isclass(totest.EvolutionWarning))
        self.assertTrue(issubclass(totest.EvolutionWarning,
                                   totest.POSYDONWarning))

    def test_instance_InappropriateValueWarning(self):
        self.assertTrue(isclass(totest.InappropriateValueWarning))
        self.assertTrue(issubclass(totest.InappropriateValueWarning,
                                   totest.POSYDONWarning))

    def test_instance_IncompletenessWarning(self):
        self.assertTrue(isclass(totest.IncompletenessWarning))
        self.assertTrue(issubclass(totest.IncompletenessWarning,
                                   totest.POSYDONWarning))

    def test_instance_InterpolationWarning(self):
        self.assertTrue(isclass(totest.InterpolationWarning))
        self.assertTrue(issubclass(totest.InterpolationWarning,
                                   totest.POSYDONWarning))

    def test_instance_MissingFilesWarning(self):
        self.assertTrue(isclass(totest.MissingFilesWarning))
        self.assertTrue(issubclass(totest.MissingFilesWarning,
                                   totest.POSYDONWarning))

    def test_instance_OverwriteWarning(self):
        self.assertTrue(isclass(totest.OverwriteWarning))
        self.assertTrue(issubclass(totest.OverwriteWarning,
                                   totest.POSYDONWarning))

    def test_instance_ReplaceValueWarning(self):
        self.assertTrue(isclass(totest.ReplaceValueWarning))
        self.assertTrue(issubclass(totest.ReplaceValueWarning,
                                   totest.POSYDONWarning))

    def test_instance_UnsupportedModelWarning(self):
        self.assertTrue(isclass(totest.UnsupportedModelWarning))
        self.assertTrue(issubclass(totest.UnsupportedModelWarning,
                                   totest.POSYDONWarning))

    def test_instance_POSYDONWarning_subclasses(self):
        self.assertIsInstance(totest._POSYDONWarning_subclasses, (dict))

    def test_instance_get_POSYDONWarning_class(self):
        self.assertTrue(isroutine(totest._get_POSYDONWarning_class))

    def test_instance_POSYDON_WARNINGS_REGISTRY(self):
        self.assertIsInstance(totest._POSYDON_WARNINGS_REGISTRY, (dict))

    def test_instance_get_stats(self):
        self.assertTrue(isroutine(totest.get_stats))

    def test_instance_print_stats(self):
        self.assertTrue(isroutine(totest.print_stats))

    def test_instance_apply_POSYDON_filter(self):
        self.assertTrue(isroutine(totest._apply_POSYDON_filter))

    def test_instance_issue_warn(self):
        self.assertTrue(isroutine(totest._issue_warn))

    def test_instance_Caught_POSYDON_Warnings(self):
        self.assertTrue(isclass(totest._Caught_POSYDON_Warnings))

    def test_instance_CAUGHT_POSYDON_WARNINGS(self):
        self.assertIsInstance(totest._CAUGHT_POSYDON_WARNINGS,
                              totest._Caught_POSYDON_Warnings)

    def test_instance_Catch_POSYDON_Warnings(self):
        self.assertTrue(isclass(totest.Catch_POSYDON_Warnings))

    def test_instance_Pwarn(self):
        self.assertTrue(isroutine(totest.Pwarn))

    def test_instance_SetPOSYDONWarnings(self):
        self.assertTrue(isroutine(totest.SetPOSYDONWarnings))

    def test_instance_NoPOSYDONWarnings(self):
        self.assertTrue(isroutine(totest.NoPOSYDONWarnings))

    def test_instance_AllPOSYDONWarnings(self):
        self.assertTrue(isroutine(totest.AllPOSYDONWarnings))


class TestValues(unittest.TestCase):
    # check that the values fit
    def test_value_POSYDONWarning_subclasses(self):
        self.assertIn('ApproximationWarning',
                      totest._POSYDONWarning_subclasses)

    def test_value_POSYDON_WARNINGS_REGISTRY(self):
        self.assertDictEqual({}, totest._POSYDON_WARNINGS_REGISTRY)

    def test_value_CAUGHT_POSYDON_WARNINGS(self):
        self.assertFalse(totest._CAUGHT_POSYDON_WARNINGS.catch_warnings)
        self.assertTrue(totest._CAUGHT_POSYDON_WARNINGS.record)
        self.assertTrue(totest._CAUGHT_POSYDON_WARNINGS.filter_first)
        self.assertListEqual(totest._CAUGHT_POSYDON_WARNINGS.caught_warnings,
                             [])
        self.assertEqual(totest._CAUGHT_POSYDON_WARNINGS.registry,
                         totest._POSYDON_WARNINGS_REGISTRY)


class TestFunctions(unittest.TestCase):
    def tearDown(self):
        # empyt the global POSYDON warnings registry after each test
        keys = []
        for k in totest._POSYDON_WARNINGS_REGISTRY:
            keys.append(k)
        for k in keys:
            del totest._POSYDON_WARNINGS_REGISTRY[k]
        # set POSYDON warnings back to default
        totest.SetPOSYDONWarnings()

    # test functions
    def test_get_POSYDONWarning_class(self):
        # missing argument
        with self.assertRaises(TypeError):
            totest._get_POSYDONWarning_class()
        # default to POSYDONWarning
        self.assertEqual(totest._get_POSYDONWarning_class(""),
                         totest.POSYDONWarning)
        self.assertEqual(totest._get_POSYDONWarning_class(
            totest.POSYDONWarning), totest.POSYDONWarning)
        # check subclasses of POSYDONWarning
        for k,v in totest._POSYDONWarning_subclasses.items():
            self.assertEqual(totest._get_POSYDONWarning_class(k), v)
            self.assertEqual(totest._get_POSYDONWarning_class(v), v)
        # bad input
        self.assertIsNone(totest._get_POSYDONWarning_class(1))

    def test_get_stats(self):
        self.assertEqual(totest.get_stats(), totest._POSYDON_WARNINGS_REGISTRY)

    def test_print_stats(self):
        # no warnings to print
        with redirect_stdout(StringIO()) as print_out:
            totest.print_stats()
        self.assertEqual("No POSYDON warnings occured.\n",
                         print_out.getvalue())
        # add an artifical entry in the warnings registry and get the printout
        totest._POSYDON_WARNINGS_REGISTRY = {'Unit': 'Test'}
        with redirect_stdout(StringIO()) as print_out:
            totest.print_stats()
        self.assertEqual("There have been POSYDON warnings in the global "+
                         "registry:\n "+str(totest._POSYDON_WARNINGS_REGISTRY)+
                         "\n", print_out.getvalue())

    def test_apply_POSYDON_filter(self):
        # wrong arguments
        with self.assertRaises(TypeError):
            totest._apply_POSYDON_filter(warning="Test")
        with self.assertRaises(TypeError):
            totest._apply_POSYDON_filter(warning={'message': 1})
        with self.assertRaises(TypeError):
            totest._apply_POSYDON_filter(warning={'stacklevel': "Test"})
        with redirect_stdout(StringIO()) as print_out:
            totest._apply_POSYDON_filter(registry="Test")
        self.assertEqual("Reset registry, old was: Test\n",
                         print_out.getvalue())
        # check default warning
        self.assertDictEqual(dict(message="No warning"),
                             totest._apply_POSYDON_filter())
        # check that further default warnings are filtered out but added to the
        # registry
        for i in range(10):
            for k,v in totest._POSYDON_WARNINGS_REGISTRY.items():
                self.assertEqual(i+1, v)
            self.assertIsNone(totest._apply_POSYDON_filter())

    def test_issue_warn(self):
        # wrong arguments
        with self.assertRaises(TypeError):
            totest._issue_warn(warning="Test")
        with self.assertRaises(TypeError):
            totest._issue_warn(warning={'message': 1})
        with self.assertRaises(TypeError):
            totest._issue_warn(warning={'stacklevel': "Test"})
        with redirect_stdout(StringIO()) as print_out:
            # it will issue a default warning on the reset registry
            with redirect_stderr(StringIO()) as print_err:
                totest._issue_warn(registry="Test")
        self.assertEqual("Reset registry, old was: Test\n",
                         print_out.getvalue())
        # check default warning
        self.assertIn("UserWarning: No warning", print_err.getvalue())
        # check filtered warning
        self.assertIsNone(totest._issue_warn())

    def test_Pwarn(self):
        # wrong arguments
        with self.assertRaises(TypeError):
            totest.Pwarn()
        with self.assertRaises(TypeError):
            totest.Pwarn(1)
        with self.assertRaises(TypeError):
            totest.Pwarn("Unit", stacklevel="Test")
        # check output of POSYDONwarning and UserWarning
        with redirect_stderr(StringIO()) as print_err:
            totest.Pwarn("Unit test", "POSYDONWarning")
        self.assertIn("POSYDONWarning: 'Unit test'", print_err.getvalue())
        with redirect_stderr(StringIO()) as print_err:
            totest.Pwarn("Unit test")
        self.assertIn("UserWarning: Unit test", print_err.getvalue())

    def test_SetPOSYDONWarnings(self):
        totest.SetPOSYDONWarnings()
        self.assertEqual("('default', None, <class "
                         "'posydon.utils.posydonwarning.POSYDONWarning'>, "
                         "None, 0)", str(totest.warnings.filters[0]))

    def test_NoPOSYDONWarnings(self):
        totest.NoPOSYDONWarnings()
        self.assertEqual("('ignore', None, <class "
                         "'posydon.utils.posydonwarning.POSYDONWarning'>, "
                         "None, 0)", str(totest.warnings.filters[0]))

    def test_AllPOSYDONWarnings(self):
        totest.AllPOSYDONWarnings()
        self.assertEqual("('always', None, <class "
                         "'posydon.utils.posydonwarning.POSYDONWarning'>, "
                         "None, 0)", str(totest.warnings.filters[0]))


class TestPOSYDONWarning(unittest.TestCase):
    # test the POSYDONWarning class
    def setUp(self):
        # initialize an instance of the class for each test
        self.POSYDONWarning = totest.POSYDONWarning()

    def test_init(self):
        self.assertTrue(isroutine(self.POSYDONWarning.__init__))
        # check that the instance is of correct type and all code in the
        # __init__ got executed: the elements are created and initialized
        self.assertIsInstance(self.POSYDONWarning, totest.POSYDONWarning)
        self.assertEqual('', self.POSYDONWarning.message)

    def test_str(self):
        self.assertTrue(isroutine(self.POSYDONWarning.__str__))
        self.assertEqual("''", str(self.POSYDONWarning))


class TestApproximationWarning(unittest.TestCase):
    # test the ApproximationWarning class
    def setUp(self):
        # initialize an instance of the class for each test
        self.ApproximationWarning = totest.ApproximationWarning()

    def test_init(self):
        self.assertTrue(isroutine(self.ApproximationWarning.__init__))
        # check that the instance is of correct type and all code in the
        # __init__ got executed: the elements are created and initialized
        self.assertIsInstance(self.ApproximationWarning,
                              totest.ApproximationWarning)
        self.assertEqual('', self.ApproximationWarning.message)


class TestBinaryParsingWarning(unittest.TestCase):
    # test the BinaryParsingWarning class
    def setUp(self):
        # initialize an instance of the class for each test
        self.BinaryParsingWarning = totest.BinaryParsingWarning()

    def test_init(self):
        self.assertTrue(isroutine(self.BinaryParsingWarning.__init__))
        # check that the instance is of correct type and all code in the
        # __init__ got executed: the elements are created and initialized
        self.assertIsInstance(self.BinaryParsingWarning,
                              totest.BinaryParsingWarning)
        self.assertEqual('', self.BinaryParsingWarning.message)


class TestClassificationWarning(unittest.TestCase):
    # test the ClassificationWarning class
    def setUp(self):
        # initialize an instance of the class for each test
        self.ClassificationWarning = totest.ClassificationWarning()

    def test_init(self):
        self.assertTrue(isroutine(self.ClassificationWarning.__init__))
        # check that the instance is of correct type and all code in the
        # __init__ got executed: the elements are created and initialized
        self.assertIsInstance(self.ClassificationWarning,
                              totest.ClassificationWarning)
        self.assertEqual('', self.ClassificationWarning.message)


class TestEvolutionWarning(unittest.TestCase):
    # test the EvolutionWarning class
    def setUp(self):
        # initialize an instance of the class for each test
        self.EvolutionWarning = totest.EvolutionWarning()

    def test_init(self):
        self.assertTrue(isroutine(self.EvolutionWarning.__init__))
        # check that the instance is of correct type and all code in the
        # __init__ got executed: the elements are created and initialized
        self.assertIsInstance(self.EvolutionWarning, totest.EvolutionWarning)
        self.assertEqual('', self.EvolutionWarning.message)


class TestInappropriateValueWarning(unittest.TestCase):
    # test the InappropriateValueWarning class
    def setUp(self):
        # initialize an instance of the class for each test
        self.InappropriateValueWarning = totest.InappropriateValueWarning()

    def test_init(self):
        self.assertTrue(isroutine(self.InappropriateValueWarning.__init__))
        # check that the instance is of correct type and all code in the
        # __init__ got executed: the elements are created and initialized
        self.assertIsInstance(self.InappropriateValueWarning,
                              totest.InappropriateValueWarning)
        self.assertEqual('', self.InappropriateValueWarning.message)


class TestIncompletenessWarning(unittest.TestCase):
    # test the IncompletenessWarning class
    def setUp(self):
        # initialize an instance of the class for each test
        self.IncompletenessWarning = totest.IncompletenessWarning()

    def test_init(self):
        self.assertTrue(isroutine(self.IncompletenessWarning.__init__))
        # check that the instance is of correct type and all code in the
        # __init__ got executed: the elements are created and initialized
        self.assertIsInstance(self.IncompletenessWarning,
                              totest.IncompletenessWarning)
        self.assertEqual('', self.IncompletenessWarning.message)


class TestInterpolationWarning(unittest.TestCase):
    # test the InterpolationWarning class
    def setUp(self):
        # initialize an instance of the class for each test
        self.InterpolationWarning = totest.InterpolationWarning()

    def test_init(self):
        self.assertTrue(isroutine(self.InterpolationWarning.__init__))
        # check that the instance is of correct type and all code in the
        # __init__ got executed: the elements are created and initialized
        self.assertIsInstance(self.InterpolationWarning,
                              totest.InterpolationWarning)
        self.assertEqual('', self.InterpolationWarning.message)


class TestMissingFilesWarning(unittest.TestCase):
    # test the MissingFilesWarning class
    def setUp(self):
        # initialize an instance of the class for each test
        self.MissingFilesWarning = totest.MissingFilesWarning()

    def test_init(self):
        self.assertTrue(isroutine(self.MissingFilesWarning.__init__))
        # check that the instance is of correct type and all code in the
        # __init__ got executed: the elements are created and initialized
        self.assertIsInstance(self.MissingFilesWarning,
                              totest.MissingFilesWarning)
        self.assertEqual('', self.MissingFilesWarning.message)


class TestOverwriteWarning(unittest.TestCase):
    # test the OverwriteWarning class
    def setUp(self):
        # initialize an instance of the class for each test
        self.OverwriteWarning = totest.OverwriteWarning()

    def test_init(self):
        self.assertTrue(isroutine(self.OverwriteWarning.__init__))
        # check that the instance is of correct type and all code in the
        # __init__ got executed: the elements are created and initialized
        self.assertIsInstance(self.OverwriteWarning,
                              totest.OverwriteWarning)
        self.assertEqual('', self.OverwriteWarning.message)


class TestReplaceValueWarning(unittest.TestCase):
    # test the ReplaceValueWarning class
    def setUp(self):
        # initialize an instance of the class for each test
        self.ReplaceValueWarning = totest.ReplaceValueWarning()

    def test_init(self):
        self.assertTrue(isroutine(self.ReplaceValueWarning.__init__))
        # check that the instance is of correct type and all code in the
        # __init__ got executed: the elements are created and initialized
        self.assertIsInstance(self.ReplaceValueWarning,
                              totest.ReplaceValueWarning)
        self.assertEqual('', self.ReplaceValueWarning.message)


class TestUnsupportedModelWarning(unittest.TestCase):
    # test the UnsupportedModelWarning class
    def setUp(self):
        # initialize an instance of the class for each test
        self.UnsupportedModelWarning = totest.UnsupportedModelWarning()

    def test_init(self):
        self.assertTrue(isroutine(self.UnsupportedModelWarning.__init__))
        # check that the instance is of correct type and all code in the
        # __init__ got executed: the elements are created and initialized
        self.assertIsInstance(self.UnsupportedModelWarning,
                              totest.UnsupportedModelWarning)
        self.assertEqual('', self.UnsupportedModelWarning.message)


class Test_Caught_POSYDON_Warnings(unittest.TestCase):
    # test the _Caught_POSYDON_Warnings class
    def setUp(self):
        # initialize an instance of the class for each test
        self._Caught_POSYDON_Warnings = totest._Caught_POSYDON_Warnings()

    def tearDown(self):
        # empty the cache to not print remaining records
        self._Caught_POSYDON_Warnings(empty_cache=True)
        # empyt the global POSYDON warnings registry after each test
        keys = []
        for k in totest._POSYDON_WARNINGS_REGISTRY:
            keys.append(k)
        for k in keys:
            del totest._POSYDON_WARNINGS_REGISTRY[k]

    def test_init(self):
        self.assertTrue(isroutine(self._Caught_POSYDON_Warnings.__init__))
        # check that the instance is of correct type and all code in the
        # __init__ got executed: the elements are created and initialized
        self.assertIsInstance(self._Caught_POSYDON_Warnings,
                              totest._Caught_POSYDON_Warnings)
        self.assertFalse(self._Caught_POSYDON_Warnings.catch_warnings)
        self.assertListEqual([], self._Caught_POSYDON_Warnings.caught_warnings)
        self.assertTrue(self._Caught_POSYDON_Warnings.record)
        self.assertTrue(self._Caught_POSYDON_Warnings.filter_first)
        self.assertFalse(self._Caught_POSYDON_Warnings._got_called)
        self.assertDictEqual(totest._POSYDON_WARNINGS_REGISTRY,
                             self._Caught_POSYDON_Warnings.registry)
        # bad input
        with self.assertRaises(TypeError):
            totest._Caught_POSYDON_Warnings(catch_warnings="Test")
        with self.assertRaises(TypeError):
            totest._Caught_POSYDON_Warnings(record="Test")
        with self.assertRaises(TypeError):
            totest._Caught_POSYDON_Warnings(filter_first="Test")
        with self.assertRaises(TypeError):
            totest._Caught_POSYDON_Warnings(registry="Test")

    def test_str(self):
        self.assertTrue(isroutine(self._Caught_POSYDON_Warnings.__str__))
        self.assertEqual("POSYDON warnings are shown.",
                         str(self._Caught_POSYDON_Warnings))
        # check with catching
        caught_object = totest._Caught_POSYDON_Warnings(catch_warnings=True)
        self.assertEqual("POSYDON warnings will be caught and recorded. "
                         "Filters are applied before recording.",
                         str(caught_object))
        # check without recording and own registry
        test_registry = {'Unit': "Test"}
        caught_object = totest._Caught_POSYDON_Warnings(catch_warnings=True,
                                                        record=False,
                                                        registry=test_registry)
        self.assertEqual("POSYDON warnings will be caught and discarded. "
                         "Currently a private registry is used, it contains:\n"
                         "{'Unit': 'Test'}", str(caught_object))

    def test_call(self):
        self.assertTrue(isroutine(self._Caught_POSYDON_Warnings.__call__))
        # bad input
        with self.assertRaises(ValueError):
            self._Caught_POSYDON_Warnings()
        self.assertTrue(self._Caught_POSYDON_Warnings._got_called)
        with self.assertRaises(TypeError):
            self._Caught_POSYDON_Warnings(change_settings="Test")
        with self.assertRaises(TypeError):
            self._Caught_POSYDON_Warnings(change_settings={'catch_warnings':
                                                           "Test"})
        with self.assertRaises(AttributeError):
            self._Caught_POSYDON_Warnings(change_settings={'Unit': "Test"})
        # change setting to catch warnings and add a new one to the record list
        self._Caught_POSYDON_Warnings(change_settings={'catch_warnings': True},
                                      new_warning={'message': "Test"})
        self.assertTrue(self._Caught_POSYDON_Warnings.catch_warnings)
        self.assertListEqual([{'message': "Test"}],
                             self._Caught_POSYDON_Warnings.caught_warnings)

    def test_del(self):
        self.assertTrue(isroutine(self._Caught_POSYDON_Warnings.__del__))
        # create an object with a catched warning, call the destructor and
        # empty it while recording the stdout and stderr
        with redirect_stderr(StringIO()) as print_err:
            caught_object = totest._Caught_POSYDON_Warnings(
                                catch_warnings=True)
            caught_object(new_warning={'message': "Unit Test"})
            caught_object.__del__()
            caught_object(empty_cache=True)
        self.assertIn("There are still recorded warnings:",
                      print_err.getvalue())
        self.assertIn("UserWarning: Unit Test", print_err.getvalue())

    def test_got_called(self):
        self.assertTrue(isroutine(self._Caught_POSYDON_Warnings.got_called))
        self.assertFalse(self._Caught_POSYDON_Warnings.got_called())
        self._Caught_POSYDON_Warnings(empty_cache=True)
        self.assertTrue(self._Caught_POSYDON_Warnings.got_called())

    def test_has_records(self):
        self.assertTrue(isroutine(self._Caught_POSYDON_Warnings.has_records))
        self.assertFalse(self._Caught_POSYDON_Warnings.has_records())
        self._Caught_POSYDON_Warnings(change_settings={'catch_warnings': True},
                                      new_warning={'message': "Unit Test"})
        self.assertTrue(self._Caught_POSYDON_Warnings.has_records())

    def test_get_cache(self):
        self.assertTrue(isroutine(self._Caught_POSYDON_Warnings.get_cache))
        self.assertListEqual([], self._Caught_POSYDON_Warnings.get_cache())

    def test_reset_cache(self):
        self.assertTrue(isroutine(self._Caught_POSYDON_Warnings.reset_cache))
        self._Caught_POSYDON_Warnings(change_settings={'catch_warnings': True},
                                      new_warning={'message': "Unit Test"})
        self.assertListEqual([{'message': "Unit Test"}],
                             self._Caught_POSYDON_Warnings.caught_warnings)
        self._Caught_POSYDON_Warnings.reset_cache()
        self.assertListEqual([], self._Caught_POSYDON_Warnings.caught_warnings)


class TestCatch_POSYDON_Warnings(unittest.TestCase):
    # test the Catch_POSYDON_Warnings class
    def setUp(self):
        # initialize an instance of the class for each test
        self.Catch_POSYDON_Warnings = totest.Catch_POSYDON_Warnings()

    def tearDown(self):
        # empyt the global POSYDON warnings registry after each test
        keys = []
        for k in totest._POSYDON_WARNINGS_REGISTRY:
            keys.append(k)
        for k in keys:
            del totest._POSYDON_WARNINGS_REGISTRY[k]

    def test_init(self):
        self.assertTrue(isroutine(self.Catch_POSYDON_Warnings.__init__))
        # check that the instance is of correct type and all code in the
        # __init__ got executed: the elements are created and initialized
        self.assertTrue(self.Catch_POSYDON_Warnings.catch_warnings)
        self.assertTrue(self.Catch_POSYDON_Warnings.record)
        self.assertTrue(self.Catch_POSYDON_Warnings.filter_first)
        self.assertIsNone(self.Catch_POSYDON_Warnings.context_registry)
        self.assertIsNone(self.Catch_POSYDON_Warnings.python_catch)
        # check other inputs
        catch_object = totest.Catch_POSYDON_Warnings(own_registry=True,
                                                     use_python_catch=True)
        self.assertDictEqual({}, catch_object.context_registry)
        self.assertIsInstance(catch_object.python_catch,
                            totest.warnings.catch_warnings)

    def test_enter_exit(self):
        self.assertTrue(isroutine(self.Catch_POSYDON_Warnings.__enter__))

    def test_exit(self):
        self.assertTrue(isroutine(self.Catch_POSYDON_Warnings.__exit__))

    def test_context(self):
        with totest.Catch_POSYDON_Warnings() as cpw:
            self.assertEqual(totest._CAUGHT_POSYDON_WARNINGS, cpw)


if __name__ == "__main__":
    unittest.main()
