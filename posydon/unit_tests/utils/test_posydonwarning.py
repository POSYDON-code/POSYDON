"""Unit tests of posydon/utils/posydonwarning.py
"""

__authors__ = [
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>"
]

# ensure that python forgets about previous imports of posydonwarning, e.g. in
# other tests, before importing it here
from sys import modules as sys_modules
sys_modules.pop('posydon.utils.posydonwarning', None)

# import the module which will be tested
import posydon.utils.posydonwarning as totest

# import other needed code for the tests, which is not already imported in the
# module you like to test
from pytest import fixture, raises, warns
from inspect import isclass, isroutine

# define test classes collecting several test functions
class TestElements:
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = {'AllPOSYDONWarnings', 'ApproximationWarning',\
                    'BinaryParsingWarning', 'Catch_POSYDON_Warnings',\
                    'ClassificationWarning', 'EvolutionWarning',\
                    'InappropriateValueWarning', 'IncompletenessWarning',\
                    'InterpolationWarning', 'MissingFilesWarning',\
                    'NoPOSYDONWarnings', 'OverwriteWarning', 'SFHModelWarning',\
                    'POSYDONWarning','Pwarn', 'ReplaceValueWarning',\
                    'SetPOSYDONWarnings','DeprecationWarning', \
                    'UnsupportedModelWarning', '_CAUGHT_POSYDON_WARNINGS',\
                    '_Caught_POSYDON_Warnings', '_POSYDONWarning_subclasses',\
                    '_POSYDON_WARNINGS_REGISTRY', '__authors__',\
                    '__builtins__', '__cached__', '__doc__', '__file__',\
                    '__loader__', '__name__', '__package__', '__spec__',\
                    '_apply_POSYDON_filter', '_get_POSYDONWarning_class',\
                    '_issue_warn', 'copy', 'get_stats', 'print_stats', 'sys',\
                    'warnings', 'nosrc_code_format'}
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

    def test_instance_POSYDONWarning(self):
        assert isclass(totest.POSYDONWarning)
        assert issubclass(totest.POSYDONWarning, Warning)

    def test_instance_ApproximationWarning(self):
        assert isclass(totest.ApproximationWarning)
        assert issubclass(totest.ApproximationWarning, totest.POSYDONWarning)

    def test_instance_BinaryParsingWarning(self):
        assert isclass(totest.BinaryParsingWarning)
        assert issubclass(totest.BinaryParsingWarning, totest.POSYDONWarning)

    def test_instance_ClassificationWarning(self):
        assert isclass(totest.ClassificationWarning)
        assert issubclass(totest.ClassificationWarning, totest.POSYDONWarning)

    def test_instance_EvolutionWarning(self):
        assert isclass(totest.EvolutionWarning)
        assert issubclass(totest.EvolutionWarning, totest.POSYDONWarning)

    def test_instance_InappropriateValueWarning(self):
        assert isclass(totest.InappropriateValueWarning)
        assert issubclass(totest.InappropriateValueWarning,\
                          totest.POSYDONWarning)

    def test_instance_IncompletenessWarning(self):
        assert isclass(totest.IncompletenessWarning)
        assert issubclass(totest.IncompletenessWarning, totest.POSYDONWarning)

    def test_instance_InterpolationWarning(self):
        assert isclass(totest.InterpolationWarning)
        assert issubclass(totest.InterpolationWarning, totest.POSYDONWarning)

    def test_instance_MissingFilesWarning(self):
        assert isclass(totest.MissingFilesWarning)
        assert issubclass(totest.MissingFilesWarning, totest.POSYDONWarning)

    def test_instance_OverwriteWarning(self):
        assert isclass(totest.OverwriteWarning)
        assert issubclass(totest.OverwriteWarning, totest.POSYDONWarning)

    def test_instance_ReplaceValueWarning(self):
        assert isclass(totest.ReplaceValueWarning)
        assert issubclass(totest.ReplaceValueWarning, totest.POSYDONWarning)

    def test_instance_UnsupportedModelWarning(self):
        assert isclass(totest.UnsupportedModelWarning)
        assert issubclass(totest.UnsupportedModelWarning,\
                          totest.POSYDONWarning)
        
    def test_instance_SFHModelWarning(self):
        assert isclass(totest.SFHModelWarning)
        assert issubclass(totest.SFHModelWarning, totest.POSYDONWarning)
        
    def test_instance_DeprecationWarning(self):
        assert isclass(totest.DeprecationWarning)
        assert issubclass(totest.DeprecationWarning, totest.POSYDONWarning)

    def test_instance_POSYDONWarning_subclasses(self):
        assert isinstance(totest._POSYDONWarning_subclasses, (dict))

    def test_instance_get_POSYDONWarning_class(self):
        assert isroutine(totest._get_POSYDONWarning_class)

    def test_instance_POSYDON_WARNINGS_REGISTRY(self):
        assert isinstance(totest._POSYDON_WARNINGS_REGISTRY, (dict))

    def test_instance_get_stats(self):
        assert isroutine(totest.get_stats)

    def test_instance_print_stats(self):
        assert isroutine(totest.print_stats)

    def test_instance_apply_POSYDON_filter(self):
        assert isroutine(totest._apply_POSYDON_filter)

    def test_instance_issue_warn(self):
        assert isroutine(totest._issue_warn)

    def test_instance_Caught_POSYDON_Warnings(self):
        assert isclass(totest._Caught_POSYDON_Warnings)

    def test_instance_CAUGHT_POSYDON_WARNINGS(self):
        assert isinstance(totest._CAUGHT_POSYDON_WARNINGS,\
                          totest._Caught_POSYDON_Warnings)

    def test_instance_Catch_POSYDON_Warnings(self):
        assert isclass(totest.Catch_POSYDON_Warnings)

    def test_instance_Pwarn(self):
        assert isroutine(totest.Pwarn)

    def test_instance_SetPOSYDONWarnings(self):
        assert isroutine(totest.SetPOSYDONWarnings)

    def test_instance_NoPOSYDONWarnings(self):
        assert isroutine(totest.NoPOSYDONWarnings)

    def test_instance_AllPOSYDONWarnings(self):
        assert isroutine(totest.AllPOSYDONWarnings)


class TestValues:
    # check that the values fit
    def test_value_POSYDONWarning_subclasses(self):
        assert 'ApproximationWarning' in totest._POSYDONWarning_subclasses

    def test_value_POSYDON_WARNINGS_REGISTRY(self):
        assert totest._POSYDON_WARNINGS_REGISTRY == {}

    def test_value_CAUGHT_POSYDON_WARNINGS(self):
        assert totest._CAUGHT_POSYDON_WARNINGS.catch_warnings == False
        assert totest._CAUGHT_POSYDON_WARNINGS.record
        assert totest._CAUGHT_POSYDON_WARNINGS.filter_first
        assert totest._CAUGHT_POSYDON_WARNINGS.caught_warnings == []
        assert totest._CAUGHT_POSYDON_WARNINGS.registry is\
               totest._POSYDON_WARNINGS_REGISTRY


class TestFunctions:
    @fixture
    def clear_registry(self):
        yield
        # empty the global POSYDON warnings registry after each test
        keys = []
        for k in totest._POSYDON_WARNINGS_REGISTRY:
            keys.append(k)
        for k in keys:
            del totest._POSYDON_WARNINGS_REGISTRY[k]

    @fixture
    def reset_filter(self):
        yield
        # set POSYDON warnings back to default
        totest.warnings.filterwarnings(action='ignore',\
                                       category=ResourceWarning)
        totest.warnings.filterwarnings(action='default',\
                                       category=totest.POSYDONWarning)

    # test functions
    def test_get_POSYDONWarning_class(self):
        # missing argument
        with raises(TypeError, match="missing 1 required positional argument"):
            totest._get_POSYDONWarning_class()
        # default to POSYDONWarning
        assert totest._get_POSYDONWarning_class("") == totest.POSYDONWarning
        assert totest._get_POSYDONWarning_class(totest.POSYDONWarning) ==\
               totest.POSYDONWarning
        # check subclasses of POSYDONWarning
        for (k, v) in totest._POSYDONWarning_subclasses.items():
            assert totest._get_POSYDONWarning_class(k) == v
            assert totest._get_POSYDONWarning_class(v) == v
        # bad input
        assert totest._get_POSYDONWarning_class(1) is None

    def test_get_stats(self):
        assert totest.get_stats() == totest._POSYDON_WARNINGS_REGISTRY

    def test_print_stats(self, capsys, clear_registry):
        # empty the global POSYDON warnings registry before this test
        keys = []
        for k in totest._POSYDON_WARNINGS_REGISTRY:
            keys.append(k)
        for k in keys:
            del totest._POSYDON_WARNINGS_REGISTRY[k]
        # no warnings to print
        totest.print_stats()
        assert "No POSYDON warnings occured.\n" == capsys.readouterr().out
        # add an artifical entry in the warnings registry and get the printout
        totest._POSYDON_WARNINGS_REGISTRY = {'Unit': 'Test'}
        totest.print_stats()
        assert "There have been POSYDON warnings in the global registry:\n "\
               + str(totest._POSYDON_WARNINGS_REGISTRY) + "\n" ==\
               capsys.readouterr().out

    def test_apply_POSYDON_filter(self, capsys, clear_registry, reset_filter):
        # wrong arguments
        with raises(TypeError, match="warning must be a dictionary"):
            totest._apply_POSYDON_filter(warning="Test")
        with raises(TypeError, match="message must be a string"):
            totest._apply_POSYDON_filter(warning={'message': 1})
        with raises(TypeError, match="stacklevel must be an integer"):
            totest._apply_POSYDON_filter(warning={'stacklevel': "Test"})
        with raises(TypeError, match="registry must be a dictionary or None"):
            totest._apply_POSYDON_filter(registry="Test")
        # check default warning
        assert totest._apply_POSYDON_filter() == dict(message="No warning")
        # check that further default warnings are filtered out but added to the
        # registry
        for i in range(10):
            for (k, v) in totest._POSYDON_WARNINGS_REGISTRY.items():
                assert v == i + 1
            assert totest._apply_POSYDON_filter() is None
        # check usage of python filter with a python warning
        totest.warnings.filterwarnings(action='ignore',\
                                       category=ResourceWarning)
        assert totest._apply_POSYDON_filter(warning={'message': "Test",\
                                                     'category':\
                                                     ResourceWarning}) is None
        # check the route of an always filter
        totest.warnings.filterwarnings(action='always',\
                                       category=ResourceWarning)
        for i in range(10):
            assert totest._apply_POSYDON_filter(warning={'message': "Test"\
                                                                    +str(i),\
                                                         'category':\
                                                         ResourceWarning}) ==\
                   {'message': "Test"+str(i), 'category': ResourceWarning}
            for (k, v) in totest._POSYDON_WARNINGS_REGISTRY.items():
                if "ResourceWarning" in k:
                    assert v == i
        # check the route of an error filter
        totest.warnings.filterwarnings(action='error',\
                                       category=ResourceWarning)
        assert totest._apply_POSYDON_filter(warning={'message': "Test"+str(i),\
                                                     'category':\
                                                     ResourceWarning}) ==\
               {'message': "Test"+str(i), 'category': ResourceWarning}

    def test_issue_warn(self, capsys, monkeypatch, recwarn, clear_registry):
        def mock_apply_POSYDON_filter(warning=dict(message="No warning"),\
                                      registry=None):
            return None

        # wrong arguments
        with raises(TypeError, match="warning must be a dictionary"):
            totest._issue_warn(warning="Test")
        with raises(TypeError, match="message must be a string"):
            totest._issue_warn(warning={'message': 1})
        with raises(TypeError, match="stacklevel must be an integer"):
            totest._issue_warn(warning={'stacklevel': "Test"})
        with raises(TypeError, match="registry must be a dictionary or None"):
            totest._issue_warn(registry="Test")
        # check default warning
        with warns(UserWarning, match="No warning") as warn1:
            totest._issue_warn()
        assert len(warn1) == 1
        # check filtered warning
        monkeypatch.setattr(totest, "_apply_POSYDON_filter", mock_apply_POSYDON_filter)
        totest._issue_warn()
        assert len(recwarn) == 0

    def test_Pwarn(self):
        # missing argument
        with raises(TypeError, match="missing 1 required positional argument"):
            totest.Pwarn()
        # wrong arguments
        with raises(TypeError, match="message must be a string"):
            totest.Pwarn(1)
        with raises(TypeError, match="stacklevel must be an integer"):
            totest.Pwarn("Unit", stacklevel="Test")
        # check output of POSYDONwarning
        with warns(totest.POSYDONWarning, match="Unit test") as warn1:
            totest.Pwarn("Unit test", "POSYDONWarning")
        assert len(warn1) == 1
        # check output of FutureWarning
        with warns(FutureWarning, match="Unit test") as warn2:
            totest.Pwarn("Unit test", FutureWarning)
        assert len(warn2) == 1
        # check output of UserWarning (default if unspecified)
        with warns(UserWarning, match="Unit test") as warn3:
            totest.Pwarn("Unit test")
        assert len(warn3) == 1

    def test_SetPOSYDONWarnings(self):
        totest.SetPOSYDONWarnings(action="once")
        assert str(totest.warnings.filters[0]) == "('once', None, <class "\
               + "'posydon.utils.posydonwarning.POSYDONWarning'>, None, 0)"
        totest.SetPOSYDONWarnings()
        assert str(totest.warnings.filters[0]) == "('default', None, <class "\
               + "'posydon.utils.posydonwarning.POSYDONWarning'>, None, 0)"
        # non POSYDON warnings have no effect
        totest.SetPOSYDONWarnings(category=UserWarning)
        assert str(totest.warnings.filters[0]) == "('default', None, <class "\
               + "'posydon.utils.posydonwarning.POSYDONWarning'>, None, 0)"

    def test_NoPOSYDONWarnings(self, reset_filter):
        totest.NoPOSYDONWarnings()
        assert str(totest.warnings.filters[0]) == "('ignore', None, <class "\
               + "'posydon.utils.posydonwarning.POSYDONWarning'>, None, 0)"

    def test_AllPOSYDONWarnings(self, reset_filter):
        totest.AllPOSYDONWarnings()
        assert str(totest.warnings.filters[0]) == "('always', None, <class "\
               + "'posydon.utils.posydonwarning.POSYDONWarning'>, None, 0)"


class TestPOSYDONWarning:
    @fixture
    def POSYDONWarning(self):
        # initialize an instance of the class with defaults
        return totest.POSYDONWarning()

    # test the POSYDONWarning class
    def test_init(self, POSYDONWarning):
        assert isroutine(POSYDONWarning.__init__)
        # check that the instance is of correct type and all code in the
        # __init__ got executed: the elements are created and initialized
        assert isinstance(POSYDONWarning, totest.POSYDONWarning)
        assert POSYDONWarning.message == ''

    def test_str(self, POSYDONWarning):
        assert isroutine(POSYDONWarning.__str__)
        assert str(POSYDONWarning) == "''"


class TestApproximationWarning:
    @fixture
    def ApproximationWarning(self):
        # initialize an instance of the class with defaults
        return totest.ApproximationWarning()

    # test the ApproximationWarning class
    def test_init(self, ApproximationWarning):
        assert isroutine(ApproximationWarning.__init__)
        # check that the instance is of correct type and all code in the
        # __init__ got executed: the elements are created and initialized
        assert isinstance(ApproximationWarning, totest.ApproximationWarning)
        assert ApproximationWarning.message == ''


class TestBinaryParsingWarning:
    @fixture
    def BinaryParsingWarning(self):
        # initialize an instance of the class with defaults
        return totest.BinaryParsingWarning()

    # test the BinaryParsingWarning class
    def test_init(self, BinaryParsingWarning):
        assert isroutine(BinaryParsingWarning.__init__)
        # check that the instance is of correct type and all code in the
        # __init__ got executed: the elements are created and initialized
        assert isinstance(BinaryParsingWarning, totest.BinaryParsingWarning)
        assert BinaryParsingWarning.message == ''


class TestClassificationWarning:
    @fixture
    def ClassificationWarning(self):
        # initialize an instance of the class with defaults
        return totest.ClassificationWarning()

    # test the ClassificationWarning class
    def test_init(self, ClassificationWarning):
        assert isroutine(ClassificationWarning.__init__)
        # check that the instance is of correct type and all code in the
        # __init__ got executed: the elements are created and initialized
        assert isinstance(ClassificationWarning, totest.ClassificationWarning)
        assert ClassificationWarning.message == ''


class TestEvolutionWarning:
    @fixture
    def EvolutionWarning(self):
        # initialize an instance of the class with defaults
        return totest.EvolutionWarning()

    # test the EvolutionWarning class
    def test_init(self, EvolutionWarning):
        assert isroutine(EvolutionWarning.__init__)
        # check that the instance is of correct type and all code in the
        # __init__ got executed: the elements are created and initialized
        assert isinstance(EvolutionWarning, totest.EvolutionWarning)
        assert EvolutionWarning.message == ''


class TestInappropriateValueWarning:
    @fixture
    def InappropriateValueWarning(self):
        # initialize an instance of the class with defaults
        return totest.InappropriateValueWarning()

    # test the InappropriateValueWarning class
    def test_init(self, InappropriateValueWarning):
        assert isroutine(InappropriateValueWarning.__init__)
        # check that the instance is of correct type and all code in the
        # __init__ got executed: the elements are created and initialized
        assert isinstance(InappropriateValueWarning, totest.InappropriateValueWarning)
        assert InappropriateValueWarning.message == ''


class TestIncompletenessWarning:
    @fixture
    def IncompletenessWarning(self):
        # initialize an instance of the class with defaults
        return totest.IncompletenessWarning()

    # test the IncompletenessWarning class
    def test_init(self, IncompletenessWarning):
        assert isroutine(IncompletenessWarning.__init__)
        # check that the instance is of correct type and all code in the
        # __init__ got executed: the elements are created and initialized
        assert isinstance(IncompletenessWarning, totest.IncompletenessWarning)
        assert IncompletenessWarning.message == ''


class TestInterpolationWarning:
    @fixture
    def InterpolationWarning(self):
        # initialize an instance of the class with defaults
        return totest.InterpolationWarning()

    # test the InterpolationWarning class
    def test_init(self, InterpolationWarning):
        assert isroutine(InterpolationWarning.__init__)
        # check that the instance is of correct type and all code in the
        # __init__ got executed: the elements are created and initialized
        assert isinstance(InterpolationWarning, totest.InterpolationWarning)
        assert InterpolationWarning.message == ''


class TestMissingFilesWarning:
    @fixture
    def MissingFilesWarning(self):
        # initialize an instance of the class with defaults
        return totest.MissingFilesWarning()

    # test the MissingFilesWarning class
    def test_init(self, MissingFilesWarning):
        assert isroutine(MissingFilesWarning.__init__)
        # check that the instance is of correct type and all code in the
        # __init__ got executed: the elements are created and initialized
        assert isinstance(MissingFilesWarning, totest.MissingFilesWarning)
        assert MissingFilesWarning.message == ''


class TestOverwriteWarning:
    @fixture
    def OverwriteWarning(self):
        # initialize an instance of the class with defaults
        return totest.OverwriteWarning()

    # test the OverwriteWarning class
    def test_init(self, OverwriteWarning):
        assert isroutine(OverwriteWarning.__init__)
        # check that the instance is of correct type and all code in the
        # __init__ got executed: the elements are created and initialized
        assert isinstance(OverwriteWarning, totest.OverwriteWarning)
        assert OverwriteWarning.message == ''


class TestReplaceValueWarning:
    @fixture
    def ReplaceValueWarning(self):
        # initialize an instance of the class with defaults
        return totest.ReplaceValueWarning()

    # test the ReplaceValueWarning class
    def test_init(self, ReplaceValueWarning):
        assert isroutine(ReplaceValueWarning.__init__)
        # check that the instance is of correct type and all code in the
        # __init__ got executed: the elements are created and initialized
        assert isinstance(ReplaceValueWarning, totest.ReplaceValueWarning)
        assert ReplaceValueWarning.message == ''


class TestUnsupportedModelWarning:
    @fixture
    def UnsupportedModelWarning(self):
        # initialize an instance of the class with defaults
        return totest.UnsupportedModelWarning()

    # test the UnsupportedModelWarning class
    def test_init(self, UnsupportedModelWarning):
        assert isroutine(UnsupportedModelWarning.__init__)
        # check that the instance is of correct type and all code in the
        # __init__ got executed: the elements are created and initialized
        assert isinstance(UnsupportedModelWarning,\
                          totest.UnsupportedModelWarning)
        assert UnsupportedModelWarning.message == ''


class TestSFHModelWarning:
    @fixture
    def SFHModelWarning(self):
        # initialize an instance of the class with defaults
        return totest.SFHModelWarning()

    # test the SFHModelWarning class
    def test_init(self, SFHModelWarning):
        assert isroutine(SFHModelWarning.__init__)
        # check that the instance is of correct type and all code in the
        # __init__ got executed: the elements are created and initialized
        assert isinstance(SFHModelWarning, totest.SFHModelWarning)
        assert SFHModelWarning.message == ''


class TestDeprecationWarning:
    @fixture
    def DeprecationWarning(self):
        # initialize an instance of the class with defaults
        return totest.DeprecationWarning()

    # test the DeprecationWarning class
    def test_init(self, DeprecationWarning):
        assert isroutine(DeprecationWarning.__init__)
        # check that the instance is of correct type and all code in the
        # __init__ got executed: the elements are created and initialized
        assert isinstance(DeprecationWarning, totest.DeprecationWarning)
        assert DeprecationWarning.message == ''


class Test_Caught_POSYDON_Warnings:
    @fixture
    def clear_registry(self):
        yield
        # empty the global POSYDON warnings registry after each test
        keys = []
        for k in totest._POSYDON_WARNINGS_REGISTRY:
            keys.append(k)
        for k in keys:
            del totest._POSYDON_WARNINGS_REGISTRY[k]

    @fixture
    def _Caught_POSYDON_Warnings(self, clear_registry):
        # initialize an instance of the class with defaults
        _Caught_POSYDON_Warnings = totest._Caught_POSYDON_Warnings()
        yield _Caught_POSYDON_Warnings
        # empty the cache to not print remaining records
        _Caught_POSYDON_Warnings.caught_warnings = []

    @fixture
    def test_dict(self):
        # a dictionary as a registry for testing
        return {'Unit': "Test"}

    # test the _Caught_POSYDON_Warnings class
    def test_init(self, _Caught_POSYDON_Warnings):
        assert isroutine(_Caught_POSYDON_Warnings.__init__)
        # check that the instance is of correct type and all code in the
        # __init__ got executed: the elements are created and initialized
        assert isinstance(_Caught_POSYDON_Warnings,\
                          totest._Caught_POSYDON_Warnings)
        assert _Caught_POSYDON_Warnings.catch_warnings == False
        assert _Caught_POSYDON_Warnings.caught_warnings == []
        assert _Caught_POSYDON_Warnings.record
        assert _Caught_POSYDON_Warnings.filter_first
        assert _Caught_POSYDON_Warnings._got_called == False
        assert _Caught_POSYDON_Warnings.registry is\
               totest._POSYDON_WARNINGS_REGISTRY
        # bad input
        with raises(TypeError, match="catch_warnings must be a boolean"):
            totest._Caught_POSYDON_Warnings(catch_warnings="Test")
        with raises(TypeError, match="record must be a boolean"):
            totest._Caught_POSYDON_Warnings(record="Test")
        with raises(TypeError, match="filter_first must be a boolean"):
            totest._Caught_POSYDON_Warnings(filter_first="Test")
        with raises(TypeError, match="registry must be a dictionary"):
            totest._Caught_POSYDON_Warnings(registry="Test")

    def test_str(self, _Caught_POSYDON_Warnings, test_dict):
        assert isroutine(_Caught_POSYDON_Warnings.__str__)
        assert str(_Caught_POSYDON_Warnings) == "POSYDON warnings are shown."
        # check with different setups
        test_cases = [{'catch_warnings': True, 'record': True,\
                       'filter_first': True, 'registry': None,\
                       'str': "POSYDON warnings will be caught and recorded. "\
                              +"Filters are applied before recording."},\
                      {'catch_warnings': True, 'record': True,\
                       'filter_first': False, 'registry': None,\
                       'str': "POSYDON warnings will be caught and recorded."\
                      },\
                      {'catch_warnings': True, 'record': False,\
                       'filter_first': True, 'registry': None,\
                       'str': "POSYDON warnings will be caught and discarded."\
                      },\
                      {'catch_warnings': True, 'record': False,\
                       'filter_first': False, 'registry': None,\
                       'str': "POSYDON warnings will be caught and discarded."\
                      },
                      {'catch_warnings': False, 'record': True,\
                       'filter_first': True, 'registry': test_dict,\
                       'str': "Currently a private registry is used, it "\
                              +"contains:\n{'Unit': 'Test'}"}]
        for tc in test_cases:
            caught = totest._Caught_POSYDON_Warnings(catch_warnings=\
                                                     tc['catch_warnings'],\
                                                     record=tc['record'],\
                                                     filter_first=\
                                                     tc['filter_first'],\
                                                     registry=tc['registry'])
            assert tc['str'] in str(caught)
        # check artifical caught_warnings
        for i in range(4):
            _Caught_POSYDON_Warnings.caught_warnings = i * [test_dict]
            if i==1:
                assert "There is 1 warning recorded." in\
                       str(_Caught_POSYDON_Warnings)
            elif i>1:
                assert "There are "+str(i)+" warnings recorded.".format(i) in\
                       str(_Caught_POSYDON_Warnings)
            else:
                assert "recorded" not in str(_Caught_POSYDON_Warnings)

    def test_call(self, _Caught_POSYDON_Warnings, test_dict):
        assert isroutine(_Caught_POSYDON_Warnings.__call__)
        assert _Caught_POSYDON_Warnings._got_called == False
        # bad input
        with raises(ValueError, match="Nothing to do: either empty_cache has "\
                                      +"to be True or new_warning/"\
                                      +"change_settings needs to be set."):
            _Caught_POSYDON_Warnings()
        assert _Caught_POSYDON_Warnings._got_called
        _Caught_POSYDON_Warnings.caught_warnings = [test_dict]
        _Caught_POSYDON_Warnings(empty_cache=True)
        assert _Caught_POSYDON_Warnings.caught_warnings == []
        # change_settings
        with raises(TypeError, match="change_settings has to be a dict"):
            _Caught_POSYDON_Warnings(change_settings="Test")
        for s in _Caught_POSYDON_Warnings.__dict__:
            if s == 'caught_warnings':
                _Caught_POSYDON_Warnings(change_settings={s: "Test"})
            elif s == 'registry':
                _Caught_POSYDON_Warnings(change_settings={s: "Test"})
                _Caught_POSYDON_Warnings(change_settings={s: None})
            else:
                with raises(TypeError, match="has to be a"):
                    _Caught_POSYDON_Warnings(change_settings={s: "Test"})
        with raises(AttributeError, match="unknown to "\
                                          +"_Caught_POSYDON_Warnings"):
            _Caught_POSYDON_Warnings(change_settings=test_dict)
        # change setting to catch warnings and add a new one to the record list
        _Caught_POSYDON_Warnings(change_settings={'catch_warnings': True},\
                                 new_warning={'message': "Test"})
        assert _Caught_POSYDON_Warnings.catch_warnings
        assert _Caught_POSYDON_Warnings.caught_warnings ==\
               [{'message': "Test"}]
        # unset catching: recorded warnings are issued and list is emptied
        with warns(UserWarning, match="Test"):
            _Caught_POSYDON_Warnings(change_settings={'catch_warnings': False})
        assert _Caught_POSYDON_Warnings.caught_warnings == []
        # change setting to catch warnings and add a new one to the record list
        _Caught_POSYDON_Warnings(change_settings={'catch_warnings': True,\
                                                  'filter_first': False},\
                                 new_warning={'message': "Test",\
                                              'stacklevel': -1})
        assert _Caught_POSYDON_Warnings.catch_warnings
        assert _Caught_POSYDON_Warnings.filter_first == False
        assert _Caught_POSYDON_Warnings.caught_warnings ==\
               [{'message': "Test", 'stacklevel': -1}]
        # change setting to record warnings and try to add a new one
        _Caught_POSYDON_Warnings(change_settings={'record': False},\
                                 new_warning={'message': "no record"})
        assert _Caught_POSYDON_Warnings.record == False
        assert _Caught_POSYDON_Warnings.caught_warnings ==\
               [{'message': "Test", 'stacklevel': -1}] # still has old record
        # unset catching: recorded warnings are issued and list is emptied
        with warns(UserWarning, match="Test") as winfo:
            _Caught_POSYDON_Warnings(change_settings={'catch_warnings': False})
        assert "posydonwarning.py" in winfo._list[0].filename
        assert _Caught_POSYDON_Warnings.caught_warnings == []

    def test_del(self, capsys, _Caught_POSYDON_Warnings):
        assert isroutine(_Caught_POSYDON_Warnings.__del__)
        # create an object with a catched warning, call the destructor and
        # empty it while recording the stdout and stderr
        cpw = totest._Caught_POSYDON_Warnings(catch_warnings=True)
        cpw(new_warning={'message': "Unit Test"})
        with warns(UserWarning, match="Unit Test"):
            cpw.__del__()
        assert capsys.readouterr().err ==\
               "There are still recorded warnings:\n"
        assert cpw.catch_warnings == False
        cpw.caught_warnings[0]['message'] += "Test"
        cpw.filter_first = False
        with warns(UserWarning, match="Unit TestTest"):
            cpw.__del__()
        assert capsys.readouterr().err ==\
               "There are still recorded warnings:\n"
        cpw.caught_warnings = []

    def test_got_called(self, _Caught_POSYDON_Warnings):
        assert isroutine(_Caught_POSYDON_Warnings.got_called)
        assert _Caught_POSYDON_Warnings.got_called() == False
        _Caught_POSYDON_Warnings._got_called = True
        assert _Caught_POSYDON_Warnings.got_called()

    def test_has_records(self, _Caught_POSYDON_Warnings, test_dict):
        assert isroutine(_Caught_POSYDON_Warnings.has_records)
        assert _Caught_POSYDON_Warnings.has_records() == False
        _Caught_POSYDON_Warnings.caught_warnings = [test_dict]
        assert _Caught_POSYDON_Warnings.has_records()

    def test_get_cache(self, _Caught_POSYDON_Warnings, test_dict):
        assert isroutine(_Caught_POSYDON_Warnings.get_cache)
        assert _Caught_POSYDON_Warnings.get_cache() == []
        _Caught_POSYDON_Warnings.caught_warnings = [test_dict]
        assert _Caught_POSYDON_Warnings.get_cache() == [test_dict]
        # clear the cache
        assert _Caught_POSYDON_Warnings.get_cache(empty_cache=True) ==\
               [test_dict]
        assert _Caught_POSYDON_Warnings.get_cache() == []

    def test_reset_cache(self, _Caught_POSYDON_Warnings, test_dict):
        assert isroutine(_Caught_POSYDON_Warnings.reset_cache)
        _Caught_POSYDON_Warnings.caught_warnings = [test_dict]
        _Caught_POSYDON_Warnings.reset_cache()
        assert _Caught_POSYDON_Warnings.caught_warnings == []


class TestCatch_POSYDON_Warnings:
    @fixture
    def clear_registry(self):
        yield
        # empty the global POSYDON warnings registry after each test
        keys = []
        for k in totest._POSYDON_WARNINGS_REGISTRY:
            keys.append(k)
        for k in keys:
            del totest._POSYDON_WARNINGS_REGISTRY[k]

    @fixture
    def Catch_POSYDON_Warnings(self, clear_registry):
        # initialize an instance of the class with defaults
        return totest.Catch_POSYDON_Warnings()

    # test the Catch_POSYDON_Warnings class
    def test_init(self, Catch_POSYDON_Warnings):
        assert isroutine(Catch_POSYDON_Warnings.__init__)
        # check that the instance is of correct type and all code in the
        # __init__ got executed: the elements are created and initialized
        assert Catch_POSYDON_Warnings.catch_warnings
        assert Catch_POSYDON_Warnings.record
        assert Catch_POSYDON_Warnings.filter_first
        assert Catch_POSYDON_Warnings.context_registry is None
        assert Catch_POSYDON_Warnings.python_catch is None
        # check other inputs
        cpw = totest.Catch_POSYDON_Warnings(own_registry=True,\
                                            use_python_catch=True)
        assert cpw.context_registry == {}
        assert isinstance(cpw.python_catch, totest.warnings.catch_warnings)

    def test_enter_exit(self, Catch_POSYDON_Warnings):
        assert isroutine(Catch_POSYDON_Warnings.__enter__)

    def test_exit(self, Catch_POSYDON_Warnings):
        assert isroutine(Catch_POSYDON_Warnings.__exit__)

    def test_context(self, capsys):
        with totest.Catch_POSYDON_Warnings() as cpw:
            assert cpw is totest._CAUGHT_POSYDON_WARNINGS
            assert cpw.catch_warnings
            assert cpw.record
            assert cpw.filter_first
            assert cpw.registry is totest._POSYDON_WARNINGS_REGISTRY
        with totest.Catch_POSYDON_Warnings(use_python_catch=True) as cpw:
            assert cpw is totest._CAUGHT_POSYDON_WARNINGS
            assert cpw.catch_warnings
            assert cpw.record
            assert cpw.filter_first
            assert cpw.registry is totest._POSYDON_WARNINGS_REGISTRY
