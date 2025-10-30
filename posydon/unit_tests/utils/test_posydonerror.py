"""Unit tests of posydon/utils/posydonerror.py
"""

__authors__ = [
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>"
]

# import the module which will be tested
import posydon.utils.posydonerror as totest

# import other needed code for the tests, which is not already imported in the
# module you like to test
from pytest import fixture, raises
from inspect import isclass, isroutine

@fixture
def artificial_object():
    # create a dict as test object
    return {'Test': 'object'}

# define test classes collecting several test functions
class TestElements:
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = {'ClassificationError', 'FlowError', 'GridError',\
                    'MatchingError', 'ModelError', 'NumericalError',\
                    'POSYDONError', '__authors__', '__builtins__',\
                    '__cached__', '__doc__', '__file__', '__loader__',\
                    '__name__', '__package__', '__spec__'}
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

    def test_instance_POSYDONError(self):
        assert isclass(totest.POSYDONError)
        assert issubclass(totest.POSYDONError, Exception)
        with raises(totest.POSYDONError, match="Test"):
            raise totest.POSYDONError("Test")
    
    def test_instance_ClassificationError(self):
        assert isclass(totest.ClassificationError)
        assert issubclass(totest.ClassificationError, totest.POSYDONError)
        with raises(totest.ClassificationError, match="Test"):
            raise totest.ClassificationError("Test")

    def test_instance_FlowError(self):
        assert isclass(totest.FlowError)
        assert issubclass(totest.FlowError, totest.POSYDONError)
        with raises(totest.FlowError, match="Test"):
            raise totest.FlowError("Test")

    def test_instance_GridError(self):
        assert isclass(totest.GridError)
        assert issubclass(totest.GridError, totest.POSYDONError)
        with raises(totest.GridError, match="Test"):
            raise totest.GridError("Test")

    def test_instance_MatchingError(self):
        assert isclass(totest.MatchingError)
        assert issubclass(totest.MatchingError, totest.POSYDONError)
        with raises(totest.MatchingError, match="Test"):
            raise totest.MatchingError("Test")

    def test_instance_ModelError(self):
        assert isclass(totest.ModelError)
        assert issubclass(totest.ModelError, totest.POSYDONError)
        with raises(totest.ModelError, match="Test"):
            raise totest.ModelError("Test")

    def test_instance_NumericalError(self):
        assert isclass(totest.NumericalError)
        assert issubclass(totest.NumericalError, totest.POSYDONError)
        with raises(totest.NumericalError, match="Test"):
            raise totest.NumericalError("Test")


class TestPOSYDONError:
    @fixture
    def POSYDONError(self):
        # initialize an instance of the class with defaults
        return totest.POSYDONError()

    @fixture
    def POSYDONError_position(self):
        # initialize an instance of the class with a positional argument
        return totest.POSYDONError("test message on position")

    @fixture
    def POSYDONError_key(self):
        # initialize an instance of the class with a message via key
        return totest.POSYDONError(message="test message with key")

    # test the POSYDONError class
    def test_init(self, POSYDONError, POSYDONError_position, POSYDONError_key):
        assert isroutine(POSYDONError.__init__)
        # check that the instance is of correct type and all code in the
        # __init__ got executed: the elements are created and initialized
        assert isinstance(POSYDONError, totest.POSYDONError)
        # check defaults
        assert POSYDONError.message == ""
        # test a passed message via positional argument
        assert POSYDONError_position.message == "test message on position"
        # test a passed message via key
        assert POSYDONError_key.message == "test message with key"
        # test requests on input parameters
        with raises(TypeError, match="message must be a string"):
            error_object = totest.POSYDONError(message=artificial_object)

    def test_str(self, POSYDONError_position):
        assert isroutine(POSYDONError_position.__str__)
        assert str(POSYDONError_position) == "test message on position"
