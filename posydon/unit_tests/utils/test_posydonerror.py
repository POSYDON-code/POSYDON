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

@fixture
def BinaryStar():
    # initialize a BinaryStar instance, which is a required argument
    return totest.BinaryStar()

@fixture
def SingleStar():
    # initialize a SingleStar instance, which is a required argument
    return totest.SingleStar()

# define test classes collecting several test functions
class TestElements:
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = ['BinaryStar', 'ClassificationError', 'FlowError', 'GridError', 'MatchingError',\
                    'ModelError', 'NumericalError', 'POSYDONError',\
                    'SingleStar', '__authors__', '__builtins__', '__cached__',\
                    '__doc__', '__file__', '__loader__', '__name__',\
                    '__package__', '__spec__', 'copy',\
                    'initial_condition_message']
        assert dir(totest) == elements, "There might be added or removed "\
                                        + "objects without an update on the "\
                                        + "unit test."

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

    def test_instance_initial_condition_message(self):
        assert isroutine(totest.initial_condition_message)


class TestFunctions:
    # test functions
    def test_initial_condition_message(self, BinaryStar, artificial_object):
        with raises(TypeError, match="The binary must be a BinaryStar object"):
            message = totest.initial_condition_message(binary=\
                                                       artificial_object)
        assert "Failed Binary Initial Conditions" in\
               totest.initial_condition_message(binary=BinaryStar)
        with raises(TypeError, match="is not iterable"):
            message = totest.initial_condition_message(binary=BinaryStar,\
                                                       ini_params=1)
        with raises(TypeError, match="can only concatenate str"):
            message = totest.initial_condition_message(binary=BinaryStar,\
                                                       ini_params=[1,2])
        assert totest.initial_condition_message(binary=BinaryStar, ini_params=\
                                                ["a: 1\n", "b: 2\n"]) ==\
               "a: 1\nb: 2\n"


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

    @fixture
    def POSYDONError_object(self, artificial_object):
        # initialize an instance of the class with an artifical object in a
        # list via key
        return totest.POSYDONError(objects=[artificial_object])

    @fixture
    def POSYDONError_SingleStar(self, SingleStar):
        # initialize an instance of the class with a SingleStar object
        return totest.POSYDONError(objects=SingleStar)

    @fixture
    def POSYDONError_BinaryStar(self, BinaryStar):
        # initialize an instance of the class with a BinaryStar object
        return totest.POSYDONError(objects=BinaryStar)

    @fixture
    def POSYDONError_List(self, artificial_object, SingleStar, BinaryStar):
        # initialize an instance of the class with a list of objects
        return totest.POSYDONError(objects=[artificial_object, SingleStar,\
                                            BinaryStar])

    # test the POSYDONError class
    def test_init(self, POSYDONError, POSYDONError_position, POSYDONError_key,\
                  POSYDONError_object, artificial_object):
        assert isroutine(POSYDONError.__init__)
        # check that the instance is of correct type and all code in the
        # __init__ got executed: the elements are created and initialized
        assert isinstance(POSYDONError, totest.POSYDONError)
        # check defaults
        assert POSYDONError.message == ""
        assert POSYDONError.objects is None
        # test a passed message via positional argument
        assert POSYDONError_position.message == "test message on position"
        assert POSYDONError_position.objects is None
        # test a passed message via key
        error_object = totest.POSYDONError(message="test message with key")
        assert POSYDONError_key.message == "test message with key"
        assert POSYDONError_key.objects is None
        # test a passed object
        assert POSYDONError_object.message == ""
        assert POSYDONError_object.objects == [artificial_object]
        # test requests on input parameters
        with raises(TypeError, match="message must be a string"):
            error_object = totest.POSYDONError(message=artificial_object)
        with raises(TypeError, match="objects must be None, a list, a "\
                                     +"SingleStar object, or a BinaryStar "\
                                     +"object"):
            error_object = totest.POSYDONError(objects=artificial_object)


    def test_str(self, POSYDONError_position, POSYDONError_object,\
                 POSYDONError_SingleStar, POSYDONError_BinaryStar,\
                 POSYDONError_List):
        assert isroutine(POSYDONError_position.__str__)
        assert str(POSYDONError_position) == "\ntest message on position"
        # test passed objects
        assert str(POSYDONError_object) == "\n"
        assert "OBJECT #(<class 'posydon.binary_evol.singlestar.SingleStar'>)"\
               in str(POSYDONError_SingleStar)
        assert "OBJECT #(<class 'posydon.binary_evol.binarystar.BinaryStar'>)"\
               in str(POSYDONError_BinaryStar)
        assert "OBJECT #1" not in str(POSYDONError_List)
        assert "OBJECT #2 (<class 'posydon.binary_evol.singlestar.SingleStar"\
               + "'>)" in str(POSYDONError_List)
        assert "OBJECT #3 (<class 'posydon.binary_evol.binarystar.BinaryStar"\
               + "'>)" in str(POSYDONError_List)
