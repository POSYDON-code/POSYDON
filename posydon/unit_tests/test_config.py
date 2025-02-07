"""Unit tests of posydon/config.py
"""

__authors__ = [
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>"
]

# import the module which will be tested
import posydon.config as totest
# aliases
os = totest.os

# import other needed code for the tests, which is not already imported in the
# module you like to test
from pytest import raises
from inspect import isroutine


# define test classes collecting several test functions
class TestElements:
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = ['PATH_TO_POSYDON', 'PATH_TO_POSYDON_DATA', '__builtins__',\
                    '__cached__', '__doc__', '__file__', '__loader__',\
                    '__name__', '__package__', '__spec__', 'ensure_path',\
                    'load_dotenv', 'os']
        assert dir(totest) == elements, "There might be added or removed "\
                                        + "objects without an update on the "\
                                        + "unit test."

    def test_instance_ensure_path(self):
        assert isroutine(totest.ensure_path)

    def test_instance_PATH_TO_POSYDON_DATA(self):
        assert isinstance(totest.PATH_TO_POSYDON_DATA, (str, bytes,\
                                                        os.PathLike))

    def test_instance_PATH_TO_POSYDON(self):
        assert isinstance(totest.PATH_TO_POSYDON, (str, bytes, os.PathLike))


class TestValues:
    # check that the values fit
    def test_value_PATH_TO_POSYDON(self):
        assert totest.PATH_TO_POSYDON is not None

    def test_value_PATH_TO_POSYDON_DATA(self):
        assert totest.PATH_TO_POSYDON_DATA is not None


class TestFunctions:
    # test functions
    def test_ensure_path(self):
        # missing argument
        with raises(TypeError, match="missing 1 required positional "\
                                     +"argument: 'name'"):
            totest.ensure_path()
        # bad input
        with raises(TypeError, match="'name' has to be a string."):
            totest.ensure_path(None)
        # bad input
        with raises(NameError, match="DOES_NOT_EXIST is not defined in the "\
                                     +"environment."):
            totest.ensure_path("DOES_NOT_EXIST")
        # bad input
        for k,i in os.environ.items():
            # search for first item, which is no path to check the error
            if not os.path.isdir(i):
                with raises(NotADirectoryError, match=f"{i} given in {k} is "\
                                                      +"an invalid path."):
                    totest.ensure_path(k)
                break
        # example
        for n in ["PATH_TO_POSYDON", "PATH_TO_POSYDON_DATA"]:
            assert totest.ensure_path(n) == os.getenv(n)
