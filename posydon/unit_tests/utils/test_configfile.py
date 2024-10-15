"""Unit tests of posydon/utils/configfile.py
"""

__authors__ = [
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>"
]

# import the module which will be tested
import posydon.utils.configfile as totest

# import other needed code for the tests, which is not already imported in the
# module you like to test
from pytest import fixture, raises
from inspect import isclass, isroutine
from ast import AST, parse

# define test classes collecting several test functions
class TestElements:
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = ['ConfigFile', 'VariableKey', '__authors__', '__builtins__',
                    '__cached__', '__doc__', '__file__', '__loader__',
                    '__name__', '__package__', '__spec__', 'ast',
                    'configparser', 'copy', 'json', 'np', 'operator', 'os',
                    'parse_inifile']
        assert dir(totest) == elements, "There might be added or removed "+\
               "objects without an update on the unit test."

    def test_instance_ConfigFile(self):
        assert isclass(totest.ConfigFile)

    def test_instance_parse_inifile(self):
        assert isroutine(totest.parse_inifile)

    def test_instance_VariableKey(self):
        assert isclass(totest.VariableKey)


class TestFunctions:
    @fixture
    def test_path(self, tmp_path):
        # a temporary path for testing
        return totest.os.path.join(tmp_path, "test.ini")

    @fixture
    def test_ini(self, test_path):
        # a temporary ini file for testing"
        with open(test_path, "w") as test_file:
            test_file.write("[run_parameters]\n")
            test_file.write("test_bool1 = True\n")
            test_file.write("test_bool2 = False\n")
            test_file.write("test_None = None\n")
            test_file.write("\n[mesa_inlists]\n")
            test_file.write("test_float = 0.1\n")
            test_file.write("test_int = 10\n")
            test_file.write("test_list = ['Unit Test1', 'Unit Test2']\n")
            test_file.write("test_str1 = 'Unit Test'\n")
            test_file.write("test_str2 = 'Unit,Test'\n")
            test_file.write("test_BinOp = ${test_int}+1\n")
            test_file.write("\n[mesa_extras]\n")
            test_file.write("test_else = print('1')\n")
            test_file.write("\n[slurm]\n")
            test_file.write("test_exception = ast.parse('X=1', mode='eval')\n")
        return

    # test functions
    def test_parse_inifile(self, test_path, test_ini, monkeypatch):
        def mock_parse(source, filename='<unknown>', mode='exec', *,\
                       type_comments=False, feature_version=None, optimize=-1):
            mock_node = AST()
            mock_node.body = parse(source, mode='eval')
            return mock_node
        # missing argument
        with raises(TypeError, match="missing 1 required positional"+\
                    " argument: 'inifile'"):
            totest.parse_inifile()
        # bad input
        with raises(TypeError, match="'int' object is not iterable"):
            totest.parse_inifile(1)
        # read test ini
        rp, s, mi, me = totest.parse_inifile(test_path)
        assert rp == {'MESA_DIR'.lower(): totest.os.environ['MESA_DIR'],\
                      'test_bool1'.lower(): True, 'test_bool2'.lower(): False,\
                      'test_None'.lower(): None}
        assert s == {'MESA_DIR'.lower(): totest.os.environ['MESA_DIR'],\
                     'test_exception'.lower():\
                     ["ast.parse('X=1'", "mode='eval')"]}
        assert mi == {'MESA_DIR'.lower(): totest.os.environ['MESA_DIR'],\
                      'test_float'.lower(): 0.1, 'test_int'.lower(): 10,\
                      'test_list'.lower(): ['Unit Test1', 'Unit Test2'],\
                      'test_str1'.lower(): 'Unit Test',\
                      'test_str2'.lower(): ['Unit', 'Test'],\
                      'test_BinOp'.lower(): 11}
        assert me == {'MESA_DIR'.lower(): totest.os.environ['MESA_DIR'],\
                      'test_else'.lower(): "print('1')"}
        with monkeypatch.context() as mp:
            # replace parse to check iteration
            mp.setattr(totest.ast, "parse", mock_parse)
            # replace content of test.ini to check Error
            with open(test_path, "w") as test_file:
                test_file.write("[Unit Test]\n")
                test_file.write("test_value = 'Test'\n")
            with raises(KeyError, match="'run_parameters'"):
                totest.parse_inifile(test_path)


class TestConfigFile:
    @fixture
    def ConfigFile(self):
        # initialize an instance of the class with defaults
        ConfigFile = totest.ConfigFile()
        return ConfigFile

    @fixture
    def test_path(self, tmp_path):
        # a temporary path for testing
        return totest.os.path.join(tmp_path, "test.json")

    @fixture
    def test_json(self, test_path):
        # a temporary ini file for testing"
        with open(test_path, "w") as test_file:
            test_file.write('{\n    "Unit": "Test"\n}\n')
        return

    @fixture
    def test_entries(self):
        # entries for testing
        return {'Unit1': "Test", 'Unit2': "Again"}

    @fixture
    def ConfigFile_filled(self, test_entries):
        # initialize an instance of the class with defaults
        ConfigFile = totest.ConfigFile()
        ConfigFile.update(test_entries)
        return ConfigFile

    # test the ConfigFile class
    def test_init(self, ConfigFile, tmp_path):
        assert isroutine(ConfigFile.__init__)
        # check that the instance is of correct type and all code in the
        # __init__ got executed: the elements are created and initialized
        assert isinstance(ConfigFile, totest.ConfigFile)
        assert ConfigFile.entries == {}
        assert ConfigFile.path is None
        # error on directory
        with raises(IsADirectoryError, match="Is a directory: '"+\
                                             str(tmp_path)+"'"):
            totest.ConfigFile(path=tmp_path)
        # non-existing path
        test_path = totest.os.path.join(tmp_path, "does_not_exist.test")
        test_ConfigFile = totest.ConfigFile(test_path)
        assert test_ConfigFile.path == test_path

    def test_deepcopy(self, ConfigFile):
        assert isroutine(ConfigFile.deepcopy)
        ConfigFile.entries['Unit'] = "Test"
        ConfigFile.path = "Test_path"
        test_ConfigFile = ConfigFile.deepcopy()
        # not same object
        assert test_ConfigFile is not ConfigFile
        # but same content
        assert test_ConfigFile.entries == ConfigFile.entries
        assert test_ConfigFile.path == ConfigFile.path

    def test_serialize(self, ConfigFile):
        assert isroutine(ConfigFile._serialize)
        # serialize numpy array
        np_array = totest.np.array([2, 1, 0])
        assert ConfigFile._serialize(data=np_array) == [2, 1, 0]
        # other data don't give a return
        assert ConfigFile._serialize(data="Test") is None

    def test_save(self, ConfigFile, tmp_path):
        assert isroutine(ConfigFile.save)
        # bad input
        with raises(ValueError, match="No path passed."):
            ConfigFile.save()
        with raises(PermissionError, match="JSON file not saved: overwrite "+\
                                           "not permitted."):
            ConfigFile.save(path=tmp_path, overwrite=False)
        test_path = totest.os.path.join(tmp_path, "save.test")
        ConfigFile.path = test_path
        ConfigFile.entries['Unit'] = "Test"
        ConfigFile.save()
        assert totest.os.path.isfile(test_path)

    def test_load(self, ConfigFile, test_path, test_json):
        assert isroutine(ConfigFile.load)
        # bad input
        with raises(ValueError, match="No path passed."):
            ConfigFile.load()
        # load mock file with path as keyword
        ConfigFile.load(path=test_path)
        assert ConfigFile.entries == {'Unit': "Test"}
        # set path and entries
        ConfigFile.path = test_path
        ConfigFile.entries['Unit'] = "ToReplace"
        # try to load mock file without permission to overwrite
        with raises(PermissionError,\
                    match="Not allowed to update the entries"):
            ConfigFile.load()
        assert ConfigFile.entries == {'Unit': "ToReplace"}
        # load mock file and overwrite entry
        ConfigFile.load(can_update=True)
        assert ConfigFile.entries == {'Unit': "Test"}

    def test_getattr(self, ConfigFile):
        assert isroutine(ConfigFile.__getattr__)
        # bad input
        with raises(KeyError, match='Test'):
            ConfigFile.Test
        # add a value to entries and get it
        ConfigFile.entries['Unit'] = "Test"
        assert ConfigFile.Unit == "Test"

    def test_getitem(self, ConfigFile):
        assert isroutine(ConfigFile.__getitem__)
        # bad input
        with raises(KeyError, match='Test'):
            ConfigFile['Test']
        # add a value to entries and get it
        ConfigFile.entries['Unit'] = "Test"
        assert ConfigFile['Unit'] == "Test"

    def test_setitem(self, ConfigFile):
        assert isroutine(ConfigFile.__setitem__)
        # add a value to entries and get it
        ConfigFile['Unit'] = "Test"
        assert ConfigFile.entries['Unit'] == "Test"

    def test_delitem(self, ConfigFile_filled):
        assert isroutine(ConfigFile_filled.__delitem__)
        with raises(KeyError, match="'Unit'"):
            del ConfigFile_filled['Unit']
        # delete entries
        del ConfigFile_filled['Unit1']
        del ConfigFile_filled['Unit2']
        assert ConfigFile_filled.entries == {}

    def test_iter(self, ConfigFile_filled, test_entries):
        assert isroutine(ConfigFile_filled.__iter__)
        # loop over entries
        iteration = 0
        for key in ConfigFile_filled:
            iteration += 1
            assert ConfigFile_filled.entries[key] == test_entries[key]
        # check that it was iterated over all keys in test_entries
        assert iteration == len(test_entries)

    def test_update(self, ConfigFile, test_entries):
        assert isroutine(ConfigFile.update)
        ConfigFile.update(test_entries)
        assert ConfigFile.entries == test_entries

    def test_keys(self, ConfigFile_filled, test_entries):
        assert isroutine(ConfigFile_filled.keys)
        assert ConfigFile_filled.keys() == test_entries.keys()

    def test_values(self, ConfigFile_filled, test_entries):
        assert isroutine(ConfigFile_filled.values)
        assert isinstance(ConfigFile_filled.values(),\
                          type(test_entries.values()))
        iteration = 0
        for v in ConfigFile_filled.values():
            iteration += 1
            assert v in test_entries.values()
        # check that it was iterated over all values in test_entries
        assert iteration == len(test_entries.values())

    def test_items(self, ConfigFile_filled, test_entries):
        assert isroutine(ConfigFile_filled.items)
        assert ConfigFile_filled.items() == test_entries.items()

    def test_repr(self, ConfigFile_filled):
        assert isroutine(ConfigFile_filled.__repr__)
        assert repr(ConfigFile_filled) == "Unit1: Test\nUnit2: Again\n"

    def test_contains(self, ConfigFile_filled, test_entries):
        assert isroutine(ConfigFile_filled.__contains__)
        for key in test_entries:
            assert key in ConfigFile_filled

    def test_len(self, ConfigFile_filled, test_entries):
        assert isroutine(ConfigFile_filled.__len__)
        assert len(ConfigFile_filled) == len(test_entries)


#class TestVariableKey:
#    # test the VariableKey class: looks to be not used as long as using python3
