"""Unit tests of posydon/popsyn/io.py
"""

__authors__ = [
    "Elizabeth Teng <elizabethteng@u.northwestern.edu>"
]

# import the module which will be tested
import posydon.popsyn.io as totest
# aliases
np = totest.np
pd = totest.pd

# import other needed code for the tests, which is not already imported in the
# module you like to test
from pytest import fixture, raises, warns, approx
from inspect import isroutine, isclass
import textwrap
import os
import errno
import pprint
from posydon.binary_evol.simulationproperties import SimulationProperties

# define test classes collecting several test functions
class TestElements:
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = ['BINARYPROPERTIES_DTYPES', \#'OBJECT_FIXED_SUB_DTYPES', \
                    'STARPROPERTIES_DTYPES', 'EXTRA_BINARY_COLUMNS_DTYPES', \
                    'EXTRA_STAR_COLUMNS_DTYPES', 'SCALAR_NAMES_DTYPES', \
                    'clean_binary_history_df', 'clean_binary_oneline_df', \
                    'parse_inifile', 'simprop_kwargs_from_ini', \
                    'binarypop_kwargs_from_ini', 'create_run_script_text', \
                    'create_merge_script_text', \
                    '__builtins__', '__cached__', '__doc__', '__file__',\
                    '__loader__', '__name__', '__package__', '__spec__', \
                    'ConfigParser', 'ast', 'importlib', 'os', 'errno', \
                    'pprint', 'np', 'pd','SimulationProperties']
        totest_elements = set(dir(totest))
        missing_in_test = set(elements) - totest_elements
        assert len(missing_in_test) == 0, "There are missing objects in "\
                                          +f"{totest.__name__}: "\
                                          +f"{missing_in_test}. Please "\
                                          +"check, whether they have been "\
                                          +"removed on purpose and update "\
                                          +"this unit test."
        new_in_test = totest_elements - set(elements)
        assert len(new_in_test) == 0, "There are new objects in "\
                                      +f"{totest.__name__}: {new_in_test}. "\
                                      +"Please check, whether they have been "\
                                      +"added on purpose and update this "\
                                      +"unit test."
 
class TestFunctions:
    
    @fixture
    def simple_ini(tmp_path):
        file_path = os.path.join(tmp_path, "test.ini")
        with open(file_path, "w") as f:
            f.write("[section]\nkey=value\n")
        return file_path
    
    @fixture
    def multi_ini(tmp_path):
        file1 = os.path.join(tmp_path, "a.ini")
        file2 = os.path.join(tmp_path, "b.ini")
        with open(file1, "w") as f:
            f.write("[section]\nkey1=value1\n")
        with open(file2, "w") as f:
            f.write("[section]\nkey2=value2\n")
        return [file1, file2]

    @fixture
    def empty_ini(tmp_path):
        file_path = os.path.join(tmp_path, "empty.ini")
        with open(file_path, "w") as f:
            f.write("")
        return file_path

    @fixture
    def textfile(tmp_path):
        file_path = os.path.join(tmp_path, "textfile.txt")
        with open(file_path, "w") as f:
            f.write("test")
        return file_path

    
#     def test_clean_binary_history_df():
#         pass
        
#     def test_clean_binary_oneline_df():
#         pass
        
    def test_parse_inifile():
        # missing argument
        with raises(TypeError, match="missing 1 required positional argument: 'path'"):
            totest.parse_inifile()
        # bad input
        with raises(FileNotFoundError, match="[Errno 2] No such file or directory: '/none.ini'"):
            totest.parse_inifile('/nonexistent.ini')
        with raises(FileNotFoundError, match="[Errno 2] No such file or directory: '/none.ini'"):
            totest.parse_inifile([simple_ini,'/nonexistent.ini'])
        with raises(MissingSectionHeaderError, match="File contains no section headers"):
            totest.parse_inifile(textfile)
        with raises(ValueError, match="Path must be a string or list of strings."):
            totest.parse_inifile(0)
            
        # example: single inifile
        parser = parse_inifile(simple_ini)
        assert isinstance(parser, ConfigParser)
        assert parser.has_section("section")
        assert parser.get("section", "key") == "value"
        
        # example: multiple inifiles
        parser = parse_inifile(multi_ini)
        assert parser.has_option("section", "key1")
        assert parser.has_option("section", "key2")

        
#     def test_simprop_kwargs_from_ini():
#         pass
        
#     def test_binarypop_kwargs_from_ini():
#         pass
        
    def test_create_run_script_text():
        # missing argument
        with raises(TypeError, match="missing 1 required positional argument: 'ini_file'"):
            totest.create_run_script_text()
        # bad input
        with raises(NameError, match="name 'testfile' is not defined"):
            totest.create_run_script_text(testfile.ini)
        # example        
        out = textwrap.dedent("""\
            from posydon.popsyn.binarypopulation import BinaryPopulation
            from posydon.popsyn.io import binarypop_kwargs_from_ini
            from posydon.utils.common_functions import convert_metallicity_to_string
            import argparse
            if __name__ == "__main__":
                parser = argparse.ArgumentParser()
                parser.add_argument('metallicity', type=float)
                args = parser.parse_args()
                ini_kw = binarypop_kwargs_from_ini('testfile.ini')
                ini_kw['metallicity'] = args.metallicity
                str_met = convert_metallicity_to_string(args.metallicity)
                ini_kw['temp_directory'] = str_met+'_Zsun_' + ini_kw['temp_directory']
                synpop = BinaryPopulation(**ini_kw)
                synpop.evolve()""")
        assert totest.create_run_script_text('testfile.ini') == out
        
        
    def test_create_merge_script_text():
        # missing argument
        with raises(TypeError, match="missing 1 required positional argument: 'ini_file'"):
            totest.create_merge_script_text()
        # bad input
        with raises(NameError, match="name 'testfile' is not defined"):
            totest.create_merge_script_text(testfile.ini)
        # example        
        out = textwrap.dedent("""\
            from posydon.popsyn.binarypopulation import BinaryPopulation
            from posydon.popsyn.io import binarypop_kwargs_from_ini
            from posydon.utils.common_functions import convert_metallicity_to_string
            import argparse
            import os
            if __name__ == "__main__":
                parser = argparse.ArgumentParser()
                parser.add_argument("metallicity", type=float)
                args = parser.parse_args()
                ini_kw = binarypop_kwargs_from_ini("testfile.ini")
                ini_kw["metallicity"] = args.metallicity
                str_met = convert_metallicity_to_string(args.metallicity)
                ini_kw["temp_directory"] = str_met+"_Zsun_" + ini_kw["temp_directory"]
                synpop = BinaryPopulation(**ini_kw)
                path_to_batch = ini_kw["temp_directory"]
                tmp_files = [os.path.join(path_to_batch, f) for f in os.listdir(path_to_batch) if os.path.isfile(os.path.join(path_to_batch, f))]
                tmp_files = sorted(tmp_files, key=lambda x: int(x.split(".")[-1]))
                synpop.combine_saved_files(str_met+ "_Zsun_population.h5", tmp_files)
                print("done")
                if len(os.listdir(path_to_batch)) == 0:
                    os.rmdir(path_to_batch)""")
        assert totest.create_merge_script_text('testfile.ini') == out

        
        
        
        
        