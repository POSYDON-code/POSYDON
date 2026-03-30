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

import ast
import errno
import importlib
import os
import pprint
import textwrap
from configparser import ConfigParser, MissingSectionHeaderError

# import other needed code for the tests, which is not already imported in the
# module you like to test
from pytest import approx, fixture, raises, warns


# define test classes collecting several test functions
class TestElements:
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = ['BINARYPROPERTIES_DTYPES', 'OBJECT_FIXED_SUB_DTYPES',
                    'STARPROPERTIES_DTYPES', 'EXTRA_BINARY_COLUMNS_DTYPES',
                    'EXTRA_STAR_COLUMNS_DTYPES', 'SCALAR_NAMES_DTYPES',
                    'clean_binary_history_df', 'clean_binary_oneline_df',
                    'parse_inifile', 'simprop_kwargs_from_ini',
                    'binarypop_kwargs_from_ini',
                    '__builtins__', '__cached__', '__doc__', '__file__',
                    '__loader__', '__name__', '__package__', '__spec__',
                    'ConfigParser', 'ast', 'importlib', 'os', 'errno',
                    'pprint', 'np', 'pd',]
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
    def simple_ini(self,tmp_path):
        file_path = os.path.join(tmp_path, "test.ini")
        with open(file_path, "w") as f:
            f.write("[section]\nkey=value\n")
        return file_path

    @fixture
    def multi_ini(self,tmp_path):
        file1 = os.path.join(tmp_path, "a.ini")
        file2 = os.path.join(tmp_path, "b.ini")
        with open(file1, "w") as f:
            f.write("[section]\nkey1=value1\n")
        with open(file2, "w") as f:
            f.write("[section]\nkey2=value2\n")
        return [file1, file2]

    @fixture
    def textfile(self,tmp_path):
        file_path = os.path.join(tmp_path, "textfile.txt")
        with open(file_path, "w") as f:
            f.write("test")
        return file_path

    @fixture
    def sim_ini(self,tmp_path):
        ini_content = """
        [flow]
        import = ['posydon.binary_evol.flow_chart', 'flow_chart']
        absolute_import = None

        [step_HMS_HMS]
        import = ['posydon.binary_evol.MESA.step_mesa', 'MS_MS_step']
        absolute_import = None
        interpolation_method = 'linear3c_kNN'
        save_initial_conditions = True
        verbose = False

        [extra_hooks]
        import_1 = ['posydon.binary_evol.simulationproperties', 'TimingHooks']
        absolute_import_1 = None
        kwargs_1 = {}
        import_2 = ['posydon.binary_evol.simulationproperties', 'StepNamesHooks']
        absolute_import_2 = None
        kwargs_2 = {}
        """
        file_path = os.path.join(tmp_path, "sim.ini")
        with open(file_path, "w") as f:
            f.write(ini_content)
        return file_path

    @fixture
    def binpop_ini(self, tmp_path):
        ini_content = """
        [BinaryPopulation_options]
        use_MPI = False
        metallicity = [0.02]
        number_of_binaries = 1
        temp_directory = 'tmp'

        [BinaryStar_output]
        extra_columns = {}
        only_select_columns = []
        scalar_names = []

        [SingleStar_1_output]
        include_S1 = False

        [SingleStar_2_output]
        include_S2 = False
        """
        file_path = os.path.join(tmp_path, "binpop.ini")
        with open(file_path, "w") as f:
            f.write(ini_content)
        return file_path

    @fixture
    def binpop_ini_mpi(self, tmp_path):
        ini_content = """
        [BinaryPopulation_options]
        use_MPI = True
        metallicity = [0.02]
        number_of_binaries = 1
        temp_directory = 'tmp'

        [BinaryStar_output]
        extra_columns = {}
        only_select_columns = []
        scalar_names = []

        [SingleStar_1_output]
        include_S1 = False

        [SingleStar_2_output]
        include_S2 = False
        """
        file_path = os.path.join(tmp_path, "binpop_mpi.ini")
        with open(file_path, "w") as f:
            f.write(ini_content)
        return file_path

    @fixture
    def binpop_ini_stars(self, tmp_path):
        ini_content = """
        [BinaryPopulation_options]
        use_MPI = False
        metallicity = [0.02]
        number_of_binaries = 1
        temp_directory = 'tmp'

        [BinaryStar_output]
        extra_columns = {}
        only_select_columns = []
        scalar_names = []

        [SingleStar_1_output]
        include_S1 = True
        only_select_columns = [
            'state',
            'mass',
            'log_R']

        [SingleStar_2_output]
        include_S2 = True
        only_select_columns = [
            'log_L',
            'lg_mdot']
        """
        file_path = os.path.join(tmp_path, "binpop_stars.ini")
        with open(file_path, "w") as f:
            f.write(ini_content)
        return file_path

    @fixture
    def history_df(self):
        data = {
            'state': ['disrupted'],
            'time': [1.23],
            'S1_mass': [10.0],
            'S2_spin': [0.3]
        }
        return pd.DataFrame(data)

    @fixture
    def oneline_df(self):
        data = {
            'state_i': ['detached', 'detached'],
            'state_f': ['contact', 'merged'],
            'mass_i': [1.4, 2.1],
            'mass_f': [1.3, 2.0],
            'S1_spin_i': [0.5, 0.6],
            'S1_spin_f': [0.7, 0.8],
            'S1_SN_type': ['CCSN', 'NaN'],
            'S2_mass_i': [5.0, 6.0],
            'S2_mass_f': [7.0, 8.0],
            'S2_kick': [123.0, 456.0],
        }
        df = pd.DataFrame(data)
        return df

    def test_clean_binary_history_df(self, history_df):
        extra_binary = {'extra_binary': 'int32'}
        extra_S1 = {}
        extra_S2 = {}

        clean_df = totest.clean_binary_history_df(
            history_df,
            extra_binary_dtypes_user=extra_binary,
            extra_S1_dtypes_user=extra_S1,
            extra_S2_dtypes_user=extra_S2
        )
        assert isinstance(clean_df, pd.DataFrame)
        assert clean_df.dtypes['time'] == np.dtype('float64')
        assert clean_df.dtypes['S1_mass'] == np.dtype('float64')
        assert clean_df.dtypes['S2_spin'] == np.dtype('float64')
        assert clean_df.dtypes['state'] == np.dtype('O')

    def test_clean_binary_oneline_df(self, oneline_df):
        cleaned_df = totest.clean_binary_oneline_df(oneline_df)
        assert isinstance(cleaned_df, pd.DataFrame)
        assert cleaned_df['mass_i'].dtype == np.float64
        assert cleaned_df['S1_spin_i'].dtype == np.float64
        assert cleaned_df['state_i'].dtype == 'object'
        assert cleaned_df['state_f'].dtype == 'object'
        assert cleaned_df['S1_SN_type'].dtype == 'object'
        assert cleaned_df['S2_kick'].dtype == np.float64
        assert cleaned_df.loc[0, 'mass_i'] == 1.4
        assert cleaned_df.loc[1, 'state_f'] == 'merged'

    def test_parse_inifile(self,simple_ini,multi_ini,textfile):
        # missing argument
        with raises(TypeError, match="missing 1 required positional argument: 'path'"):
            totest.parse_inifile()
        # bad input
        with raises(FileNotFoundError):
            totest.parse_inifile('nonexistent.ini')
        with raises(FileNotFoundError):
            totest.parse_inifile([simple_ini,'nonexistent.ini'])
        with raises(MissingSectionHeaderError, match="File contains no section headers"):
            totest.parse_inifile(textfile)
        with raises(ValueError, match="Path must be a string or list of strings."):
            totest.parse_inifile(0)

        # example: single inifile
        parser = totest.parse_inifile(simple_ini)
        assert isinstance(parser, ConfigParser)
        assert parser.has_section("section")
        assert parser.get("section", "key") == "value"

        # example: multiple inifiles
        parser = totest.parse_inifile(multi_ini)
        assert parser.has_option("section", "key1")
        assert parser.has_option("section", "key2")


    def test_simprop_kwargs_from_ini(self,monkeypatch,sim_ini,tmp_path):
        # example
        dummy_cls = type('DummyClass', (), {})()

        # Patch importlib.import_module to return dummy modules with dummy classes
        def dummy_import_module(name, package=None):
            class DummyModule:
                pass
            setattr(DummyModule, 'TimingHooks', dummy_cls)
            setattr(DummyModule, 'StepNamesHooks', dummy_cls)
            setattr(DummyModule, 'flow_chart', dummy_cls)
            setattr(DummyModule, 'MS_MS_step', dummy_cls)
            return DummyModule()

        monkeypatch.setattr(importlib, "import_module", dummy_import_module)

        simkwargs = totest.simprop_kwargs_from_ini(sim_ini)

        # Check keys exist
        assert 'flow' in simkwargs
        assert 'step_HMS_HMS' in simkwargs
        assert 'extra_hooks' in simkwargs

        # Check classes mapped to dummy_cls
        assert simkwargs['flow'][0] is dummy_cls
        assert simkwargs['step_HMS_HMS'][0] is dummy_cls

        # extra_hooks is a list of tuples (class, kwargs)
        hooks = simkwargs['extra_hooks']
        assert isinstance(hooks, list)
        assert hooks[0][0] is dummy_cls
        assert hooks[0][1] == {}
        assert hooks[1][0] is dummy_cls
        assert hooks[1][1] == {}

        # test with 'only' parameter
        simkwargs_only = totest.simprop_kwargs_from_ini(sim_ini, only='step_HMS_HMS')
        assert 'step_HMS_HMS' in simkwargs_only
        assert 'flow' not in simkwargs_only

        # absolute imports
        dummy_code = """
class MyDummyClass:
    def __init__(self):
        self.value = 42
"""
        dummy_path = os.path.join(tmp_path, "dummy.py")
        with open(dummy_path, "w") as f:
            f.write(dummy_code)
        ini_content = f"""
        [flow]
        import = ['builtins', 'int']
        absolute_import = ['{dummy_path}', 'MyDummyClass']
        """
        ini_path = os.path.join(tmp_path, "sim_abs_import.ini")
        with open(ini_path, "w") as f:
            f.write(ini_content)
        simkwargs = totest.simprop_kwargs_from_ini(str(ini_path))
        dummy_class = simkwargs['flow'][0]
        assert dummy_class.__name__ == "MyDummyClass"
        instance = dummy_class()
        assert instance.value == 42


    def test_binarypop_kwargs_from_ini(self,monkeypatch,binpop_ini,
                                       binpop_ini_mpi,binpop_ini_stars):
        # bad configuration: MPI and job array
        monkeypatch.setenv("SLURM_ARRAY_JOB_ID", "123")
        with raises(ValueError, match="MPI must be turned off for job arrays."):
            totest.binarypop_kwargs_from_ini(binpop_ini_mpi)

        # example: include S1 and S2
        monkeypatch.setenv("SLURM_ARRAY_JOB_ID", "456")
        monkeypatch.setenv("SLURM_ARRAY_TASK_ID", "4")
        monkeypatch.setenv("SLURM_ARRAY_TASK_MIN", "2")
        monkeypatch.setenv("SLURM_ARRAY_TASK_COUNT", "10")
        binkwargs = totest.binarypop_kwargs_from_ini(binpop_ini_stars)
        assert binkwargs["include_S1"] is True
        assert "only_select_columns" in binkwargs["S1_kwargs"]
        assert "S2_kwargs" in binkwargs
        assert "log_L" in binkwargs["S2_kwargs"]["only_select_columns"]

        # example: environment variables
        binkwargs = totest.binarypop_kwargs_from_ini(binpop_ini)
        assert binkwargs["JOB_ID"] == 456
        assert binkwargs["RANK"] == 2  # 4 - 2
        assert binkwargs["size"] == 10
        assert isinstance(binkwargs, dict)
        assert binkwargs["metallicity"] == [0.02]
        assert binkwargs["comm"] is None

        # example: no Job ID, no MPI
        monkeypatch.delenv('SLURM_ARRAY_JOB_ID', raising=False)
        monkeypatch.delenv('SLURM_ARRAY_TASK_ID', raising=False)
        monkeypatch.delenv('SLURM_ARRAY_TASK_MIN', raising=False)
        monkeypatch.delenv('SLURM_ARRAY_TASK_COUNT', raising=False)
        binkwargs = totest.binarypop_kwargs_from_ini(binpop_ini)
        assert binkwargs['RANK'] is None
        assert binkwargs['size'] is None
        assert binkwargs['comm'] is None
