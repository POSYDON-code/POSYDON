"""Unit tests of posydon/grids/io.py
"""

__authors__ = [
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>"
]

# import the module which will be tested
import posydon.grids.io as totest

# aliases
os = totest.os

from inspect import isclass, isroutine

# import other needed code for the tests, which is not already imported in the
# module you like to test
from pytest import fixture, raises, warns

from posydon.utils.posydonwarning import (
    InappropriateValueWarning,
    MissingFilesWarning,
    ReplaceValueWarning,
)


# define test classes collecting several test functions
class TestElements:
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = {'BINARY_OUTPUT_FILE', 'EEP_FILE_EXTENSIONS', 'GridReader',\
                    'POSYDON_FORMAT_OPTIONS', 'Pwarn', 'RunReader',\
                    'SINGLE_OUTPUT_FILE', '__authors__', '__builtins__',\
                    '__cached__', '__doc__', '__file__', '__loader__',\
                    '__name__', '__package__', '__spec__', 'glob', 'gzip',\
                    'initial_values_from_dirname', 'os',\
                    'print_meta_contents', 'read_MESA_data_file',\
                    'read_initial_values'}
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

    def test_instance_POSYDON_FORMAT_OPTIONS(self):
        assert isinstance(totest.POSYDON_FORMAT_OPTIONS, dict),\
               "POSYDON_FORMAT_OPTIONS is of type: "\
               + str(type(totest.POSYDON_FORMAT_OPTIONS))

    def test_instance_BINARY_OUTPUT_FILE(self):
        assert isinstance(totest.BINARY_OUTPUT_FILE, str),\
               "BINARY_OUTPUT_FILE is of type: "\
               + str(type(totest.BINARY_OUTPUT_FILE))

    def test_instance_SINGLE_OUTPUT_FILE(self):
        assert isinstance(totest.SINGLE_OUTPUT_FILE, str),\
               "SINGLE_OUTPUT_FILE is of type: "\
               + str(type(totest.SINGLE_OUTPUT_FILE))

    def test_instance_EEP_FILE_EXTENSIONS(self):
        assert isinstance(totest.EEP_FILE_EXTENSIONS, list),\
               "EEP_FILE_EXTENSIONS is of type: "\
               + str(type(totest.EEP_FILE_EXTENSIONS))

    def test_instance_RunReader(self):
        assert isclass(totest.RunReader)

    def test_instance_GridReader(self):
        assert isclass(totest.GridReader)

    def test_instance_print_meta_contents(self):
        assert isroutine(totest.print_meta_contents)

    def test_instance_read_initial_values(self):
        assert isroutine(totest.read_initial_values)

    def test_instance_initial_values_from_dirname(self):
        assert isroutine(totest.initial_values_from_dirname)


class TestValues:
    # check that the values fit
    def test_value_POSYDON_FORMAT_OPTIONS(self):
        for k in ["ignored folders", "grid metadata", "run metadata"]:
            assert k in totest.POSYDON_FORMAT_OPTIONS.keys()
        for o in totest.POSYDON_FORMAT_OPTIONS:
            assert isinstance(totest.POSYDON_FORMAT_OPTIONS[o], list)

    def test_value_BINARY_OUTPUT_FILE(self):
        assert totest.BINARY_OUTPUT_FILE == "out.txt"

    def test_value_SINGLE_OUTPUT_FILE(self):
        assert totest.SINGLE_OUTPUT_FILE == "out_star1_formation_step0.txt"

    def test_value_EEP_FILE_EXTENSIONS(self):
        assert len(totest.EEP_FILE_EXTENSIONS) > 0
        assert ".data.eep" in totest.EEP_FILE_EXTENSIONS


class TestFunctions:
    @fixture
    def grid_csv_files(self, tmp_path):
        # a temporary grid.csv file for testing
        path = os.path.join(tmp_path, "grid.csv")
        with open(path, "w") as test_file:
            test_file.write("unit,test\n")
            test_file.write("1.2345678901234567890123456789012345678901234"\
                            +"5678901234567890123456789,77777777.88888888\n")
            for i in range(3,6):
                test_file.write(f"{i}.0,{i}.0\n")
        # additionally create a zipped file
        try:
            os.system(f"gzip -1 -k {path}")
        except:
            raise RuntimeError("Please check that you have `gzip` installed "\
                               +"and up to date.")
        return path

    @fixture
    def inlist_grid_points_files(self, tmp_path):
        # a temporary inlist_grid_points file for testing
        path = os.path.join(tmp_path, "inlist_grid_points")
        with open(path, "w") as test_file:
            test_file.write("Unit = TEST\n")
        try:
            os.system(f"gzip -1 {path}")
        except:
            raise RuntimeError("Please check that you have `gzip` installed "\
                               +"and up to date.")
        # second, non-zipped file
        with open(path, "w") as test_file:
            test_file.write("Unittest\n")
            test_file.write("Test = TEST = test\n")
            test_file.write("m1 = 1.0\n")
            test_file.write("m2 = 0.2\n")
            test_file.write("initial_period_in_days = 1.0d+3\n")
            test_file.write("initial_z = 0.01\n")
        return path

    # test functions
    def test_print_meta_contents(self, grid_csv_files, capsys):
        # missing argument
        with raises(TypeError, match="missing 1 required positional "\
                                     +"argument: 'path'"):
            totest.print_meta_contents()
        # bad input
        with raises(AttributeError, match="'NoneType' object has no "\
                                          +"attribute 'endswith'"):
            totest.print_meta_contents(None)
        # examples
        totest.print_meta_contents(grid_csv_files)
        captured_output = capsys.readouterr()
        assert "CONTENTS OF METADATA FILE:" in captured_output.out
        assert 80 * "-" in captured_output.out
        assert "unit,test" in captured_output.out
        assert "1.23456789" + 6 * "0123456789" + ",77777777." in\
               captured_output.out
        assert ".8" not in captured_output.out
        for i in range(3,6):
            assert f"{i}.0,{i}.0" in captured_output.out
        totest.print_meta_contents(grid_csv_files+".gz", max_lines=0)
        captured_output = capsys.readouterr()
        assert captured_output.out == ""
        totest.print_meta_contents(grid_csv_files+".gz", max_lines=2,\
                                   max_chars_per_line="wrap")
        captured_output = capsys.readouterr()
        assert "unit,test" in captured_output.out
        assert "77777777.88888888" in captured_output.out
        assert "3.0,3.0" not in captured_output.out

    def test_read_initial_values(self, tmp_path, inlist_grid_points_files):
        # missing argument
        with raises(TypeError, match="missing 1 required positional "\
                                     +"argument: 'mesa_dir'"):
            totest.read_initial_values()
        # bad input
        with raises(TypeError, match="expected str, bytes or os.PathLike "\
                                     +"object, not NoneType"):
            totest.read_initial_values(None)
        # examples
        assert totest.read_initial_values("./") is None
        with warns(InappropriateValueWarning, match="Multiple `=`, skipping "\
                                                    +"line in"):
            iv = totest.read_initial_values(tmp_path)
            assert iv['star_1_mass'] == 1.0
            assert iv['star_2_mass'] == 0.2
            assert iv['period_days'] == 1.0e+3
            assert iv['initial_z'] == 0.01
        os.remove(inlist_grid_points_files)
        with raises(ValueError, match="could not convert string to float: "\
                                      +"'TEST'"):
            totest.read_initial_values(tmp_path)

    def test_initial_values_from_dirname(self):
        # missing argument
        with raises(TypeError, match="missing 1 required positional "\
                                     +"argument: 'mesa_dir'"):
            totest.initial_values_from_dirname()
        # bad input
        with raises(TypeError, match="expected str, bytes or os.PathLike "\
                                     +"object, not NoneType"):
            totest.initial_values_from_dirname(None)
        with raises(AssertionError):
            totest.initial_values_from_dirname("./UnitTest")
        with raises(IndexError, match="list index out of range"):
            totest.initial_values_from_dirname("./initial_mass_1.0")
        # examples: single star v1
        assert totest.initial_values_from_dirname("./v1/initial_mass_1.0") ==\
               ('1.0',)
        # examples: single star v2
        assert totest.initial_values_from_dirname(\
               "./v2/initial_z_1.0_initial_mass_2.0") == ('2.0', '1.0')
        # examples: binary v1
        assert totest.initial_values_from_dirname(\
               "./v1/m1_1.1_m2_1.2_initial_period_in_days_1.0") ==\
               ('1.1', '1.2', '1.0')
        # examples: binary v2
        assert totest.initial_values_from_dirname("./v2/"\
               +"initial_z_1.0_m1_2.1_m2_2.2_initial_period_in_days_2.0") ==\
               ('2.1', '2.2', '2.0', '1.0')
        # examples: binary v2 with eccentricity
        assert totest.initial_values_from_dirname(\
               "initial_z_1.0_m1_2.1_m2_2.2_initial_period_in_days_2.0_initial_eccentricity_0.5") ==\
                ('2.1', '2.2', '2.0', '1.0', '0.5')

class TestRunReader:
    @fixture
    def RunReader(self, tmp_path):
        # initialize an instance of the class with defaults
        with warns(MissingFilesWarning, match="No relevant files found in"):
            return totest.RunReader(tmp_path)

    @fixture
    def Run_files(self, tmp_path):
        # a temporary files for testing
        for fn in ["binary_history.data", "final_star1.mod",\
                   "final_star2.mod", totest.BINARY_OUTPUT_FILE,\
                   totest.SINGLE_OUTPUT_FILE+".gz",\
                   totest.POSYDON_FORMAT_OPTIONS["run metadata"][0]]:
            file_path = os.path.join(tmp_path, fn)
            with open(file_path, "w") as test_file:
                test_file.write(f"Test: {file_path}\n")
        for d in ["LOGS1", "LOGS2"]:
            dir_path = os.path.join(tmp_path, d)
            os.mkdir(dir_path)
            for fn in ["history.data.gz", "final_profile.data",\
                       "initial_profile.data"]:
                file_path = os.path.join(dir_path, fn)
                with open(file_path, "w") as test_file:
                    test_file.write(f"Test: {file_path}\n")
        dir_path = os.path.join(tmp_path, "LOGS_files")
        os.mkdir(dir_path)
        for fn in ["LOGS1", "LOGS2"]:
            file_path = os.path.join(dir_path, fn)
            with open(file_path, "w") as test_file:
                test_file.write(f"Test: {file_path}\n")
        return

    # test the RunReader class
    def test_init(self, tmp_path, RunReader, capsys):
        assert isroutine(RunReader.__init__)
        # check that the instance is of correct type and all code in the
        # __init__ got executed: the elements are created and initialized
        assert isinstance(RunReader, totest.RunReader)
        assert RunReader.fmt == "posydon"
        assert RunReader.path == tmp_path
        assert RunReader.verbose == False
        assert RunReader.verbose_maxlines == 10
        assert RunReader.binary == True
        assert RunReader.history1_path is None
        assert RunReader.history2_path is None
        assert RunReader.binary_history_path is None
        assert RunReader.final_profile1_path is None
        assert RunReader.final_profile2_path is None
        assert RunReader.initial_profile1_path is None
        assert RunReader.initial_profile2_path is None
        assert RunReader.final_star1_path is None
        assert RunReader.final_star2_path is None
        assert RunReader.out_txt_path is None
        assert RunReader.metadata_files == []
        # missing argument
        with raises(TypeError, match="missing 1 required positional "\
                                     +"argument: 'path'"):
            test_RunReader = totest.RunReader()
        # bad input
        with raises(ValueError, match="Format test not supported."):
            test_RunReader = totest.RunReader("./", fmt="test")
        # examples
        assert capsys.readouterr().out == ""
        test_RunReader = totest.RunReader("./", fmt="posydon_single",\
                                          binary=False, verbose=True,\
                                          verbose_maxlines=1)
        assert test_RunReader.fmt == "posydon_single"
        assert test_RunReader.path == "./"
        assert test_RunReader.verbose == True
        assert test_RunReader.verbose_maxlines == 1
        assert test_RunReader.binary == False

    def test_read_posydon_format(self, RunReader, Run_files, capsys):
        assert isroutine(RunReader.read_posydon_format)
        # examples: default
        RunReader.read_posydon_format()
        assert capsys.readouterr().out == ""
        # examples: single
        RunReader.metadata_files = []
        RunReader.binary = False
        RunReader.verbose = True
        RunReader.read_posydon_format()
        captured_output = capsys.readouterr()
        assert "Reading run in {}".format(RunReader.path) in captured_output.out
        assert captured_output.out.count(": Yes") == 10
        assert "METADATA FILES: 2" in captured_output.out
        # examples: LOGS? are files instead of directories
        RunReader.path = os.path.join(RunReader.path, "LOGS_files")
        RunReader.verbose = False
        RunReader.read_posydon_format()

    def test_report(self, RunReader, capsys):
        assert isroutine(RunReader.report)
        RunReader.report()
        captured_output = capsys.readouterr()
        assert 80 * "-" in captured_output.out
        assert "DATA FOUND" in captured_output.out
        for l in ["MESA screen output", "History of Star 1",\
                  "History of Star 2", "Binary history",\
                  "Final profile of Star 1", "Final profile of Star 2",\
                  "Initial profile of Star 1", "Initial profile of Star 2",\
                  "Final model of Star 1", "Final model of Star 2"]:
            assert l + (25-len(l)) * " " + ": No" in captured_output.out
        assert "METADATA FILES:" in captured_output.out
        RunReader.metadata_files.append("./Test.csv")
        with raises(FileNotFoundError):
            RunReader.report()


class TestGridReader:
    @fixture
    def Run_path(self, tmp_path):
        # a temporary path for testing
        return os.path.join(tmp_path, "test_run_initial_z_0.01_m1_1.0"\
                                      +"_m2_0.2_initial_period_in_days_0.1")

    @fixture
    def Grid_files(self, Run_path):
        # a temporary files for testing
        path = Run_path
        os.mkdir(path)
        for fn in [totest.BINARY_OUTPUT_FILE]:
            file_path = os.path.join(path, fn)
            with open(file_path, "w") as test_file:
                test_file.write(f"Test: {file_path}\n")
        path += "_index_0"
        os.mkdir(path)
        for fn in [totest.BINARY_OUTPUT_FILE]:
            file_path = os.path.join(path, fn)
            with open(file_path, "w") as test_file:
                test_file.write(f"Test: {file_path}\n")
        return path

    @fixture
    def GridReader(self, tmp_path, Grid_files):
        # initialize an instance of the class with defaults
        with warns(ReplaceValueWarning, match="substitutes run in"):
            return totest.GridReader(str(tmp_path))

    @fixture
    def Empty_dir(self, tmp_path):
        dir_path = os.path.join(tmp_path, "empty")
        os.mkdir(dir_path)
        return dir_path

    # test the GridReader class
    def test_init(self, tmp_path, Run_path, GridReader):
        assert GridReader.fmt == "posydon"
        assert GridReader.path == str(tmp_path)
        assert GridReader.verbose == False
        assert GridReader.verbose_maxlines == 10
        assert GridReader.metadata_files == []
        assert isinstance(GridReader.runs, list)
        assert GridReader.binary == True
        # it depend on the system, which of the two directories will get read
        # first and therefore overwritten by the other
        assert (Run_path in GridReader.input_folders) or\
               (Run_path+"_index_0" in GridReader.input_folders)
        # missing argument
        with raises(TypeError, match="missing 1 required positional "\
                                     +"argument: 'path'"):
            test_GridReader = totest.GridReader()
        # bad input
        with raises(ValueError, match="Format test not supported."):
            test_GridReader = totest.GridReader(str(tmp_path), fmt="test")

    def test_get_input_folders(self, Empty_dir, GridReader, capsys):
        assert isroutine(GridReader._get_input_folders)
        GridReader.verbose = True
        GridReader.path = (GridReader.path,)
        input_folders = GridReader._get_input_folders()
        for f in GridReader.input_folders:
            assert f in input_folders
        assert GridReader.metadata_files == []
        captured_output = capsys.readouterr()
        for p in GridReader.path:
            assert "Searching for MESA runs in `{}`".format(p) in\
                   captured_output.out
        GridReader.verbose = False
        # bad input
        GridReader.path = Empty_dir
        with raises(ValueError, match="No folders found in "\
                                      +"{}".format(Empty_dir)):
            GridReader._get_input_folders()
        GridReader.path = None
        with raises(ValueError, match="MESA path must be a string or list"):
            GridReader._get_input_folders()

    def test_read_posydon_format(self, GridReader, Empty_dir, capsys):
        assert isroutine(GridReader.read_posydon_format)
        GridReader.verbose = True
        GridReader.read_posydon_format()
        assert len(GridReader.runs) == 2
        assert isinstance(GridReader.runs[0], totest.RunReader)
        captured_output = capsys.readouterr()
        assert "Reading grid in {}".format(GridReader.path) in\
               captured_output.out
        GridReader.verbose = False
        GridReader.input_folders = list(GridReader.input_folders)+[Empty_dir,\
                           totest.POSYDON_FORMAT_OPTIONS["ignored folders"][0]]
        with warns(MissingFilesWarning):
            GridReader.read_posydon_format()

    def test_infer_history_columns(self, GridReader):
        assert isroutine(GridReader.infer_history_columns)
        # missing argument
        with raises(TypeError, match="missing 3 required positional "\
                                     +"arguments: 'BH_cols', 'H1_cols', and "\
                                     +"'H2_cols'"):
            GridReader.infer_history_columns()
        # examples
        assert GridReader.infer_history_columns([], [], []) == []
        GridReader.runs.append(GridReader.runs[0])
        for bh in [None, "unit.test"]:
            for h1 in [None, "unit.test"]:
                for h2 in [None, "unit.test"]:
                    GridReader.runs[0].binary_history_path = bh
                    GridReader.runs[0].history1_path = h1
                    GridReader.runs[0].history2_path = h2
                    assert GridReader.infer_history_columns([], [], []) == []
