"""Unit tests of posydon/utils/gridutils.py
"""

__authors__ = [
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>"
]

# import the module which will be tested
import posydon.utils.gridutils as totest

# import other needed code for the tests, which is not already imported in the
# module you like to test
from pytest import fixture, raises, warns, approx
from inspect import isclass, isroutine
from posydon.utils.posydonwarning import MissingFilesWarning

# define test classes collecting several test functions
class TestElements:
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = ['LG_MTRANSFER_RATE_THRESHOLD', 'Msun', 'Pwarn', 'Rsun',\
                    'T_merger_P', 'T_merger_a', '__authors__', '__builtins__',\
                    '__cached__', '__doc__', '__file__', '__loader__',\
                    '__name__', '__package__', '__spec__', 'add_field',\
                    'beta_gw', 'cgrav', 'clean_inlist_file', 'clight',\
                    'convert_output_to_table', 'find_index_nearest_neighbour',\
                    'find_nearest', 'fix_He_core', 'get_cell_edges',\
                    'get_final_proposed_points', 'get_new_grid_name', 'gzip',\
                    'join_lists', 'kepler3_a', 'np', 'os', 'pd',\
                    'read_EEP_data_file', 'read_MESA_data_file', 'secyear']
        assert dir(totest) == elements, "There might be added or removed "+\
               "objects without an update on the unit test."

    def test_instance_join_lists(self):
        assert isroutine(totest.join_lists)

    def test_instance_read_MESA_data_file(self):
        assert isroutine(totest.read_MESA_data_file)

    def test_instance_read_EEP_data_file(self):
        assert isroutine(totest.read_EEP_data_file)

    def test_instance_fix_He_core(self):
        assert isroutine(totest.fix_He_core)

    def test_instance_add_field(self):
        assert isroutine(totest.add_field)

    def test_instance_get_cell_edges(self):
        assert isroutine(totest.get_cell_edges)

    def test_instance_find_nearest(self):
        assert isroutine(totest.find_nearest)

    def test_instance_find_index_nearest_neighbour(self):
        assert isroutine(totest.find_index_nearest_neighbour)

    def test_instance_get_final_proposed_points(self):
        assert isroutine(totest.get_final_proposed_points)

    def test_instance_T_merger_P(self):
        assert isroutine(totest.T_merger_P)

    def test_instance_beta_gw(self):
        assert isroutine(totest.beta_gw)

    def test_instance_kepler3_a(self):
        assert isroutine(totest.kepler3_a)

    def test_instance_T_merger_a(self):
        assert isroutine(totest.T_merger_a)

    def test_instance_convert_output_to_table(self):
        assert isroutine(totest.convert_output_to_table)

    def test_instance_clean_inlist_file(self):
        assert isroutine(totest.clean_inlist_file)

    def test_instance_get_new_grid_name(self):
        assert isroutine(totest.get_new_grid_name)


class TestFunctions:
    @fixture
    def data_path(self, tmp_path):
        # a temporary data path for testing
        return totest.os.path.join(tmp_path, "history.data")

    @fixture
    def no_path(self, tmp_path):
        # a path which does not exist for testing
        return totest.os.path.join(tmp_path, "does_not_exist.test")

    @fixture
    def MESA_data(self):
        # mock data
        return totest.np.array([(1, 2, 3.3), (1, -2, -3.3)],\
               dtype=[('COL1', '<f8'), ('COL2', '<f8'), ('COL3', '<f8')])

    @fixture
    def MESA_data_file(self, data_path, MESA_data):
        # a temporary data file for testing
        with open(data_path, "w") as test_file:
            for i in range(5):
                test_file.write("Test HEADER{}\n".format(i+1))
            for k in MESA_data.dtype.names:
                test_file.write("{:<5}".format(k))
            for j in range(2):
                for i,v in enumerate(MESA_data[j]):
                    if v==int(v):
                        v = int(v)
                    if i==0:
                        test_file.write("\n{:>4}".format(v))
                    else:
                        test_file.write(" {:>4}".format(v))
        return

    @fixture
    def EEP_data_file(self, data_path, MESA_data):
        # a temporary data file for testing
        with open(data_path, "w") as test_file:
            for i in range(11):
                test_file.write("Test HEADER{}\n".format(i+1))
            for k in MESA_data.dtype.names:
                test_file.write("{:<5}".format(k))
            for j in range(2):
                for i,v in enumerate(MESA_data[j]):
                    if v==int(v):
                        v = int(v)
                    if i==0:
                        test_file.write("\n{:>4}".format(v))
                    else:
                        test_file.write(" {:>4}".format(v))
        return

    @fixture
    def MESA_BH_data(self):
        # mock data
        return totest.np.array([(1.0, 15.0, 30.0, -5.0),\
                                (1.0, 15.0, 30.0, -4.0)],\
               dtype=[('period_days', '<f8'), ('star_1_mass', '<f8'),\
                      ('star_2_mass', '<f8'), ('lg_mtransfer_rate', '<f8')])

    @fixture
    def BH_path(self, tmp_path):
        # a temporary binary history path for testing
        return totest.os.path.join(tmp_path, "binary_history.data")

    @fixture
    def BH_file(self, BH_path, MESA_BH_data):
        # a temporary binary history file for testing
        with open(BH_path, "w") as test_file:
            for i in range(5):
                test_file.write("Test HEADER{}\n".format(i+1))
            for k in MESA_BH_data.dtype.names:
                test_file.write("{:<18}".format(k))
            for j in range(2):
                for i,v in enumerate(MESA_BH_data[j]):
                    if v==int(v):
                        v = int(v)
                    if i==0:
                        test_file.write("\n{:>17}".format(v))
                    else:
                        test_file.write(" {:>17}".format(v))
        return

    @fixture
    def MESA_BH_data2(self):
        # mock data
        return totest.np.array([(100.0, 15.0, 30.0, -5.0),\
                                (100.0, 15.0, 30.0, -4.0)],\
               dtype=[('period_days', '<f8'), ('star_1_mass', '<f8'),\
                      ('star_2_mass', '<f8'), ('lg_mtransfer_rate', '<f8')])

    @fixture
    def BH_path2(self, tmp_path):
        # a temporary binary history path for testing
        return totest.os.path.join(tmp_path, "binary_history2.data")

    @fixture
    def BH_file2(self, BH_path2, MESA_BH_data2):
        # a temporary binary history file for testing
        with open(BH_path2, "w") as test_file:
            for i in range(5):
                test_file.write("Test HEADER{}\n".format(i+1))
            for k in MESA_BH_data2.dtype.names:
                test_file.write("{:<18}".format(k))
            for j in range(2):
                for i,v in enumerate(MESA_BH_data2[j]):
                    if v==int(v):
                        v = int(v)
                    if i==0:
                        test_file.write("\n{:>17}".format(v))
                    else:
                        test_file.write(" {:>17}".format(v))
        return

    @fixture
    def MESA_SH_data(self):
        # mock data
        return totest.np.array([(1.0, 4.0, 0.0, 0.0), (2.0, 4.0, 5.0, 4.0)],\
               dtype=[('log_L', '<f8'), ('log_Teff', '<f8'),\
                      ('he_core_mass', '<f8'), ('c_core_mass', '<f8')])

    @fixture
    def H1_path(self, tmp_path):
        # a temporary star1 history path for testing
        return totest.os.path.join(tmp_path, "history1.data")

    @fixture
    def H1_file(self, H1_path, MESA_SH_data):
        # a temporary star1 history file for testing
        with open(H1_path, "w") as test_file:
            for i in range(5):
                test_file.write("Test HEADER{}\n".format(i+1))
            for k in MESA_SH_data.dtype.names:
                test_file.write("{:<13}".format(k))
            for j in range(2):
                for i,v in enumerate(MESA_SH_data[j]):
                    if v==int(v):
                        v = int(v)
                    if i==0:
                        test_file.write("\n{:>12}".format(v))
                    else:
                        test_file.write(" {:>12}".format(v))
        return

    @fixture
    def H2_path(self, tmp_path):
        # a temporary star2 history path for testing
        return totest.os.path.join(tmp_path, "history2.data")

    @fixture
    def H2_file(self, H2_path, MESA_SH_data):
        # a temporary star2 history file for testing
        with open(H2_path, "w") as test_file:
            for i in range(5):
                test_file.write("Test HEADER{}\n".format(i+1))
            for k in MESA_SH_data.dtype.names:
                test_file.write("{:<13}".format(k))
            for j in range(2):
                for i,v in enumerate(MESA_SH_data[j]):
                    if v==int(v):
                        v = int(v)
                    if i==0:
                        test_file.write("\n{:>12}".format(v))
                    else:
                        test_file.write(" {:>12}".format(v))
        return

    @fixture
    def out_path(self, tmp_path):
        # a temporary out path for testing
        return totest.os.path.join(tmp_path, "out.txt")

    @fixture
    def out_file(self, out_path):
        # a temporary out file for testing
        with open(out_path, "w") as test_file:
            for i in range(5):
                test_file.write("Test HEADER{}\n".format(i+1))
        return

    @fixture
    def out_gz_path(self, tmp_path):
        # a temporary out path for testing
        return totest.os.path.join(tmp_path, "out_gz.txt")

    @fixture
    def out_gz_file(self, out_gz_path):
        # a temporary out file for testing
        with open(out_gz_path, "w") as test_file:
            for i in range(5):
                test_file.write("Test HEADER{}\n".format(i+1))
            test_file.write("model is overflowing at ZAMS\n")
        totest.os.system(f"gzip -1 {out_gz_path}")
        return

    @fixture
    def out_path2(self, tmp_path):
        # a temporary out path for testing
        return totest.os.path.join(tmp_path, "out2.txt")

    @fixture
    def out_file2(self, out_path2):
        # a temporary out file for testing
        with open(out_path2, "w") as test_file:
            for i in range(5):
                test_file.write("Test HEADER{}\n".format(i+1))
            test_file.write("TURNING ON CE\n")
        return

    @fixture
    def out_path3(self, tmp_path):
        # a temporary out path for testing
        return totest.os.path.join(tmp_path, "out3.txt")

    @fixture
    def out_file3(self, out_path3):
        # a temporary out file for testing
        with open(out_path3, "w") as test_file:
            for i in range(5):
                test_file.write("Test HEADER{}\n".format(i+1))
            test_file.write("Terminate because accretor (r-rl)/rl > "+\
                            "accretor_overflow_terminate\n")
        return

    @fixture
    def out_path4(self, tmp_path):
        # a temporary out path for testing
        return totest.os.path.join(tmp_path, "out4.txt")

    @fixture
    def out_file4(self, out_path4):
        # a temporary out file for testing
        with open(out_path4, "w") as test_file:
            for i in range(5):
                test_file.write("Test HEADER{}\n".format(i+1))
            test_file.write("TURNING ON CE\n")
            test_file.write("Terminate because accretor (r-rl)/rl > "+\
                            "accretor_overflow_terminate\n")
        return

    @fixture
    def out_path5(self, tmp_path):
        # a temporary out path for testing
        return totest.os.path.join(tmp_path, "out5.txt")

    @fixture
    def out_file5(self, out_path5):
        # a temporary out file for testing
        with open(out_path5, "w") as test_file:
            for i in range(5):
                test_file.write("Test HEADER{}\n".format(i+1))
            test_file.write("Terminate due to primary depleting carbon\n")
        return

    @fixture
    def out_path6(self, tmp_path):
        # a temporary out path for testing
        return totest.os.path.join(tmp_path, "out6.txt")

    @fixture
    def out_file6(self, out_path6):
        # a temporary out file for testing
        with open(out_path6, "w") as test_file:
            for i in range(5):
                test_file.write("Test HEADER{}\n".format(i+1))
            test_file.write("TURNING ON CE\n")
            test_file.write("Terminate due to primary depleting carbon\n")
        return

    @fixture
    def ini_path(self, tmp_path):
        # a temporary ini path for testing
        return totest.os.path.join(tmp_path, "test.ini")

    @fixture
    def ini_file(self, ini_path):
        # a temporary ini file for testing
        with open(ini_path, "w") as test_file:
            test_file.write("[TEST INI SECTION]\n")
            test_file.write("test_float = 1.0\n")
            test_file.write("test_string = 'unit'\n")
            test_file.write("!test_comment = 0\n")
        return

    @fixture
    def ini_path2(self, tmp_path):
        # a temporary ini path for testing
        return totest.os.path.join(tmp_path, "test2.ini")

    @fixture
    def ini_file2(self, ini_path2):
        # a temporary ini file for testing
        with open(ini_path2, "w") as test_file:
            test_file.write("TEST INI SECTION 1&\n")
            test_file.write("test_int = 1\n")
            test_file.write("TEST INI SECTION 2&\n")
            test_file.write("test_empty = ''\n")
        return

    # test functions
    def test_join_lists(self):
        # missing argument
        with raises(TypeError, match="missing 2 required positional"+\
                    " arguments: 'A' and 'B'"):
            totest.join_lists()
        # bad input
        with raises(TypeError, match="'int' object is not iterable"):
            totest.join_lists(1, [])
        with raises(TypeError, match="'int' object is not iterable"):
            totest.join_lists([], 1)
        # test cases
        assert totest.join_lists([1, 2, 3], [4, 5]) == [1, 2, 3, 4, 5]
        assert totest.join_lists([1, 2, 3], [4, 1]) == [1, 2, 3, 4]
        assert totest.join_lists([1, 2, 1], [4, 1]) == [1, 2, 1, 4]

    def test_read_MESA_data_file(self, no_path, data_path,\
                                 MESA_data_file, MESA_data):
        # missing argument
        with raises(TypeError, match="missing 2 required positional"+\
                    " arguments: 'path' and 'columns'"):
            totest.read_MESA_data_file()
        # path is None or non-existing
        assert totest.read_MESA_data_file(None, ['Test']) is None
        assert totest.read_MESA_data_file(no_path, ['Test']) is None
        # read full test file
        assert totest.np.array_equal(MESA_data,\
               totest.read_MESA_data_file(data_path, MESA_data.dtype.names))
        # read columns individually from test file
        for k in MESA_data.dtype.names:
            assert totest.np.array_equal(MESA_data[[k]],\
                   totest.read_MESA_data_file(data_path, [k]))
        # read column pairs from test file
        k0 = MESA_data.dtype.names[0]
        for k in MESA_data.dtype.names:
            if k!=k0:
                assert totest.np.array_equal(MESA_data[[k0, k]],\
                       totest.read_MESA_data_file(data_path, [k0, k]))
        # warning for non readable data
        with open(data_path, "w") as test_file:
            test_file.write(10*"Test\n")
        with warns(MissingFilesWarning, match="Problems with reading file "+\
                   data_path):
            assert totest.read_MESA_data_file(data_path,\
                   MESA_data.dtype.names) is None

    def test_read_EEP_data_file(self, no_path, data_path,\
                                EEP_data_file, MESA_data):
        # missing argument
        with raises(TypeError, match="missing 2 required positional"+\
                    " arguments: 'path' and 'columns'"):
            totest.read_EEP_data_file()
        # path is None or non-existing
        assert totest.read_EEP_data_file(None, ['Test']) is None
        with warns(MissingFilesWarning, match="Problems with reading file "+\
                   no_path):
            assert totest.read_EEP_data_file(no_path, ['Test']) is None
        # read full test file
        assert totest.np.array_equal(MESA_data,\
               totest.read_EEP_data_file(data_path, MESA_data.dtype.names))
        # read columns individually from test file
        for k in MESA_data.dtype.names:
            assert totest.np.array_equal(MESA_data[[k]],\
                   totest.read_EEP_data_file(data_path, [k]))
        # read column pairs from test file
        k0 = MESA_data.dtype.names[0]
        for k in MESA_data.dtype.names:
            if k!=k0:
                assert totest.np.array_equal(MESA_data[[k0, k]],\
                       totest.read_EEP_data_file(data_path, [k0, k]))
        # warning for non readable data
        with open(data_path, "w") as test_file:
            test_file.write(15*"Test\n")
        with warns(match="Problems with reading file "+data_path):
            assert totest.read_EEP_data_file(data_path,\
                   MESA_data.dtype.names) is None

    def test_fix_He_core(self, MESA_data):
        # missing argument
        with raises(TypeError, match="missing 1 required positional"+\
                    " argument: 'history'"):
            totest.fix_He_core()
        # history is None
        assert totest.fix_He_core(None) is None
        # history is ndarray without the required columns
        assert totest.np.array_equal(MESA_data, totest.fix_He_core(MESA_data))
        # history is ndarray with the required columns and corrects them
        test_history = totest.np.array([(1, 2, 3, 4), (1, 2, 4, 3),\
                                        (2, 1, 3, 4), (2, 1, 4, 3)],\
                                       dtype=[('he_core_mass', '<f8'),\
                                              ('co_core_mass', '<f8'),\
                                              ('he_core_radius', '<f8'),\
                                              ('co_core_radius', '<f8')])
        fixed_history = totest.np.array([(2, 2, 4, 4), (2, 2, 4, 3),\
                                         (2, 1, 4, 4), (2, 1, 4, 3)],\
                                        dtype=[('he_core_mass', '<f8'),\
                                               ('co_core_mass', '<f8'),\
                                               ('he_core_radius', '<f8'),\
                                               ('co_core_radius', '<f8')])
        assert totest.np.array_equal(fixed_history,\
               totest.fix_He_core(test_history))

    def test_add_field(self, MESA_data):
        # missing argument
        with raises(TypeError, match="missing 2 required positional"+\
                    " arguments: 'a' and 'descr'"):
            totest.add_field()
        # add to empty ndarray
        with raises(ValueError, match="'a' must be a structured numpy array"):
            totest.add_field(totest.np.array([]), [('new', '<f8')])
        # add to test data
        extended_ndarray = totest.add_field(MESA_data, [('new', '<f8')])
        assert MESA_data.dtype.descr+[('new', '<f8')] ==\
               extended_ndarray.dtype.descr
        assert totest.np.array_equal(MESA_data,\
               extended_ndarray[[k for k in MESA_data.dtype.names]])

    def test_get_cell_edges(self):
        # missing argument
        with raises(TypeError, match="missing 2 required positional"+\
                    " arguments: 'grid_x' and 'grid_y'"):
            totest.get_cell_edges()
        # example for x = [0.1, 0.3, 0.5, 0.7, 0.9] and y = [0.1, 0.3, 0.5]
        assert totest.np.allclose(\
               totest.get_cell_edges(totest.np.linspace(0.1, 0.9, 5),\
                                     totest.np.linspace(0.1, 0.5, 3)),\
               totest.np.meshgrid(totest.np.linspace(0.0, 1.0, 6),\
                                  totest.np.linspace(0.0, 0.6, 4)))

    def test_find_nearest(self):
        # missing argument
        with raises(TypeError, match="missing 2 required positional"+\
                    " arguments: 'val' and 'array'"):
            totest.find_nearest()
        # example for [0.1, 0.3, 0.5, 0.7, 0.9]
        test_data = totest.np.linspace(0.1, 0.9, 5)
        for v in test_data:
            assert totest.find_nearest(v+0.1, test_data) == v

    def test_find_index_nearest_neighbour(self):
        # missing argument
        with raises(TypeError, match="missing 2 required positional"+\
                    " arguments: 'array' and 'value'"):
            totest.find_index_nearest_neighbour()
        # example for [0.1, 0.3, 0.5, 0.7, 0.9]
        test_data = totest.np.linspace(0.1, 0.9, 5)
        for i,v in enumerate(test_data):
            assert totest.find_index_nearest_neighbour(test_data, v+0.1) == i

    def test_get_final_proposed_points(self, capsys):
        # missing argument
        with raises(TypeError, match="missing 4 required positional"+\
                    " arguments: 'proposed_x', 'grid_x', 'proposed_y', and"+\
                    " 'grid_y'"):
            totest.get_final_proposed_points()
        # example for grid_x = [0.1, 0.3, 0.5], grid_y = [0.1, 0.3, 0.5],
        # proposed_x = [0.19, 0.29, 0.39, 0.49, 0.59], and
        # proposed_y = [0.21, 0.31, 0.41, 0.51, 0.61]
        mx,my = totest.get_final_proposed_points(\
                 totest.np.linspace(0.19, 0.59, 5),\
                 totest.np.linspace(0.1, 0.5, 3),\
                 totest.np.linspace(0.21, 0.61, 5),\
                 totest.np.linspace(0.1, 0.5, 3))
        assert totest.np.allclose(mx, totest.np.array([0.1, 0.3, 0.3, 0.5]))
        assert totest.np.allclose(my, totest.np.array([0.3, 0.3, 0.5, 0.5]))

    def test_T_merger_P(self):
        # missing argument
        with raises(TypeError, match="missing 3 required positional"+\
                    " arguments: 'P', 'm1', and 'm2'"):
            totest.T_merger_P()
        # examples
        assert totest.T_merger_P(1.0, 15.0, 30.0) ==\
               approx(0.37210532488, abs=6e-12)
        assert totest.T_merger_P(2.0, 15.0, 30.0) ==\
               approx(2.36272153666, abs=6e-12)
        assert totest.T_merger_P(1.0, 30.0, 30.0) ==\
               approx(0.20477745195, abs=6e-12)
        assert totest.T_merger_P(1.0, 15.0, 60.0) ==\
               approx(0.22058982311, abs=6e-12)

    def test_beta_gw(self):
        # missing argument
        with raises(TypeError, match="missing 2 required positional"+\
                    " arguments: 'm1' and 'm2'"):
            totest.beta_gw()
        # examples
        assert totest.beta_gw(15.0, 30.0) ==\
               approx(3.18232660295e-69, abs=6e-81)
        assert totest.beta_gw(30.0, 30.0) ==\
               approx(8.48620427454e-69, abs=6e-81)
        assert totest.beta_gw(30.0, 60.0) ==\
               approx(2.54586128236e-68, abs=6e-80)

    def test_kepler3_a(self):
        # missing argument
        with raises(TypeError, match="missing 3 required positional"+\
                    " arguments: 'P', 'm1', and 'm2'"):
            totest.kepler3_a()
        # examples
        assert totest.kepler3_a(1.0, 15.0, 30.0) ==\
               approx(14.9643417735, abs=6e-11)
        assert totest.kepler3_a(2.0, 15.0, 30.0) ==\
               approx(23.7544118733, abs=6e-11)
        assert totest.kepler3_a(1.0, 30.0, 30.0) ==\
               approx(16.4703892879, abs=6e-11)
        assert totest.kepler3_a(1.0, 15.0, 60.0) ==\
               approx(17.7421890201, abs=6e-11)

    def test_T_merger_a(self):
        # missing argument
        with raises(TypeError, match="missing 3 required positional"+\
                    " arguments: 'a', 'm1', and 'm2'"):
            totest.T_merger_a()
        # examples
        assert totest.T_merger_a(1.0, 15.0, 30.0) ==\
               approx(7.42053829341e-06, abs=6e-18)
        assert totest.T_merger_a(2.0, 15.0, 30.0) ==\
               approx(1.18728612695e-04, abs=6e-16)
        assert totest.T_merger_a(1.0, 30.0, 30.0) ==\
               approx(2.78270186003e-06, abs=6e-18)
        assert totest.T_merger_a(1.0, 15.0, 60.0) ==\
               approx(2.22616148802e-06, abs=6e-18)

    def test_convert_output_to_table(self, no_path, out_path, out_file,\
                                     MESA_SH_data, MESA_BH_data,\
                                     MESA_BH_data2, BH_path, BH_file,\
                                     BH_path2, BH_file2, H1_path, H1_file,\
                                     H2_path, H2_file, out_gz_path,\
                                     out_gz_file, out_path2, out_file2,\
                                     out_path3, out_file3, out_path4,\
                                     out_file4, out_path5, out_file5,\
                                     out_path6, out_file6):
        # missing argument
        with raises(TypeError, match="missing 1 required positional"+\
                    " argument: 'output_file'"):
            totest.convert_output_to_table()
        # bad input
        with raises(ValueError, match="File "+no_path+" does not exist"):
            totest.convert_output_to_table(no_path)
        totest.os.system("touch "+out_path+"0")
        with raises(ValueError, match="The output file does not have any "+\
                    "data."):
            totest.convert_output_to_table(out_path+"0")
        # check warnings for non specified files
        with warns(MissingFilesWarning, match="You have not supplied a "+\
                   "binary history file to parse. This will cause all "+\
                   "binary history columns to be dashes."):
            with warns(MissingFilesWarning, match="You have not supplied a "+\
                       "star1 history file to parse. This will cause all "+\
                       "star1 history columns to be dashes."):
                with warns(MissingFilesWarning, match="You have not supplied"+\
                           " a star2 history file to parse. This will cause "+\
                           "all star2 history columns to be dashes."):
                    assert totest.convert_output_to_table(out_path,\
                           column_names=["Test"]).equals(\
                           totest.pd.DataFrame({'CE_flag': [0],\
                           'result': ["error"], 'Test': ["-"]}))
        # check zipped out file + initial RLO
        assert totest.convert_output_to_table(out_gz_path+".gz",\
               binary_history_file=BH_path, star1_history_file=H1_path,\
               star2_history_file=H2_path, column_names=["Test"]).equals(\
                totest.pd.DataFrame({'result': ["ZAMS_RLOF"], 'Test': ["-"]}))
        # check CE + existing column
        assert totest.convert_output_to_table(out_path2,\
               binary_history_file=BH_path, star1_history_file=H1_path,\
               star2_history_file=H2_path, column_names=["Test",\
               "result"]).equals(totest.pd.DataFrame({'CE_flag': [1],\
                                 'result': ["error"], 'Test': ["-"]}))
        # check contact
        assert totest.convert_output_to_table(out_path3,\
               binary_history_file=BH_path, star1_history_file=H1_path,\
               star2_history_file=H2_path, column_names=["Test"]).equals(\
                totest.pd.DataFrame({'CE_flag': [0], 'result': ["contact"],\
                                     'Test': ["-"]}))
        # check CE merger
        assert totest.convert_output_to_table(out_path4,\
               binary_history_file=BH_path, star1_history_file=H1_path,\
               star2_history_file=H2_path, column_names=["Test"]).equals(\
                totest.pd.DataFrame({'CE_flag': [1], 'result': ["CE_merger"],\
                                     'Test': ["-"]}))
        # check carbon depletion: stable_MT_to_merger
        test_table = totest.convert_output_to_table(out_path5,\
                     binary_history_file=BH_path, star1_history_file=H1_path,\
                     star2_history_file=H2_path, column_names=["Test"])
        assert test_table.at[0,'CE_flag']==0
        assert test_table.at[0,'result']=="stable_MT_to_merger"
        assert test_table.at[0,'Test']=="-"
        for i in [1,2]:
            assert test_table.at[0,f'log_L_{i}']==MESA_SH_data['log_L'][-1]
            assert test_table.at[0,f'log_T_{i}']==MESA_SH_data['log_Teff'][-1]
            assert test_table.at[0,f'He_core_{i}(Msun)']==\
                   MESA_SH_data['he_core_mass'][-1]
            assert test_table.at[0,f'C_core_{i}(Msun)']==\
                   MESA_SH_data['c_core_mass'][-1]
        assert test_table.at[0,'M_1f(Msun)']==MESA_BH_data['star_1_mass'][-1]
        assert test_table.at[0,'M_2f(Msun)']==MESA_BH_data['star_2_mass'][-1]
        assert test_table.at[0,'Porb_f(d)']==MESA_BH_data['period_days'][-1]
        assert 'tmerge(Gyr)' in test_table.columns
        # check carbon depletion: no_interaction
        test_table = totest.convert_output_to_table(out_path6,\
                     column_names=["Test"])
        assert test_table.at[0,'CE_flag']==1
        assert test_table.at[0,'result']=="no_interaction"
        assert test_table.at[0,'Test']=="-"
        for i in [1,2]:
            assert f'log_L_{i}' not in test_table.columns
            assert f'log_T_{i}' not in test_table.columns
            assert f'He_core_{i}(Msun)' not in test_table.columns
            assert f'C_core_{i}(Msun)' not in test_table.columns
        assert 'M_1f(Msun)' not in test_table.columns
        assert 'M_2f(Msun)' not in test_table.columns
        assert 'Porb_f(d)' not in test_table.columns
        assert 'tmerge(Gyr)' not in test_table.columns
        # check carbon depletion: CE_ejection_m
        test_table = totest.convert_output_to_table(out_path6,\
                     binary_history_file=BH_path, column_names=["Test"])
        assert test_table.at[0,'CE_flag']==1
        assert test_table.at[0,'result']=="CE_ejection_m"
        assert test_table.at[0,'Test']=="-"
        assert test_table.at[0,'M_1f(Msun)']==MESA_BH_data['star_1_mass'][-1]
        assert test_table.at[0,'M_2f(Msun)']==MESA_BH_data['star_2_mass'][-1]
        assert test_table.at[0,'Porb_f(d)']==MESA_BH_data['period_days'][-1]
        assert 'tmerge(Gyr)' in test_table.columns
        # check carbon depletion: stable_MT_to_wide_binary
        test_table = totest.convert_output_to_table(out_path5,\
                     binary_history_file=BH_path2, column_names=["Test"])
        assert test_table.at[0,'CE_flag']==0
        assert test_table.at[0,'result']=="stable_MT_to_wide_binary"
        assert test_table.at[0,'Test']=="-"
        assert test_table.at[0,'M_1f(Msun)']==MESA_BH_data2['star_1_mass'][-1]
        assert test_table.at[0,'M_2f(Msun)']==MESA_BH_data2['star_2_mass'][-1]
        assert test_table.at[0,'Porb_f(d)']==MESA_BH_data2['period_days'][-1]
        assert 'tmerge(Gyr)' in test_table.columns
        # check carbon depletion: CE_ejection
        test_table = totest.convert_output_to_table(out_path6,\
                     binary_history_file=BH_path2, column_names=["Test"])
        assert test_table.at[0,'CE_flag']==1
        assert test_table.at[0,'result']=="CE_ejection"
        assert test_table.at[0,'Test']=="-"
        assert test_table.at[0,'M_1f(Msun)']==MESA_BH_data2['star_1_mass'][-1]
        assert test_table.at[0,'M_2f(Msun)']==MESA_BH_data2['star_2_mass'][-1]
        assert test_table.at[0,'Porb_f(d)']==MESA_BH_data2['period_days'][-1]
        assert 'tmerge(Gyr)' in test_table.columns

    def test_clean_inlist_file(self, ini_path, ini_file, ini_path2, ini_file2):
        # missing argument
        with raises(TypeError, match="missing 1 required positional"+\
                    " argument: 'inlist'"):
            totest.clean_inlist_file()
        assert totest.clean_inlist_file(ini_path)==\
               {None: {'test_float': '1.0', 'test_string': "'unit'"}}
        assert totest.clean_inlist_file(ini_path, section="test_section")==\
               {'test_section': {'test_float': '1.0', 'test_string': "'unit'"}}
        assert totest.clean_inlist_file(ini_path2)==\
               {'TEST INI SECTION 1&': {'test_int': '1'}, 'TEST INI SECTION 2&': {}}

    def test_get_new_grid_name(self, tmp_path):
        # missing argument
        with raises(TypeError, match="missing 2 required positional"+\
                    " arguments: 'path' and 'compression'"):
            totest.get_new_grid_name()
        # get name without check for existence of compression directory
        grid_dir = totest.os.path.join(tmp_path, "test_grid")
        compression_dir = totest.os.path.join(tmp_path, "TEST_COMPRESSION")
        assert totest.get_new_grid_name(grid_dir, "TEST_COMPRESSION")==\
               totest.os.path.join(compression_dir, "test_grid.h5")
        # create compression directory
        assert totest.get_new_grid_name(grid_dir, "TEST_COMPRESSION",\
                                        create_missing_directories=True)==\
               totest.os.path.join(compression_dir, "test_grid.h5")
        assert totest.os.path.isdir(compression_dir)
        # compression directory already exists, while it would be created
        # otherwise
        assert totest.get_new_grid_name(grid_dir, "TEST_COMPRESSION",\
                                        create_missing_directories=True)==\
               totest.os.path.join(compression_dir, "test_grid.h5")

