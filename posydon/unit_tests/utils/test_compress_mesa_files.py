"""Unit tests of posydon/utils/compress_mesa_files.py

"""

__authors__ = [
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>"
]

# import the module which will be tested
import posydon.utils.compress_mesa_files as totest

# aliases
os = totest.os

from inspect import isroutine

# import other needed code for the tests, which is not already imported in the
# module you like to test
from pytest import fixture, raises, warns

from posydon.unit_tests._helper_functions_for_tests.MESA import get_MESA_dir


# define test classes collecting several test functions
class TestElements:
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = {'Pwarn', '__authors__', '__builtins__', '__cached__',
                    '__doc__', '__file__', '__loader__', '__name__',
                    '__package__', '__spec__', '_compress_MESA',
                    '_parse_commandline', 'argparse', 'compress_dir', 'os',
                    'random', 'set_up_test', 'shutil', 'sys', 'textsize',
                    'tqdm', 'get_size'}
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

    def test_instance_parse_commandline(self):
        assert isroutine(totest._parse_commandline)

    def test_instance_textsize(self):
        assert isroutine(totest.textsize)

    def test_instance_set_up_test(self):
        assert isroutine(totest.set_up_test)

    def test_instance_get_size(self):
        assert isroutine(totest.get_size)

    def test_instance_compress_dir(self):
        assert isroutine(totest.compress_dir)

    def test_instance_compress_MESA(self):
        assert isroutine(totest._compress_MESA)


class TestFunctions:
    # test functions
    def test_parse_commandline(self, monkeypatch, capsys):
        with monkeypatch.context() as mp:
            # call help
            mp.setattr(totest.sys, "argv", ['compress-mesa', '--help'])
            with raises(SystemExit, match="0"):
                totest._parse_commandline()
            captured_out = capsys.readouterr().out
            assert "usage: compress-mesa [-h] [-td TEST_DIR] [-dsr DSR] [-v] "\
                   + "[-d] mesa_dir" in captured_out
            assert "Compressing MESA files" in captured_out
            assert "positional arguments:\n  mesa_dir" in captured_out
            assert "options:\n  -h, --help" in captured_out
            assert "-td TEST_DIR, --test_dir TEST_DIR" in captured_out
            assert "-dsr DSR, --dsr DSR" in captured_out
            assert "-v, --verbose" in captured_out
            assert "-d, --debug" in captured_out
            # missing argument
            tests = [(['compress-mesa-fail'], "the following arguments are "\
                                              +"required: mesa_dir"),\
                     (['compress-mesa-fail', '--test_dir', 'test'],\
                      "the following arguments are required: mesa_dir"),\
                     (['compress-mesa-fail', 'mesa_dir', '--test_dir'],\
                      "argument -td/--test_dir: expected one argument"),\
                     (['compress-mesa-fail', 'mesa_dir', '--dsr'],\
                      "argument -dsr/--dsr: expected one argument")]
            for (t, m) in tests:
                mp.setattr(totest.sys, "argv", t)
                with raises(SystemExit, match="2"):
                    totest._parse_commandline()
                assert m in capsys.readouterr().err
            # examples
            for mesa_dir in ['.', 'unit']:
                for test_dir in [None, 'test']:
                    for dsr in [0.01, 0.1]:
                        for verbose in [False, True]:
                            for debug in [False, True]:
                                commandline_args = ['compress-mesa-test']
                                if mesa_dir is not None:
                                    commandline_args += [str(mesa_dir)]
                                if test_dir is not None:
                                    commandline_args += ['--test_dir',\
                                                         str(test_dir)]
                                if dsr != 0.01:
                                    commandline_args += ['--dsr', str(dsr)]
                                if verbose:
                                    commandline_args += ['--verbose']
                                if debug:
                                    commandline_args += ['--debug']
                                mp.setattr(totest.sys, "argv",\
                                           commandline_args)
                                assert totest._parse_commandline()\
                                       == totest.argparse.Namespace(\
                                           mesa_dir=mesa_dir,\
                                           test_dir=test_dir, dsr=dsr,\
                                           verbose=verbose, debug=debug)

    def test_textsize(self):
        # missing argument
        with raises(TypeError, match="missing 1 required positional "\
                                     +"argument: 'filesize'"):
            totest.textsize()
        # bad input
        with raises(ValueError, match="base=1 should be 1000 or 1024"):
            totest.textsize(1, base=1)
        # bad input
        for t in [-1, 0]:
            with raises(ValueError, match=f"threshold={t} should be larger "\
                                          +"than 0"):
                totest.textsize(1, threshold=t)
        # bad input
        for (t, b) in [(1001, 1000), (1025, 1024)]:
            with raises(ValueError, match=f"threshold={t} should be smaller "\
                                          +f"or equal to base={b}"):
                totest.textsize(1, base=b, threshold=t)
        # examples
        tests = [(1, "1b"), (10, "10b"), (1e+2, "100b"), (1e+3, "0.977K"),\
                 (1e+4, "9.77K"), (1e+5, "97.7K"), (1e+6, "977K"),\
                 (1e+7, "9.54M"), (1e+8, "95.4M"), (1e+9, "954M"),\
                 (1e+10, "9.31G"), (1e+11, "93.1G"), (1e+12, "931G"),\
                 (1e+13, "9.09T"), (1e+14, "90.9T"), (1e+15, "909T"),\
                 (1e+16, "8.88P"), (1e+17, "88.8P"), (1e+18, "888P"),\
                 (1e+19, "8.67E"), (1e+20, "86.7E"), (1e+21, "867E"),\
                 (1e+22, "8.47Z"), (1e+23, "84.7Z"), (1e+24, "847Z"),\
                 (1e+25, "8.27Y"), (1e+26, "82.7Y"), (1e+27, "827Y"),\
                 (1e+28, "1e+28 bytes"), (-1, "-1b")]
        for (fz, r) in tests:
            assert totest.textsize(fz) == r
        assert totest.textsize(1e+3, floatfmt=".4g") == "0.9766K"
        assert totest.textsize(1e+3, base=1000) == "1K"
        assert totest.textsize(1e+3, threshold=1024) == "1e+03b"

    def test_set_up_test(self, tmp_path, capsys):
        # missing argument
        with raises(TypeError, match="missing 1 required positional "\
                                     +"argument: 'args'"):
            totest.set_up_test()
        # bad input
        with raises(NameError, match="--test_dir needs to be specified for "\
                                     +"set_up_test"):
            test_args = totest.argparse.Namespace(test_dir=None)
            totest.set_up_test(test_args)
        # bad input
        with raises(NotADirectoryError, match=f"Directory does_not_exist "\
                                              +"does not exist."):
            test_args = totest.argparse.Namespace(test_dir='does_not_exist')
            totest.set_up_test(test_args)
        # bad input
        with raises(NameError, match="mesa_dir needs to be specified for "\
                                     +"set_up_test"):
            test_args = totest.argparse.Namespace(mesa_dir=None, test_dir='.')
            totest.set_up_test(test_args)
        # bad input
        with raises(NotADirectoryError, match=f"Directory does_not_exist "\
                                              +"does not exist."):
            test_args = totest.argparse.Namespace(mesa_dir='does_not_exist',\
                                                  test_dir='.')
            totest.set_up_test(test_args)
        # examples: copy all MESA stuff (incl. extra stuff in batch directory)
        MESA_dir = get_MESA_dir(tmp_path, 1)
        os.mkdir(os.path.join(MESA_dir, "empty_dir")) # extra directory
        with open(os.path.join(MESA_dir, "test.txt"), "w") as test_file:
            test_file.write("Unit test:\nset_up_test\n") # extra file
        with open(os.path.join(tmp_path, "test.txt"), "w") as test_file:
            test_file.write("Unit test:\nset_up_test\n") # extra file
        test_dir = os.path.join(tmp_path, "test1")
        os.mkdir(test_dir)
        test_args = totest.argparse.Namespace(mesa_dir=tmp_path,\
                                              test_dir=test_dir, dsr=1.0)
        totest.set_up_test(test_args)
        assert capsys.readouterr().out == "Created Test Directory at "\
                                          + f"{test_dir}.\n"

    def test_get_size(self, tmp_path):
        # examples: empty directory
        MESA_dir = get_MESA_dir(tmp_path, 0)
        total_size, remove_files, compress_files, n_runs, n_remove_files,\
         n_compress_files = totest.get_size(start_path=MESA_dir)
        assert total_size == 0
        assert remove_files == []
        assert compress_files == []
        assert n_runs == 0
        assert n_remove_files == 0
        assert n_compress_files == 0
        # examples: binary runs
        MESA_dir = get_MESA_dir(tmp_path, 1, n_runs=20)
        total_size, remove_files, compress_files, n_runs, n_remove_files,\
         n_compress_files = totest.get_size(start_path=MESA_dir)
        assert total_size > 0
        assert remove_files == []
        assert len(compress_files) > 0
        assert n_runs == 20
        assert n_remove_files == 0
        assert n_compress_files > 0
        # examples: with core dump files to remove
        MESA_dir = get_MESA_dir(tmp_path, 2, n_runs=20)
        for i in [0, -1]:
            with open(os.path.join(MESA_dir, os.listdir(MESA_dir)[i],\
                                   "core.2"), "w") as core_dump_file:
                core_dump_file.write("Test\nremove\n")
        total_size, remove_files, compress_files, n_runs, n_remove_files,\
         n_compress_files = totest.get_size(start_path=MESA_dir)
        assert total_size > 0
        assert len(remove_files) == 2
        assert len(compress_files) > 0
        assert n_runs == 20
        assert n_remove_files == 2
        assert n_compress_files > 0
        # examples: with link
        MESA_dir = get_MESA_dir(tmp_path, 3, n_runs=20)
        MESA_runs = os.listdir(MESA_dir)
        if len(MESA_runs)>0:
            for i in [0, -1]:
                MESA_run_dir = os.path.join(MESA_dir, MESA_runs[i])
                os.symlink(MESA_run_dir, os.path.join(MESA_dir, f"link{i}"))
                MESA_run_file = os.path.join(MESA_run_dir,\
                                             os.listdir(MESA_run_dir)[0])
                os.symlink(MESA_run_file, os.path.join(MESA_dir,\
                                                       f"link{i}.file0"))
        total_size, remove_files, compress_files, n_runs, n_remove_files,\
         n_compress_files = totest.get_size(start_path=MESA_dir)
        assert total_size > 0
        assert len(remove_files) == 0
        assert len(compress_files) > 0
        assert n_runs == 20
        assert n_remove_files == 0
        assert n_compress_files > 0

    def test_compress_dir(self, tmp_path, capsys):
        # missing argument
        with raises(TypeError, match="missing 1 required positional "\
                                     +"argument: 'args'"):
            totest.compress_dir()
        # bad input
        with raises(NameError, match="mesa_dir needs to be specified for "\
                                     +"set_up_test"):
            test_args = totest.argparse.Namespace(mesa_dir=None)
            totest.compress_dir(test_args)
        # bad input
        with raises(NotADirectoryError, match=f"Directory does_not_exist "\
                                              +"does not exist."):
            test_args = totest.argparse.Namespace(mesa_dir='does_not_exist')
            totest.compress_dir(test_args)
        # examples: empty directory
        MESA_dir = get_MESA_dir(tmp_path, 0)
        test_args = totest.argparse.Namespace(mesa_dir=MESA_dir,\
                                              verbose=False, debug=False)
        totest.compress_dir(test_args)
        # examples: binary runs
        MESA_dir = get_MESA_dir(tmp_path, 1)
        test_args = totest.argparse.Namespace(mesa_dir=MESA_dir,\
                                              verbose=False, debug=False)
        totest.compress_dir(test_args)
        # examples: single star runs
        MESA_dir = get_MESA_dir(tmp_path, -1)
        test_args = totest.argparse.Namespace(mesa_dir=MESA_dir,\
                                              verbose=False, debug=False)
        totest.compress_dir(test_args)
        # examples: with verbose and debug output
        MESA_dir = get_MESA_dir(tmp_path, 2)
        test_args = totest.argparse.Namespace(mesa_dir=MESA_dir,\
                                              verbose=True, debug=True)
        totest.compress_dir(test_args)
        captured_out = capsys.readouterr().out.split('\n')
        assert "remove" in captured_out[0]
        assert "core dump files in" in captured_out[0]
        assert "directories of" in captured_out[0]
        assert "MESA runs" in captured_out[0]
        assert "compress" in captured_out[1]
        assert "files in" in captured_out[1]
        assert "directories of" in captured_out[1]
        assert "MESA runs" in captured_out[1]
        assert "compress:" in captured_out[2]
        assert "/MESA_data_index2/" in captured_out[2]
        assert "Compressed MESA tracks" in captured_out[-3]
        assert "Original size" in captured_out[-2]
        assert "| Compressed size" in captured_out[-2]
        assert "" == captured_out[-1]
        # examples: with core dump files to remove
        MESA_dir = get_MESA_dir(tmp_path, 3)
        for i in [0, -1]:
            with open(os.path.join(MESA_dir, os.listdir(MESA_dir)[i],\
                                   "core.3"), "w") as core_dump_file:
                core_dump_file.write("Test\nremove\n")
        test_args = totest.argparse.Namespace(mesa_dir=MESA_dir,\
                                              verbose=True, debug=True)
        totest.compress_dir(test_args)
        captured_out = capsys.readouterr().out.split('\n')
        assert "remove" in captured_out[0]
        assert "core dump files in" in captured_out[0]
        assert "directories of" in captured_out[0]
        assert "MESA runs" in captured_out[0]
        assert "remove:" in captured_out[1]
        assert "/MESA_data_index3/" in captured_out[1]
        assert "/core." in captured_out[1]
        assert "compress" in captured_out[3]
        assert "files in" in captured_out[3]
        assert "directories of" in captured_out[3]
        assert "MESA runs" in captured_out[3]
        assert "compress:" in captured_out[4]
        assert "/MESA_data_index3/" in captured_out[4]
        assert "Compressed MESA tracks" in captured_out[-3]
        assert "Original size" in captured_out[-2]
        assert "| Compressed size" in captured_out[-2]
        assert "" == captured_out[-1]
        # examples: reuse last with core dump files to remove (no debug output)
        for i in [0, -1]:
            with open(os.path.join(MESA_dir, os.listdir(MESA_dir)[i],\
                                   "core.3"), "w") as core_dump_file:
                core_dump_file.write("Test\nremove\n")
        test_args = totest.argparse.Namespace(mesa_dir=MESA_dir,\
                                              verbose=True, debug=False)
        totest.compress_dir(test_args)
        captured_out = capsys.readouterr().out.split('\n')
        assert "remove" in captured_out[0]
        assert "core dump files in" in captured_out[0]
        assert "directories of" in captured_out[0]
        assert "MESA runs" in captured_out[0]
        assert "compress" in captured_out[1]
        assert "files in" in captured_out[1]
        assert "directories of" in captured_out[1]
        assert "MESA runs" in captured_out[1]
        assert "Compressed MESA tracks" in captured_out[-3]
        assert "Original size" in captured_out[-2]
        assert "| Compressed size" in captured_out[-2]
        assert "" == captured_out[-1]

    def test_compress_MESA(self, monkeypatch):
        def mock_parse_commandline():
            return self.commandline_args
        def mock_set_up_test(args):
            self.test_set_up = args
        def mock_compress_dir(args):
            self.compress_dir = args
        with monkeypatch.context() as mp:
            mp.setattr(totest, "_parse_commandline", mock_parse_commandline)
            # examples: run compress_dir
            mp.setattr(totest, "compress_dir", mock_compress_dir)
            for mesa_dir in ['.', 'unit']:
                for verbose in [False, True]:
                    for debug in [False, True]:
                        self.compress_dir = None
                        self.commandline_args = totest.argparse.Namespace(\
                         mesa_dir=mesa_dir, test_dir=None, dsr=0.01,\
                         verbose=verbose, debug=debug)
                        totest._compress_MESA()
                        assert self.compress_dir.mesa_dir == mesa_dir
                        assert self.compress_dir.verbose == verbose
                        assert self.compress_dir.debug == debug
            # examples: run set_up_test
            mp.setattr(totest, "set_up_test", mock_set_up_test)
            for mesa_dir in ['.', 'unit']:
                for test_dir in ['test', 'dir']:
                    for dsr in [0.01, 0.1]:
                        for verbose in [False, True]:
                            for debug in [False, True]:
                                self.test_set_up = None
                                self.commandline_args =\
                                 totest.argparse.Namespace(mesa_dir=mesa_dir,\
                                                           test_dir=test_dir,\
                                                           dsr=dsr,\
                                                           verbose=verbose,\
                                                           debug=debug)
                                totest._compress_MESA()
                                assert self.test_set_up.mesa_dir == mesa_dir
                                assert self.test_set_up.test_dir == test_dir
                                assert self.test_set_up.dsr == dsr
                                assert self.test_set_up.verbose == verbose
                                assert self.test_set_up.debug == debug
