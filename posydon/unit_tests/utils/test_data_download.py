"""Unit tests of posydon/utils/data_download.py

"""

__authors__ = [
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>"
]

# import the module which will be tested
import posydon.utils.data_download as totest

# aliases
os = totest.os

# import other needed code for the tests, which is not already imported in the
# module you like to test
import hashlib
from contextlib import chdir
from inspect import isclass, isroutine
from shutil import copyfile, rmtree
from unittest.mock import patch

from pytest import approx, fixture, raises, warns

from posydon.utils.posydonwarning import ReplaceValueWarning


# define test classes collecting several test functions
class TestElements:
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = {'COMPLETE_SETS', 'PATH_TO_POSYDON_DATA', 'ProgressBar',\
                    'Pwarn', 'ZENODO_COLLECTION', '__authors__',\
                    '__builtins__', '__cached__', '__doc__', '__file__',\
                    '__loader__', '__name__', '__package__', '__spec__',\
                    '_GRID_DIRS', '_INTERP_GRID_DIRS', '_INTERP_METHODS',\
                    '_dataset_installed', '_expected_paths', '_md5_of_file',\
                    '_get_posydon_data', '_parse_commandline', 'argparse',\
                    'convert_metallicity_to_string', 'data_download',\
                    'download_one_dataset', 'hashlib', 'list_datasets', 'os',\
                    'progressbar', 'tarfile', 'textwrap', 'tqdm', 'urllib'}
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

    def test_instance_ProgressBar(self):
        assert isclass(totest.ProgressBar)

    def test_instance_list_datasets(self):
        assert isroutine(totest.list_datasets)

    def test_instance_expected_paths(self):
        assert isroutine(totest._expected_paths)

    def test_instance_dataset_installed(self):
        assert isroutine(totest._dataset_installed)

    def test_instance_md5_of_file(self):
        assert isroutine(totest._md5_of_file)

    def test_instance_download_one_dataset(self):
        assert isroutine(totest.download_one_dataset)

    def test_instance_data_download(self):
        assert isroutine(totest.data_download)

    def test_instance_get_posydon_data(self):
        assert isroutine(totest._get_posydon_data)


class TestFunctions:
    @fixture
    def test_path(self, tmp_path):
        # a temporary path to POSYDON_data for testing
        return os.path.join(tmp_path, "POSYDON_data")

    @fixture
    def download_statement(self):
        # statement that the download started
        return "Downloading POSYDON data '{}' from Zenodo to "

    @fixture
    def failed_MD5_statement(self):
        # statement that MD5 verfication failed
        return "Failed to read the tar.gz file for MD5 verification"

    @fixture
    def extraction_statement(self):
        # statement that the tar extraction started
        return "Extracting POSYDON data '{}' from tar file..."

    @fixture
    def removal_statement(self):
        # statement that the tar file gets removed
        return "Removed downloaded tar file."

    # test functions
    def test_parse_commandline(self, monkeypatch):
        def mock_parse_args(parser):
            return self.commandline_args
        def mock_error(parser, message):
            raise totest.argparse.ArgumentError(None, message)
        with monkeypatch.context() as mp:
            mp.setattr(totest.argparse.ArgumentParser, "parse_args",
                       mock_parse_args)
            mp.setattr(totest.argparse.ArgumentParser, "error", mock_error)
            # bad input
            self.commandline_args = totest.argparse.Namespace(dataset='Test',\
             listedsets=None, nomd5check=False, force=False, verbose=False)
            with raises(totest.argparse.ArgumentError,\
                        match="unknown dataset, use -l to show defined sets"):
                totest._parse_commandline()
            # example
            self.commandline_args = totest.argparse.Namespace(dataset='DR1',\
             listedsets=None, nomd5check=False, force=False, verbose=False)
            assert totest._parse_commandline() == self.commandline_args

    def test_list_datasets(self, capsys):
        for v in [True, False]:
            # individual sets
            totest.list_datasets(individual_sets=True, verbose=v)
            captured_output = capsys.readouterr()
            assert "Defined individual sets are:" in captured_output.out
            if v:
                assert "more information at" in captured_output.out
            # complete sets
            totest.list_datasets(individual_sets=False, verbose=v)
            captured_output = capsys.readouterr()
            assert "Defined complete sets are:" in captured_output.out
            if v:
                assert "more information at" in captured_output.out

    def test_expected_paths(self):
        # grid data sets contain the grids and their pre-trained interpolators
        for suffix, z_str in [('2', '2e+00'), ('0.45', '4.5e-01'),
                              ('1e-3', '1e-03')]:
            expected_paths = totest._expected_paths(
                'DR2_grids_'+suffix+'Zsun')
            assert len(expected_paths)\
                   == (len(totest._GRID_DIRS)
                       +(len(totest._INTERP_GRID_DIRS)
                         *len(totest._INTERP_METHODS)))
            for grid_dir in totest._GRID_DIRS:
                assert os.path.join(grid_dir, z_str+"_Zsun.h5")\
                       in expected_paths
            for grid_dir in totest._INTERP_GRID_DIRS:
                for interp_method in totest._INTERP_METHODS:
                    assert os.path.join(grid_dir, "interpolators",
                                        interp_method,
                                        z_str+"_Zsun.pkl")\
                           in expected_paths
        # unsupported metallicities cannot be verified
        assert totest._expected_paths('DR2_grids_bogusZsun') is None
        assert totest._expected_paths('DR2_grids_5Zsun') is None
        assert totest._expected_paths('DR2_grids_Zsun') is None
        # other known data sets
        assert totest._expected_paths('auxiliary') is not None
        assert totest._expected_paths('DR1_for_v2.0.0-pre1') is not None
        assert totest._expected_paths('DR1-super_Eddington') is not None
        assert len(totest._expected_paths('DR1-super_Eddington')) == 3
        # unknown data sets cannot be verified
        assert totest._expected_paths('DR2') is None
        assert totest._expected_paths('v2_tutorial_populations') is None

    def test_dataset_installed(self, monkeypatch, test_path):
        with monkeypatch.context() as mp:
            mp.setattr(totest, "PATH_TO_POSYDON_DATA", test_path)
            # data sets with unknown content cannot be verified
            assert totest._dataset_installed('Test') == False
            assert totest._dataset_installed('v2_tutorial_populations')\
                   == False
            # partially extracted data sets are not installed
            auxiliary_paths = totest._expected_paths('auxiliary')
            for path in auxiliary_paths[:-1]:
                os.makedirs(os.path.join(test_path, os.path.dirname(path)),
                            exist_ok=True)
                with open(os.path.join(test_path, path), "a"):
                    pass
            assert totest._dataset_installed('auxiliary') == False
            # fully extracted data sets are installed
            with open(os.path.join(test_path, auxiliary_paths[-1]), "a"):
                pass
            assert totest._dataset_installed('auxiliary') == True

    def test_md5_of_file(self, tmp_path):
        content = b"Unit Test\n"
        test_file_path = os.path.join(tmp_path, "unit_test.txt")
        with open(test_file_path, "wb") as test_file:
            test_file.write(content)
        assert totest._md5_of_file(test_file_path)\
               == hashlib.md5(content).hexdigest()

    def test_download_one_dataset(self, capsys, monkeypatch, test_path,\
                                  download_statement, failed_MD5_statement,\
                                  extraction_statement, removal_statement):
        # mocks
        def failing_urlretrieve(url, filename=None, reporthook=None,\
                                data=None):
            # raise an error, if a download should not happen
            raise RuntimeError("The dataset should not be downloaded.")
        def mock_urlretrieve_empty(url, filename=None, reporthook=None,\
                                   data=None):
            # simulate a download by creating an empty file
            if isinstance(filename, str):
                with open(filename, "wb"):
                    pass
            return None
        class mock_TarFile:
            def getmembers(self):
                return []
            def extract(self, member, path='', set_attrs=True, *,\
                        numeric_owner=False, filter=None):
                return
        class mock_open:
            def __init__(self, name=None, mode='r', fileobj=None,\
                         bufsize=10240, **kwargs):
                pass
            def __enter__(self):
                return mock_TarFile()
            def __exit__(self, exc_type, exc_value, exc_traceback):
                return False

        # helpers
        directory = os.path.dirname(test_path)
        filepath = os.path.join(directory, "POSYDON_data.tar.gz")
        partpath = filepath+".part"
        pristine_archive = os.path.join(directory,
                                        "POSYDON_data_pristine.tar.gz")

        def mock_urlretrieve_archive(url, filename=None, reporthook=None,\
                                     data=None):
            # simulate a download by providing the correct archive
            copyfile(pristine_archive, filename)
            return None

        def build_archive():
            # create a small tar.gz archive at filepath
            os.makedirs(test_path, exist_ok=True)
            with open(os.path.join(test_path, "test.txt"), "w")\
                 as test_file:
                test_file.write("Unit Test\n")
            with chdir(directory):
                if os.system("tar -czf POSYDON_data.tar.gz POSYDON_data")\
                   != 0:
                    raise RuntimeError("Please check that you have `tar` "\
                                       +"installed and up to date.")
            rmtree(test_path)

        def clean_up():
            # remove leftovers of previous test scenarios
            for leftover in [filepath, partpath]:
                if os.path.exists(leftover):
                    os.remove(leftover)
            if os.path.exists(test_path):
                rmtree(test_path)

        # bad input
        with raises(TypeError, match="'dataset' should be a string."):
            totest.download_one_dataset(dataset=None)
        # bad input
        with raises(KeyError, match="The dataset 'Test' is not defined."):
            totest.download_one_dataset(dataset='Test')
        # bad input
        with monkeypatch.context() as mp:
            mock_ZENODO_COLLECTION = {'Test': {'data' : None, 'md5': None}}
            mp.setattr(totest, "ZENODO_COLLECTION", mock_ZENODO_COLLECTION)
            with raises(ValueError, match="The dataset 'Test' has no "\
                                          +"publication yet."):
                totest.download_one_dataset(dataset='Test')
        # bad input
        with monkeypatch.context() as mp:
            mock_PATH_TO_POSYDON_DATA = "./DOES_NOT_EXIST/POSYDON_data"
            mp.setattr(totest, "PATH_TO_POSYDON_DATA",\
                       mock_PATH_TO_POSYDON_DATA)
            with raises(NotADirectoryError, match="PATH_TO_POSYDON_DATA does "\
                                                  +"not refer to a valid "\
                                                  +"directory."):
                totest.download_one_dataset(dataset='DR1_for_v2.0.0-pre1')

        # installed data sets get skipped
        with monkeypatch.context() as mp:
            mp.setattr(totest, "PATH_TO_POSYDON_DATA", test_path)
            mp.setattr(totest, "ZENODO_COLLECTION",
                       {'auxiliary': totest.ZENODO_COLLECTION['auxiliary']})
            mp.setattr(totest.urllib.request, "urlretrieve",
                       failing_urlretrieve)
            for path in totest._expected_paths('auxiliary'):
                os.makedirs(os.path.join(test_path, os.path.dirname(path)),
                            exist_ok=True)
                with open(os.path.join(test_path, path), "a"):
                    pass
            totest.download_one_dataset(dataset='auxiliary')
            captured_output = capsys.readouterr()
            assert "POSYDON data 'auxiliary' is already present, skipping."\
                   in captured_output.out
            clean_up()
            # forced downloads happen anyway
            mp.setattr(totest.urllib.request, "urlretrieve",
                       mock_urlretrieve_empty)
            mp.setattr(totest.tarfile, "open", mock_open)
            totest.download_one_dataset(dataset='auxiliary', MD5_check=False,
                                        force=True)
            captured_output = capsys.readouterr()
            assert download_statement.format('auxiliary')\
                   in captured_output.out
            assert extraction_statement.format('auxiliary')\
                   in captured_output.out
            clean_up()

        # downloads without an available MD5 checksum cannot be verified
        with monkeypatch.context() as mp:
            mp.setattr(totest, "PATH_TO_POSYDON_DATA", test_path)
            mock_ZENODO_COLLECTION = {'Test': {'data': "POSYDON_data.tar.gz",
                                               'md5': None}}
            mp.setattr(totest, "ZENODO_COLLECTION", mock_ZENODO_COLLECTION)
            mp.setattr(totest.urllib.request, "urlretrieve",
                       mock_urlretrieve_empty)
            mp.setattr(totest.tarfile, "open", mock_open)
            with warns(ReplaceValueWarning, match="MD5 undefined, skip MD5 "\
                                                  +"check."):
                totest.download_one_dataset(dataset='Test')
            captured_output = capsys.readouterr()
            assert download_statement.format('Test') in captured_output.out
            assert failed_MD5_statement not in captured_output.out
            assert extraction_statement.format('Test') in captured_output.out
            assert removal_statement not in captured_output.out
            clean_up()
            # incomplete leftover downloads get removed and restarted
            with open(partpath, "wb"):
                pass
            totest.download_one_dataset(dataset='Test', verbose=True)
            captured_output = capsys.readouterr()
            assert "Removing incomplete download 'POSYDON_data.tar.gz"\
                   +".part'" in captured_output.out
            assert download_statement.format('Test') in captured_output.out
            assert extraction_statement.format('Test') in captured_output.out
            clean_up()

        # corrupted fresh downloads get removed and raise an error
        with monkeypatch.context() as mp:
            mp.setattr(totest, "PATH_TO_POSYDON_DATA", test_path)
            mock_ZENODO_COLLECTION = {'Test': {'data': "POSYDON_data.tar.gz",
                                               'md5': "Unit"}}
            mp.setattr(totest, "ZENODO_COLLECTION", mock_ZENODO_COLLECTION)
            mp.setattr(totest.urllib.request, "urlretrieve",
                       mock_urlretrieve_empty)
            mp.setattr(totest.tarfile, "open", mock_open)
            with raises(ValueError, match="MD5 verification failed!."):
                totest.download_one_dataset(dataset='Test')
            captured_output = capsys.readouterr()
            assert download_statement.format('Test') in captured_output.out
            assert extraction_statement.format('Test') not in\
                captured_output.out
            assert failed_MD5_statement not in captured_output.out
            assert not os.path.exists(filepath)
            clean_up()
            # unreadable files cannot be verified, but continue anyway;
            # here, the archive vanishes before its removal
            def failing_md5(unused_filepath):
                os.remove(filepath)
                raise OSError("unreadable file")
            mp.setattr(totest, "_md5_of_file", failing_md5)
            totest.download_one_dataset(dataset='Test')
            captured_output = capsys.readouterr()
            assert failed_MD5_statement in captured_output.out
            assert extraction_statement.format('Test') in captured_output.out
            assert removal_statement not in captured_output.out
            clean_up()

        # existing complete archives are verified and extracted instead of
        # being downloaded again
        build_archive()
        copyfile(filepath, pristine_archive)
        original_md5 = totest._md5_of_file(filepath)
        with monkeypatch.context() as mp:
            mp.setattr(totest, "PATH_TO_POSYDON_DATA", test_path)
            mock_ZENODO_COLLECTION = {'Test': {'data': "POSYDON_data.tar.gz",
                                               'md5': original_md5}}
            mp.setattr(totest, "ZENODO_COLLECTION", mock_ZENODO_COLLECTION)
            mp.setattr(totest.urllib.request, "urlretrieve",
                       failing_urlretrieve)
            # correct checksum: reuse the existing archive
            totest.download_one_dataset(dataset='Test', verbose=True)
            captured_output = capsys.readouterr()
            assert "Verifying existing archive 'POSYDON_data.tar.gz'"\
                   in captured_output.out
            assert "MD5 verified." in captured_output.out
            assert download_statement.format('Test') not in\
                captured_output.out
            assert extraction_statement.format('Test') in captured_output.out
            assert removal_statement in captured_output.out
            assert os.path.exists(os.path.join(test_path, "test.txt"))
            clean_up()
            # without verbose output
            build_archive()
            totest.download_one_dataset(dataset='Test')
            captured_output = capsys.readouterr()
            assert "Verifying existing archive" not in captured_output.out
            clean_up()
            # wrong checksum: replace the existing archive by a fresh download
            build_archive()
            with open(filepath, "ab") as corrupted_archive:
                corrupted_archive.write(b"corrupted\n")
            mp.setattr(totest.urllib.request, "urlretrieve",
                       mock_urlretrieve_archive)
            totest.download_one_dataset(dataset='Test', verbose=True)
            captured_output = capsys.readouterr()
            assert "The existing archive did not pass the MD5 verification,"\
                   in captured_output.out
            assert download_statement.format('Test') in captured_output.out
            assert "MD5 verified." in captured_output.out
            assert extraction_statement.format('Test') in captured_output.out
            assert removal_statement in captured_output.out
            assert os.path.exists(os.path.join(test_path, "test.txt"))
            clean_up()
        os.remove(pristine_archive)

    def test_data_download(self, capsys, monkeypatch):
        def mock_download_one_dataset(**kwargs):
            self.kwargs = kwargs
            return
        # bad input
        with raises(TypeError, match="'set_name' should be a string."):
            totest.data_download(set_name=None)
        # bad input
        with raises(KeyError, match="The dataset 'Test' is not defined."):
            totest.data_download(set_name='Test')
        # skip real download: do nothing instead
        with monkeypatch.context() as mp:
            # test single dataset
            mock_ZENODO_COLLECTION = {'Test': {}, 'Test2': {}}
            mp.setattr(totest, "ZENODO_COLLECTION", mock_ZENODO_COLLECTION)
            mp.setattr(totest, "download_one_dataset", mock_download_one_dataset)
            for v in [True, False]:
                self.kwargs = None
                totest.data_download(set_name='Test', verbose=v)
                assert self.kwargs['dataset'] == 'Test'
                assert self.kwargs['MD5_check'] == True
                assert self.kwargs['verbose'] == v
                assert self.kwargs['force'] == False
                if v:
                    assert "You are downloading a single data set, which "\
                           + "might not contain all the data needed.\n"\
                           == capsys.readouterr().out
            # test forced download of a single dataset
            self.kwargs = None
            totest.data_download(set_name='Test', force=True)
            assert self.kwargs['dataset'] == 'Test'
            assert self.kwargs['MD5_check'] == True
            assert self.kwargs['verbose'] == False
            assert self.kwargs['force'] == True
            # test complete set
            mock_COMPLETE_SETS = {'Unit': ['Test', 'Test2']}
            mp.setattr(totest, "COMPLETE_SETS", mock_COMPLETE_SETS)
            self.kwargs = None
            totest.data_download(set_name='Unit')
            assert self.kwargs['dataset'] == 'Test2'
            assert self.kwargs['MD5_check'] == True
            assert self.kwargs['verbose'] == False
            assert self.kwargs['force'] == False

    def test_get_posydon_data(self, monkeypatch):
        def mock_parse_commandline():
            return self.commandline_args
        def mock_list_datasets(**kwargs):
            self.list_printed = kwargs
        def mock_data_download(**kwargs):
            self.downloaded = kwargs
        with monkeypatch.context() as mp:
            mp.setattr(totest, "_parse_commandline", mock_parse_commandline)
            mp.setattr(totest, "list_datasets", mock_list_datasets)
            mp.setattr(totest, "data_download", mock_data_download)
            # examples
            for l in ['complete', 'individual']:
                for v in [True, False]:
                    self.list_printed = None
                    self.commandline_args = totest.argparse.Namespace(\
                     dataset='DR1', listedsets=l, nomd5check=False,
                     force=False, verbose=v)
                    totest._get_posydon_data()
                    assert self.list_printed['individual_sets']\
                           == (l == 'individual')
                    assert self.list_printed['verbose'] == v
            # examples
            for n in [True, False]:
                for f in [True, False]:
                    for v in [True, False]:
                        for d in ['DR1', 'DR2']:
                            self.downloaded = None
                            self.commandline_args = totest.argparse.Namespace(\
                             dataset=d, listedsets=None, nomd5check=n,
                             force=f, verbose=v)
                            totest._get_posydon_data()
                            assert self.downloaded['set_name'] == d
                            assert self.downloaded['MD5_check'] == (not n)
                            assert self.downloaded['verbose'] == v
                            assert self.downloaded['force'] == f


class TestProgressBar:
    @fixture
    def ProgressBar(self):
        # initialize an instance of the class with defaults
        return totest.ProgressBar()

    # test the ProgressBar class
    def test_init(self, ProgressBar):
        assert isroutine(ProgressBar.__init__)
        # check that the instance is of correct type and all code in the
        # __init__ got executed: the elements are created and initialized
        assert isinstance(ProgressBar, totest.ProgressBar)
        assert ProgressBar.pbar is None
        assert isinstance(ProgressBar.widgets, list)

    def test_call(self, ProgressBar):
        assert isroutine(ProgressBar.__call__)
        # missing argument
        with raises(TypeError, match="missing 3 required positional "\
                                     +"arguments: 'block_num', 'block_size', "\
                                     +"and 'total_size'"):
            ProgressBar()
        # bad input
        with raises(TypeError, match="'<' not supported between instances of "\
                                     +"'str' and 'int'"):
            ProgressBar("Test", 1, 1)
            # the progressbar starts before the error
        # hence, tearDown for pbar needed
        ProgressBar.pbar = None
        with raises(TypeError, match="'<' not supported between instances of "\
                                     +"'str' and 'int'"):
            ProgressBar(1, "Test", 1)
            # the progressbar starts before the error
        # hence, tearDown for pbar needed
        ProgressBar.pbar = None
        with raises(TypeError, match="'>' not supported between instances of "\
                                     +"'int' and 'str'"):
            ProgressBar(1, 1, "Test")
            # the progressbar starts before the error
        # hence, tearDown for pbar needed
        ProgressBar.pbar = None
        for i in range(9):
            ProgressBar(i, 1, 8)
            assert ProgressBar.pbar.percentage == approx(i*12.5)
        ProgressBar(9, 1, 8)
        assert ProgressBar.pbar.percentage == approx(100.0)
