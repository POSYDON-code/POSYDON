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
from pytest import fixture, raises, warns, approx
from inspect import isclass, isroutine
from shutil import rmtree
from contextlib import chdir
from unittest.mock import patch
from posydon.utils.posydonwarning import ReplaceValueWarning

# define test classes collecting several test functions
class TestElements:
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = {'COMPLETE_SETS', 'PATH_TO_POSYDON_DATA', 'ProgressBar',\
                    'Pwarn', 'ZENODO_COLLECTION', '__authors__',\
                    '__builtins__', '__cached__', '__doc__', '__file__',\
                    '__loader__', '__name__', '__package__', '__spec__',\
                    '_get_posydon_data', '_parse_commandline', 'argparse',\
                    'data_download', 'download_one_dataset', 'hashlib',\
                    'list_datasets', 'os', 'progressbar', 'tarfile',\
                    'textwrap', 'tqdm', 'urllib'}
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
             listedsets=None, nomd5check=False, verbose=False)
            with raises(totest.argparse.ArgumentError,\
                        match="unknown dataset, use -l to show defined sets"):
                totest._parse_commandline()
            # example
            self.commandline_args = totest.argparse.Namespace(dataset='v1',\
             listedsets=None, nomd5check=False, verbose=False)
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

    def test_download_one_dataset(self, capsys, monkeypatch, test_path,\
                                  download_statement, failed_MD5_statement,\
                                  extraction_statement, removal_statement):
        def mock_urlretrieve(url, filename=None, reporthook=None, data=None):
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
        def mock_urlretrieve2(url, filename=None, reporthook=None, data=None):
            os.mkdir(test_path)
            if (isinstance(filename, str) and (len(filename)>=8)):
                test_ID = filename[-8]
            else:
                test_ID = ""
            with open(os.path.join(test_path, "test"+test_ID+".txt"), "w")\
                 as test_file:
                test_file.write("Unit Test\n")
            with chdir(os.path.join(test_path,"..")):
                os.system("tar -czf POSYDON_data"+test_ID\
                          +".tar.gz POSYDON_data")
            rmtree(test_path)
            return None

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
                totest.download_one_dataset(dataset='v1_for_v2.0.0-pre1')
        # bad input
        with monkeypatch.context() as mp:
            mock_ZENODO_COLLECTION = {'Test': {'data' : "./", 'md5': "Unit"}}
            mp.setattr(totest, "ZENODO_COLLECTION", mock_ZENODO_COLLECTION)
            with raises(FileExistsError, match="POSYDON data already exists "\
                                               +"at"):
                totest.download_one_dataset(dataset='Test')
        # skip real download: do nothing instead
        with monkeypatch.context() as mp:
            mp.setattr(totest, "PATH_TO_POSYDON_DATA", test_path)
            mock_ZENODO_COLLECTION = {'Test': {'data': "POSYDON_data.tar.gz",
                                               'md5': "Unit"}}
            mp.setattr(totest, "ZENODO_COLLECTION", mock_ZENODO_COLLECTION)
            mp.setattr(totest.urllib.request, "urlretrieve", mock_urlretrieve)
            mp.setattr(totest.tarfile, "open", mock_open)
            # mocked download: fails MD5check, no tar file to remove
            totest.download_one_dataset(dataset='Test', verbose=True)
            captured_output = capsys.readouterr()
            assert download_statement.format('Test') in captured_output.out
            assert failed_MD5_statement in captured_output.out
            assert extraction_statement.format('Test') in captured_output.out
            assert removal_statement not in captured_output.out
            # without MD5 check
            totest.download_one_dataset(dataset='Test', MD5_check=False)
            captured_output = capsys.readouterr()
            assert download_statement.format('Test') in captured_output.out
            assert failed_MD5_statement not in captured_output.out
            assert extraction_statement.format('Test') in captured_output.out
            assert removal_statement not in captured_output.out
            # without MD5 check
            mock_ZENODO_COLLECTION['Test']['md5'] = None
            with warns(ReplaceValueWarning, match="MD5 undefined, skip MD5 "\
                                                  +"check."):
                totest.download_one_dataset(dataset='Test')
            captured_output = capsys.readouterr()
            assert download_statement.format('Test') in captured_output.out
            assert failed_MD5_statement not in captured_output.out
            assert extraction_statement.format('Test') in captured_output.out
            assert removal_statement not in captured_output.out
        # skip real download: create mock file instead
        with monkeypatch.context() as mp:
            mp.setattr(totest, "PATH_TO_POSYDON_DATA", test_path)
            mock_ZENODO_COLLECTION = {'Test0': {'data': "POSYDON_data0.tar.gz",
                                                'md5': "Unit"},
                                      'Test1': {'data': "POSYDON_data1.tar.gz",
                                                'md5': "Unit"},
                                      'Test2': {'data': "POSYDON_data2.tar.gz",
                                                'md5': "Unit"}}
            mp.setattr(totest, "ZENODO_COLLECTION", mock_ZENODO_COLLECTION)
            mp.setattr(totest.urllib.request, "urlretrieve", mock_urlretrieve2)
            # mocked download: fails MD5check, which removes tar file and
            #                  causes a FileNotFoundError
            with raises(FileNotFoundError, match="No such file or directory:"):
                totest.download_one_dataset(dataset='Test0', verbose=True)
            captured_output = capsys.readouterr()
            assert download_statement.format('Test0') in captured_output.out
            assert failed_MD5_statement in captured_output.out
            assert extraction_statement.format('Test0') in captured_output.out
            assert removal_statement not in captured_output.out
            # return expected hash for testfile
            with patch("hashlib.md5") as p_md5:
                p_md5.return_value.hexdigest.return_value\
                 = mock_ZENODO_COLLECTION['Test1']['md5']
                totest.download_one_dataset(dataset='Test1')
                captured_output = capsys.readouterr()
                assert download_statement.format('Test1')\
                       in captured_output.out
                assert failed_MD5_statement not in captured_output.out
                assert extraction_statement.format('Test1')\
                       in captured_output.out
                assert removal_statement not in captured_output.out
                assert os.path.exists(test_path)
                assert os.path.exists(os.path.join(test_path, "test1.txt"))
                rmtree(test_path) # removed extracted data for next test
                # with verification output
                totest.download_one_dataset(dataset='Test2', verbose=True)
                captured_output = capsys.readouterr()
                assert download_statement.format('Test2')\
                       in captured_output.out
                assert "MD5 verified" in captured_output.out
                assert extraction_statement.format('Test2')\
                       in captured_output.out
                assert removal_statement in captured_output.out
                assert os.path.exists(test_path)
                assert os.path.exists(os.path.join(test_path, "test2.txt"))
#                rmtree(test_path) # removed extracted data for next test

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
                if v:
                    assert "You are downloading a single data set, which "\
                           + "might not contain all the data needed.\n"\
                           == capsys.readouterr().out
            # test complete set
            mock_COMPLETE_SETS = {'Unit': ['Test', 'Test2']}
            mp.setattr(totest, "COMPLETE_SETS", mock_COMPLETE_SETS)
            self.kwargs = None
            totest.data_download(set_name='Unit')
            assert self.kwargs['dataset'] == 'Test2'
            assert self.kwargs['MD5_check'] == True
            assert self.kwargs['verbose'] == False

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
                     dataset='v1', listedsets=l, nomd5check=False, verbose=v)
                    totest._get_posydon_data()
                    assert self.list_printed['individual_sets']\
                           == (l == 'individual')
                    assert self.list_printed['verbose'] == v
            # examples
            for n in [True, False]:
                for v in [True, False]:
                    for d in ['v1', 'v2']:
                        self.downloaded = None
                        self.commandline_args = totest.argparse.Namespace(\
                         dataset=d, listedsets=None, nomd5check=n, verbose=v)
                        totest._get_posydon_data()
                        assert self.downloaded['set_name'] == d
                        assert self.downloaded['MD5_check'] == (not n)
                        assert self.downloaded['verbose'] == v


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
