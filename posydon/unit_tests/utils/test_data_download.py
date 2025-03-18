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
from pytest import fixture, raises, approx
from inspect import isclass, isroutine
from shutil import rmtree
from contextlib import chdir
from unittest.mock import patch

# define test classes collecting several test functions
class TestElements:
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = ['PATH_TO_POSYDON_DATA', 'ProgressBar', '__authors__',\
                    '__builtins__', '__cached__', '__doc__', '__file__',\
                    '__loader__', '__name__', '__package__', '__spec__',\
                    'data_download', 'data_url', 'file', 'hashlib',\
                    'original_md5', 'os', 'progressbar', 'tarfile', 'tqdm',\
                    'urllib']
        assert dir(totest) == elements, "There might be added or removed "\
                                        + "objects without an update on the "\
                                        + "unit test."

    def test_instance_PATH_TO_POSYDON_DATA(self):
        assert isinstance(totest.PATH_TO_POSYDON_DATA, (str, bytes,\
                                                        os.PathLike))

    def test_instance_file(self):
        assert isinstance(totest.file, (str, bytes, os.PathLike))

    def test_instance_data_url(self):
        assert isinstance(totest.data_url, (str, bytes, os.PathLike))

    def test_instance_original_md5(self):
        assert isinstance(totest.original_md5, (str, bytes, os.PathLike))

    def test_instance_ProgressBar(self):
        assert isclass(totest.ProgressBar)

    def test_instance_data_download(self):
        assert isroutine(totest.data_download)


class TestValues:
    # check that the values fit
    def test_value_PATH_TO_POSYDON_DATA(self):
        assert "POSYDON_data" in totest.PATH_TO_POSYDON_DATA

    def test_value_file(self):
        assert "POSYDON_data.tar.gz" in totest.file

    def test_value_data_url(self):
        assert totest.data_url == "https://zenodo.org/record/14205146/files/"\
                                  + "POSYDON_data.tar.gz"

    def test_value_original_md5(self):
        assert totest.original_md5 == "cf645a45b9b92c2ad01e759eb1950beb"


class TestFunctions:
    @fixture
    def test_path(self, tmp_path):
        # a temporary path to POSYDON_data for testing
        return os.path.join(tmp_path, "POSYDON_data")

    @fixture
    def download_statement(self):
        # statement that the download started
        return "Downloading POSYDON data from Zenodo to PATH_TO_POSYDON_DATA="

    @fixture
    def failed_MD5_statement(self):
        # statement that MD5 verfication failed
        return "Failed to read the tar.gz file for MD5 verification"

    @fixture
    def extraction_statement(self):
        # statement that the tar extraction started
        return "Extracting POSYDON data from tar file..."

    @fixture
    def removal_statement(self):
        # statement that the tar file gets removed
        return "Removed downloaded tar file."

    # test functions
    def test_data_download(self, capsys, monkeypatch, test_path,\
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
        with raises(TypeError, match="path should be string, bytes, "\
                                     +"os.PathLike or integer"):
            totest.data_download(file={})
        # bad input
        with monkeypatch.context() as mp:
            mock_environ = totest.os.environ.copy()
            mock_environ.pop("PATH_TO_POSYDON_DATA")
            mp.setattr(totest.os, "environ", mock_environ)
            with raises(NameError, match="You must define the "\
                                         +"PATH_TO_POSYDON_DATA environment "\
                                         +"variable before downloading "\
                                         +"POSYDON datasets"):
                totest.data_download(file=test_path+".tar.gz")
        # bad input
        with raises(FileExistsError, match="POSYDON data already exists at "\
                                           +"./"):
            totest.data_download(file="./")
        # skip real download: do nothing instead
        with monkeypatch.context() as mp:
            mp.setattr(totest.urllib.request, "urlretrieve", mock_urlretrieve)
            mp.setattr(totest.tarfile, "open", mock_open)
            # mocked download: fails MD5check, no tar file to remove
            totest.data_download(file=test_path+".tar.gz", verbose=True)
            captured_output = capsys.readouterr()
            assert download_statement in captured_output.out
            assert failed_MD5_statement in captured_output.out
            assert extraction_statement in captured_output.out
            assert removal_statement not in captured_output.out
            # without MD5 check
            totest.data_download(file=test_path+".tar.gz", MD5_check=False)
            captured_output = capsys.readouterr()
            assert download_statement in captured_output.out
            assert failed_MD5_statement not in captured_output.out
            assert extraction_statement in captured_output.out
            assert removal_statement not in captured_output.out
        # skip real download: create mock file instead
        with monkeypatch.context() as mp:
            mp.setattr(totest.urllib.request, "urlretrieve", mock_urlretrieve2)
            # mocked download: fails MD5check, which removes tar file and
            #                  causes a FileNotFoundError
            with raises(FileNotFoundError, match="No such file or directory:"):
                totest.data_download(file=test_path+"0.tar.gz", verbose=True)
            captured_output = capsys.readouterr()
            assert download_statement in captured_output.out
            assert failed_MD5_statement in captured_output.out
            assert extraction_statement in captured_output.out
            assert removal_statement not in captured_output.out
            # return expected hash for testfile
            with patch("hashlib.md5") as p_md5:
                p_md5.return_value.hexdigest.return_value = totest.original_md5
                totest.data_download(file=test_path+"1.tar.gz")
                captured_output = capsys.readouterr()
                assert download_statement in captured_output.out
                assert failed_MD5_statement not in captured_output.out
                assert extraction_statement in captured_output.out
                assert removal_statement not in captured_output.out
                assert os.path.exists(test_path)
                assert os.path.exists(os.path.join(test_path, "test1.txt"))
                rmtree(test_path) # removed extracted data for next test
                # with verification output
                totest.data_download(file=test_path+"2.tar.gz", verbose=True)
                captured_output = capsys.readouterr()
                assert download_statement in captured_output.out
                assert "MD5 verified" in captured_output.out
                assert extraction_statement in captured_output.out
                assert removal_statement in captured_output.out
                assert os.path.exists(test_path)
                assert os.path.exists(os.path.join(test_path, "test2.txt"))
#                rmtree(test_path) # removed extracted data for next test


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
