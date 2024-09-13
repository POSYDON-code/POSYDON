"""Unit tests of posydon/utils/data_download.py
"""

__authors__ = [
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>"
]

# import the unittest module and the module which will be tested
import unittest
import posydon.utils.data_download as totest

# import other needed code for the tests
from contextlib import redirect_stdout
from inspect import isclass, isroutine
from io import StringIO
from shutil import rmtree
import sys

DO_DOWNLOAD = False

# define test classes
class TestElements(unittest.TestCase):
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = ['PATH_TO_POSYDON_DATA', 'ProgressBar', '__authors__',
                    '__builtins__', '__cached__', '__doc__', '__file__',
                    '__loader__', '__name__', '__package__', '__spec__',
                    'data_download', 'data_url', 'file', 'hashlib',
                    'original_md5', 'os', 'progressbar', 'tarfile', 'tqdm',
                    'urllib']
        self.assertListEqual(dir(totest), elements,
                             msg="There might be added or removed objects "
                                 "without an update on the unit test.")

    def test_instance_PATH_TO_POSYDON_DATA(self):
        self.assertIsInstance(totest.PATH_TO_POSYDON_DATA, (str, bytes,
                                                           totest.os.PathLike))

    def test_instance_file(self):
        self.assertIsInstance(totest.file, (str, bytes, totest.os.PathLike))

    def test_instance_data_url(self):
        self.assertIsInstance(totest.data_url, (str, bytes,
                                                totest.os.PathLike))

    def test_instance_original_md5(self):
        self.assertIsInstance(totest.original_md5, (str, bytes,
                                                    totest.os.PathLike))

    def test_instance_ProgressBar(self):
        self.assertTrue(isclass(totest.ProgressBar))

    def test_instance_data_download(self):
        self.assertTrue(isroutine(totest.data_download))


class TestValues(unittest.TestCase):
    # check that the values fit
    def test_value_PATH_TO_POSYDON_DATA(self):
        self.assertIn("POSYDON_data", totest.PATH_TO_POSYDON_DATA)

    def test_value_file(self):
        self.assertIn("POSYDON_data.tar.gz", totest.file)

    def test_value_data_url(self):
        self.assertEqual("https://zenodo.org/record/6655751/files/"
                         "POSYDON_data.tar.gz", totest.data_url)

    def test_value_original_md5(self):
        self.assertIn("8873544d9a568ebb85bccffbf1bdcd99", totest.original_md5)


class TestFunctions(unittest.TestCase):
    # test functions
    def test_data_download(self):
        # bad input
        with self.assertRaises(TypeError):
            totest.data_download(file={})
        with redirect_stdout(StringIO()) as print_out:
            totest.data_download(file="./", verbose=True)
        self.assertEqual("POSYDON data alraedy exists at ./\n",
                         print_out.getvalue())

    @unittest.skipUnless(DO_DOWNLOAD, "Usually, skip the test on the actual "
                                      "download.")
    def test_data_download(self):
        # real download: may add skip option
        tmp_path = "./tmp"
        test_path = totest.os.path.join(tmp_path, "POSYDON_data")
        if not totest.os.path.exists(test_path):
            totest.os.makedirs(test_path)
        totest.PATH_TO_POSYDON_DATA = test_path
        with redirect_stdout(StringIO()) as print_out:
            totest.data_download(file=test_path+".tar.gz", verbose=True)
        self.assertIn("Downloading POSYDON data from Zenodo to "
                      f"PATH_TO_POSYDON_DATA={test_path}",
                      print_out.getvalue())
        self.assertIn("MD5 verified.", print_out.getvalue())
        self.assertIn("Extracting POSYDON data from tar file...",
                      print_out.getvalue())
        self.assertIn("Removed downloaded tar file.", print_out.getvalue())
        with self.assertRaises(FileNotFoundError):
            totest.original_md5 += '_'
            with redirect_stdout(StringIO()) as print_out:
                totest.data_download(file=test_path+".tar.gz", verbose=True)
        if totest.original_md5[-1]=='_':
            totest.original_md5 = totest.original_md5[:-1]
        self.assertIn("Failed to read the tar.gz file for MD5 verificaton, "
                      "cannot guarantee file integrity (this error seems to "
                      "happen only on macOS).", print_out.getvalue())
        rmtree(tmp_path)


class TestProgressBar(unittest.TestCase):
    # test the ProgressBar class
    def setUp(self):
        # initialize an instance of the class for each test
        self.ProgressBar = totest.ProgressBar()

    def tearDown(self):
        def finish_remove(pbar):
            # finish progress and remove the bar
            pbar.finish(end='\r' + ' ' * pbar.term_width + '\r')

        # finish, remove, and reset pbar if needed
        if self.ProgressBar.pbar:
            finish_remove(self.ProgressBar.pbar)
            self.ProgressBar.pbar = None

    def test_init(self):
        self.assertTrue(isroutine(self.ProgressBar.__init__))
        # check that the instance is of correct type and all code in the
        # __init__ got executed: the elements are created and initialized
        self.assertIsInstance(self.ProgressBar, totest.ProgressBar)
        self.assertIsNone(self.ProgressBar.pbar)
        self.assertIsInstance(self.ProgressBar.widgets, list)

    def test_call(self):
        self.assertTrue(isroutine(self.ProgressBar.__call__))
        # missing argument
        with self.assertRaises(TypeError):
            self.ProgressBar()
        # bad input
        with self.assertRaises(TypeError):
            self.ProgressBar("Test", 1, 1)
            # the progressbar starts before the error, hence tearDown
        self.tearDown()
        with self.assertRaises(TypeError):
            self.ProgressBar(1, "Test", 1)
            # the progressbar starts before the error, hence tearDown
        self.tearDown()
        with self.assertRaises(TypeError):
            self.ProgressBar(1, 1, "Test")
        self.ProgressBar(1, 1, 2)
        self.assertAlmostEqual(self.ProgressBar.pbar.percentage, 50.0)


if __name__ == "__main__":
    unittest.main()
