"""Unit-testing module for the `ConfigFile` class."""

import os
import unittest
import numpy as np
from numpy.testing import assert_array_equal
from posydon.utils.common_functions import PATH_TO_POSYDON
from posydon.utils.configfile import ConfigFile


configuration = {
    "name": "Yoda",
    "age": 900,
    "dimensions": np.array([26.0, 11.1, 5.2]),
    "apprentices": ["Dooku", "Windu"],
}


class TestConfigFile(unittest.TestCase):
    """Unit-test class for the `ConfigFile` class."""

    @classmethod
    def setUpClass(cls):
        """Get path to temporary file."""
        cls.path = os.path.join(PATH_TO_POSYDON, "posydon/tests/data/tmp.cfg")

    @classmethod
    def tearDownClass(cls):
        """Ensure temporary file was deleted after tests."""
        if os.path.exists(cls.path):
            os.remove(cls.path)

    def setUp(self):
        """Create a ConfigFile object before each test."""
        self.config = ConfigFile()
        self.config.update(configuration)

    def test_IO(self):
        """Test writing/reading from disk, and serialization."""
        # save to disk and destroy the object
        self.config.save(self.path, overwrite=False)
        del self.config

        # read it and test data integrity
        self.config = ConfigFile(self.path)
        self.assertEqual(len(self.config), len(configuration))
        for key in configuration:
            original = configuration[key]
            fromdisk = self.config[key]
            if np.iterable(original):
                assert_array_equal(original, fromdisk)
            else:
                self.assertEqual(original, fromdisk)

        # delete temporary file
        os.remove(self.path)

        # check `print` and `str` methods
        print(str(self))
        print(len(str(self).split("\m")))

    def test_getters(self):
        """Test the two ways to get entries."""
        self.assertEqual(self.config.name, self.config["name"])

    def test_setters(self):
        """Test the setting of entries one-by-one."""
        self.config["name"] = "Yodah"
        self.assertEqual(self.config["name"], "Yodah")

        self.config["side"] = "light"
        self.assertTrue("side" in self.config.keys())
        del self.config["side"]
        self.assertTrue("side" not in self.config.keys())

    def test_dict_methods(self):
        """Test the `items`, `keys`, `values` methods, and iteration."""
        for key, value in self.config.items():                # check items()
            self.assertIn(key, configuration)                 # check key
            self.assertIn(key, self.config)                   # check __iter__
            self.assertIn(key, self.config.keys())            # check keys()
            if not np.iterable(value):
                self.assertIn(value, self.config.values())    # check values()

    def test_update(self):
        """Test the `update` method for altering or adding entries."""
        self.config.update({"age": 901, "starring": "Star Wars"})
        self.assertEqual(self.config.age, 901)
        self.assertEqual(self.config.starring, "Star Wars")

    def test_len(self):
        """The the `len` method applied on ConfigFile."""
        self.assertEqual(len(self.config), len(configuration))


if __name__ == "__main__":
    unittest.main()
