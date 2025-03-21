"""Unit tests of posydon/utils/datasets.py

"""

__authors__ = [
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>"
]

# import the module which will be tested
import posydon.utils.datasets as totest

# import other needed code for the tests, which is not already imported in the
# module you like to test


# define test classes collecting several test functions
class TestElements:
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = ['COMPLETE_SETS', 'ZENODO_COLLECTION', '__authors__',\
                    '__builtins__', '__cached__', '__doc__', '__file__',\
                    '__loader__', '__name__', '__package__', '__spec__']
        assert dir(totest) == elements, "There might be added or removed "\
                                        + "objects without an update on the "\
                                        + "unit test."

    def test_instance_ZENODO_COLLECTION(self):
        assert isinstance(totest.ZENODO_COLLECTION, dict)

    def test_instance_COMPLETE_SETS(self):
        assert isinstance(totest.COMPLETE_SETS, dict)


class TestValues:
    # check that the values fit
    def test_value_ZENODO_COLLECTION(self):
        # check datasets
        for k,s in totest.ZENODO_COLLECTION.items():
            assert isinstance(s, dict), f"ZENODO_COLLECTION['{k}'] should be "\
                                        + f"a dictionary, but is {type(s)}"
            # check required entries for each dataset
            for e in ['data', 'md5']:
                # entry is string or None (for not yet published datasets)
                assert e in s, f"The '{e}' entry is missing for "\
                               + f"ZENODO_COLLECTION['{k}']."
                assert isinstance(s[e], (str, type(None))),\
                       f"ZENODO_COLLECTION['{k}']['{e}'] should be a string "\
                       + f"or None, but is {type(s[e])}"
            for e in ['description', 'title', 'url']:
                # entry is string
                assert e in s, f"The '{e}' entry is missing for "\
                               + f"ZENODO_COLLECTION[{k}]."
                assert isinstance(s[e], str),\
                       f"ZENODO_COLLECTION['{k}']['{e}'] should be a string, "\
                       + f"but is {type(s[e])}"

    def test_value_COMPLETE_SETS(self):
        # check base versions
        for k in ['v1', 'v2']:
            assert k in totest.COMPLETE_SETS
        # check complete datasets
        for k,s in totest.COMPLETE_SETS.items():
            assert isinstance(s, list), f"COMPLETE_SETS['{k}'] should be a "\
                                        + f"list, but is {type(s)}"
            # check that individual datasets are defined
            for v in s:
                assert v in totest.ZENODO_COLLECTION.keys(),\
                       f"'{v}' is an unknown dataset."
