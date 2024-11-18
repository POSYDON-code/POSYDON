"""Unit tests of posydon/grids/MODELS.py
"""

__authors__ = [
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>"
]

# import the module which will be tested
import posydon.grids.MODELS as totest

# import other needed code for the tests, which is not already imported in the
# module you like to test


# define test classes collecting several test functions
class TestElements:
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = ['MODELS', 'NEUTRINO_MASS_LOSS_UPPER_LIMIT',\
                    'STATE_NS_STARMASS_UPPER_LIMIT', '__authors__',\
                    '__builtins__', '__cached__', '__doc__', '__file__',\
                    '__loader__', '__name__', '__package__', '__spec__']
        assert dir(totest) == elements, "There might be added or removed "+\
               "objects without an update on the unit test."

    def test_instance_MODELS(self):
        assert isinstance(totest.MODELS, dict), "MODELS is of type: "+\
               str(type(totest.MODELS))


class TestValues:
    # check that the values fit
    def test_value_MODELS(self):
        assert len(totest.MODELS)>0
        for m in totest.MODELS:
            assert isinstance(totest.MODELS[m], dict)
            assert 'mechanism' in totest.MODELS[m].keys()
            assert 'engine' in totest.MODELS[m].keys()
