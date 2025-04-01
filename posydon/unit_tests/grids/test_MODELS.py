"""Unit tests of posydon/grids/MODELS.py
"""

__authors__ = [
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>"
]

# import the module which will be tested
import posydon.grids.MODELS as totest

# import other needed code for the tests, which is not already imported in the
# module you like to test
from pytest import fixture, raises
from inspect import isroutine


# define test classes collecting several test functions
class TestElements:
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = {'MODELS', 'NEUTRINO_MASS_LOSS_UPPER_LIMIT',\
                    'STATE_NS_STARMASS_UPPER_LIMIT', '__authors__',\
                    '__builtins__', '__cached__', '__doc__', '__file__',\
                    '__loader__', '__name__', '__package__', '__spec__',\
                    'get_MODEL_NAME'}
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

    def test_instance_MODELS(self):
        assert isinstance(totest.MODELS, dict), "MODELS is of type: "\
                                                +str(type(totest.MODELS))

    def test_instance_get_MODEL_NAME(self):
        assert isroutine(totest.get_MODEL_NAME)


class TestValues:
    # check that the values fit
    def test_value_MODELS(self):
        assert len(totest.MODELS) > 0
        for m in totest.MODELS:
            assert isinstance(totest.MODELS[m], dict)
            assert 'mechanism' in totest.MODELS[m].keys()
            assert 'engine' in totest.MODELS[m].keys()


class TestFunctions:
    @fixture
    def bad_MODEL(self):
        return {"mechanism": "",\
                "engine": "",\
                "PISN": "",\
                "ECSN": "",\
                "conserve_hydrogen_envelope": False,\
                "max_neutrino_mass_loss":\
                 totest.NEUTRINO_MASS_LOSS_UPPER_LIMIT,\
                "max_NS_mass": totest.STATE_NS_STARMASS_UPPER_LIMIT,\
                "use_interp_values": False,\
                "use_profiles": False,\
                "use_core_masses": False,\
                "approx_at_he_depletion": False}

    # test functions
    def test_get_MODEL_NAME(self, bad_MODEL, capsys):
        # missing argument
        with raises(TypeError, match="missing 1 required positional "\
                                     +"argument: 'input_MODEL'"):
            totest.get_MODEL_NAME()
        # bad input
        with raises(TypeError, match="'NoneType' object is not subscriptable"):
            totest.get_MODEL_NAME(None)
        # examples: failed
        assert totest.get_MODEL_NAME(bad_MODEL) is None
        # examples: pre-defined models
        for n, m in totest.MODELS.items():
            for v in [True, False]:
                assert totest.get_MODEL_NAME(m, verbose=v) == n
                captured_output = capsys.readouterr().out
                if v:
                    assert "mismatch:" in captured_output
                    assert f"matched to model: {n}" in captured_output
                else:
                    assert captured_output == ""
        # examples: first model with allowed differences
        n = list(totest.MODELS.keys())[0]
        for k in ["use_interp_values", "use_profiles", "use_core_masses"]:
            m = totest.MODELS[n].copy()
            m[k] = not m[k]
            assert totest.get_MODEL_NAME(m) == n
        m = totest.MODELS[n].copy()
        m["use_test"] = "unit"
        assert totest.get_MODEL_NAME(m) == n
        m = totest.MODELS[n].copy()
        m["ECSN"] = "test"
        assert totest.get_MODEL_NAME(m) == n
