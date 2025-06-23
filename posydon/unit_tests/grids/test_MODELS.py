"""Unit tests of posydon/grids/SN_MODELS.py
"""

__authors__ = [
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>"
]

# import the module which will be tested
import posydon.grids.SN_MODELS as totest

# import other needed code for the tests, which is not already imported in the
# module you like to test
from pytest import fixture, raises
from inspect import isroutine


# define test classes collecting several test functions
class TestElements:
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = {'SN_MODELS', 'NEUTRINO_MASS_LOSS_UPPER_LIMIT',\
                    'STATE_NS_STARMASS_UPPER_LIMIT', '__authors__',\
                    '__builtins__', '__cached__', '__doc__', '__file__',\
                    '__loader__', '__name__', '__package__', '__spec__',\
                    'get_SN_MODEL_NAME', 'DEFAULT_SN_MODEL', 'get_SN_MODEL'}
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

    def test_instance_DEFAULT_SN_MODEL(self):
        assert isinstance(totest.DEFAULT_SN_MODEL, dict), "DEFAULT_SN_MODEL "\
               + "is of type: " + str(type(totest.DEFAULT_SN_MODEL))

    def test_instance_SN_MODELS(self):
        assert isinstance(totest.SN_MODELS, dict), "SN_MODELS is of type: "\
               + str(type(totest.SN_MODELS))

    def test_instance_get_SN_MODEL(self):
        assert isroutine(totest.get_SN_MODEL)

    def test_instance_get_SN_MODEL_NAME(self):
        assert isroutine(totest.get_SN_MODEL_NAME)


class TestValues:
    # check that the values fit
    def test_value_DEFAULT_SN_MODEL(self):
        for k in ['mechanism', 'engine', 'PISN', 'PISN_CO_shift',\
                  'PPI_extra_mass_loss', 'ECSN', 'conserve_hydrogen_envelope',\
                  'conserve_hydrogen_PPI', 'max_neutrino_mass_loss',\
                  'max_NS_mass', 'use_interp_values', 'use_profiles',\
                  'use_core_masses', 'allow_spin_None',\
                  'approx_at_he_depletion']:
            assert k in totest.DEFAULT_SN_MODEL.keys()

    def test_value_SN_MODELS(self):
        assert len(totest.SN_MODELS) > 0
        for m in totest.SN_MODELS:
            assert isinstance(totest.SN_MODELS[m], dict)


class TestFunctions:
    # test functions
    def test_get_SN_MODEL(self):
        # missing argument
        with raises(TypeError, match="missing 1 required positional "\
                                     +"argument: 'name'"):
            totest.get_SN_MODEL()
        # examples: undefined
        assert totest.get_SN_MODEL("Test") == totest.DEFAULT_SN_MODEL
        # examples: pre defined supernova models
        for n, m in totest.SN_MODELS.items():
            model = totest.get_SN_MODEL(n)
            for k in totest.DEFAULT_SN_MODEL.keys():
                assert k in model
            for k, v in m.items():
                assert model[k] == v

    def test_get_SN_MODEL_NAME(self, capsys):
        # missing argument
        with raises(TypeError, match="missing 1 required positional "\
                                     +"argument: 'input_SN_MODEL'"):
            totest.get_SN_MODEL_NAME()
        # bad input
        with raises(TypeError, match="argument of type 'NoneType' is not "\
                                     +"iterable"):
            totest.get_SN_MODEL_NAME(None)
        # examples: missing key
        bad_SN_MODEL = {}
        for v in [True, False]:
            assert totest.get_SN_MODEL_NAME(bad_SN_MODEL, verbose=v) is None
            captured_output = capsys.readouterr().out
            if v:
                assert "missing key: mechanism" in captured_output
            else:
                assert captured_output == ""
        # examples: wrong value
        bad_SN_MODEL = {"mechanism": ""}
        for v in [True, False]:
            assert totest.get_SN_MODEL_NAME(bad_SN_MODEL, verbose=v) is None
            captured_output = capsys.readouterr().out
            if v:
                assert "mismatch: mechanism" in captured_output
            else:
                assert captured_output == ""
        # examples: pre-defined supernova models
        for n in totest.SN_MODELS.keys():
            try:
                m = totest.get_SN_MODEL(n)
            except: # skip test as test on get_SN_MODEL should fail
                assert n in totest.SN_MODELS
                return
            for v in [True, False]:
                assert totest.get_SN_MODEL_NAME(m, verbose=v) == n
                captured_output = capsys.readouterr().out
                if v:
                    assert "mismatch:" in captured_output
                    assert f"matched to supernova model: {n}" in\
                           captured_output
                else:
                    assert captured_output == ""
        # examples: first model with allowed differences
        n = list(totest.SN_MODELS.keys())[0]
        try:
            m = totest.get_SN_MODEL(n)
        except: # skip test as test on get_SN_MODEL should fail
            assert n in totest.SN_MODELS
            return
        for k in ["use_interp_values", "use_profiles", "use_core_masses"]:
            m[k] = not m[k]
            assert totest.get_SN_MODEL_NAME(m) == n
        try:
            m = totest.get_SN_MODEL(n)
        except: # skip test as test on get_SN_MODEL should fail
            assert n in totest.SN_MODELS
            return
        m["use_test"] = "unit"
        assert totest.get_SN_MODEL_NAME(m) == n
        try:
            m = totest.get_SN_MODEL(n)
        except: # skip test as test on get_SN_MODEL should fail
            assert n in totest.SN_MODELS
            return
        m["ECSN"] = "test"
        assert totest.get_SN_MODEL_NAME(m) == n
