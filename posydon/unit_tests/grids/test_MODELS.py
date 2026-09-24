"""Unit tests of posydon/grids/SN_MODELS.py
"""

__authors__ = [
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>"
]

from inspect import isroutine

# import other needed code for the tests, which is not already imported in the
# module you like to test
from pytest import fixture, raises

# import the module which will be tested
import posydon.grids.SN_MODELS as totest


# define test classes collecting several test functions
class TestElements:
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = {'SN_MODELS', 'NEUTRINO_MASS_LOSS_UPPER_LIMIT',\
                    'STATE_NS_STARMASS_UPPER_LIMIT', '__authors__',\
                    '__builtins__', '__cached__', '__doc__', '__file__',\
                    '__loader__', '__name__', '__package__', '__spec__',\
                    'get_SN_MODEL_NAME', 'DEFAULT_SN_MODEL',\
                    'get_SN_MODEL', 'MALTSEV_MECHANISMS',\
                    'DEFAULT_MALTSEV_SN_MODEL', 'get_SN_MODEL_parameters',\
                    'missing_SN_MODEL_parameters', 'check_SN_MODELS'}
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

    def test_instance_MALTSEV_MECHANISMS(self):
        assert isinstance(totest.MALTSEV_MECHANISMS, tuple),\
               "MALTSEV_MECHANISMS is of type: "\
               + str(type(totest.MALTSEV_MECHANISMS))

    def test_instance_DEFAULT_MALTSEV_SN_MODEL(self):
        assert isinstance(totest.DEFAULT_MALTSEV_SN_MODEL, dict),\
               "DEFAULT_MALTSEV_SN_MODEL is of type: "\
               + str(type(totest.DEFAULT_MALTSEV_SN_MODEL))

    def test_instance_get_SN_MODEL_parameters(self):
        assert isroutine(totest.get_SN_MODEL_parameters)

    def test_instance_missing_SN_MODEL_parameters(self):
        assert isroutine(totest.missing_SN_MODEL_parameters)

    def test_instance_check_SN_MODELS(self):
        assert isroutine(totest.check_SN_MODELS)


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

    def test_Maltsev_parameters_not_in_DEFAULT_SN_MODEL(self):
        # a DEFAULT_SN_MODEL parameter must be present in every model and ini
        # file to match, so the Maltsev+25 extras must stay out of it
        for k in totest.DEFAULT_MALTSEV_SN_MODEL.keys():
            assert k not in totest.DEFAULT_SN_MODEL.keys()

    def test_value_DEFAULT_MALTSEV_SN_MODEL(self):
        for k in ['Maltsev25_MCO_NS_mass', 'Maltsev25_MCO_fallback_fraction',\
                  'Maltsev25_MCO_fallback_model',\
                  'Maltsev25_MCO_extrapolation_mode']:
            assert k in totest.DEFAULT_MALTSEV_SN_MODEL.keys()
        assert totest.DEFAULT_MALTSEV_SN_MODEL[
            'Maltsev25_MCO_extrapolation_mode'] == "balanced"

    def test_value_MALTSEV_MECHANISMS(self):
        for m in ["Maltsev+25-engine", "Maltsev+25-MCO-rapid"]:
            assert m in totest.MALTSEV_MECHANISMS

    def test_rapid_balanced_models(self):
        rapid_balanced = [n for n in totest.SN_MODELS
                          if totest.SN_MODELS[n].get('mechanism') ==
                          "Maltsev+25-MCO-rapid" and
                          totest.SN_MODELS[n].get(
                              'Maltsev25_MCO_extrapolation_mode') == "balanced"]
        assert len(rapid_balanced) == 4
        for n in rapid_balanced:
            model = totest.get_SN_MODEL(n)
            assert model['mechanism'] == "Maltsev+25-MCO-rapid"
            assert model['Maltsev25_MCO_extrapolation_mode'] == "balanced"
            # every Maltsev parameter is set by the model itself
            assert totest.missing_SN_MODEL_parameters(
                totest.SN_MODELS[n]) == []
            # round-trip back to the model name
            assert totest.get_SN_MODEL_NAME(model) == n

    def test_balanced_does_not_match_optimistic(self):
        n = "SN_MODEL_v2_29"
        m = totest.get_SN_MODEL(n)
        # flipping the extrapolation mode to optimistic must no longer match
        # the balanced model
        m2 = dict(m)
        m2['Maltsev25_MCO_extrapolation_mode'] = "optimistic"
        # must not match the balanced rapid model; could match another model
        # or be None, but not n
        assert totest.get_SN_MODEL_NAME(m2) != n

    def test_get_SN_MODEL_parameters(self):
        # a mechanism without extra parameters is defined by the default model
        assert totest.get_SN_MODEL_parameters("Fryer+12-delayed") ==\
               list(totest.DEFAULT_SN_MODEL.keys())
        # the Maltsev+25 mechanisms add their own parameters on top
        for m in totest.MALTSEV_MECHANISMS:
            parameters = totest.get_SN_MODEL_parameters(m)
            assert parameters == list(totest.DEFAULT_SN_MODEL.keys())\
                   + list(totest.DEFAULT_MALTSEV_SN_MODEL.keys())

    def test_missing_SN_MODEL_parameters(self):
        # nothing is required beyond the defaults for other mechanisms
        assert totest.missing_SN_MODEL_parameters(
            {'mechanism': "Fryer+12-delayed"}) == []
        # ...but the Maltsev+25 mechanisms require all of their parameters
        for m in totest.MALTSEV_MECHANISMS:
            assert totest.missing_SN_MODEL_parameters({'mechanism': m}) ==\
                   list(totest.DEFAULT_MALTSEV_SN_MODEL.keys())
            complete = {'mechanism': m, **totest.DEFAULT_MALTSEV_SN_MODEL}
            assert totest.missing_SN_MODEL_parameters(complete) == []

    def test_check_SN_MODELS(self):
        # the pre-defined models are complete
        totest.check_SN_MODELS()
        # an incomplete Maltsev model is rejected
        totest.SN_MODELS['SN_MODEL_unit_test'] = {'mechanism':
                                                  "Maltsev+25-engine"}
        try:
            with raises(ValueError, match="SN_MODEL_unit_test"):
                totest.check_SN_MODELS()
        finally:
            del totest.SN_MODELS['SN_MODEL_unit_test']

    def test_model_matching_without_Maltsev_parameters(self):
        # a non-Maltsev ini carries no Maltsev parameters, but must still
        # match a model, or the popsyn setup rejects every existing ini file
        input_SN_MODEL = totest.DEFAULT_SN_MODEL.copy()
        for k in totest.DEFAULT_MALTSEV_SN_MODEL.keys():
            assert k not in input_SN_MODEL
        assert totest.get_SN_MODEL_NAME(input_SN_MODEL) is not None

    def test_model_matching_needs_Maltsev_parameters(self):
        # a Maltsev model is only identified with all of its parameters given
        n = "SN_MODEL_v2_29"
        model = totest.get_SN_MODEL(n)
        assert totest.get_SN_MODEL_NAME(model) == n
        for k in totest.DEFAULT_MALTSEV_SN_MODEL.keys():
            incomplete = {key: v for key, v in model.items() if key != k}
            assert totest.get_SN_MODEL_NAME(incomplete) is None
        # ...and each of them is part of the model identity
        for k, v in [('Maltsev25_MCO_NS_mass', 1.5),
                     ('Maltsev25_MCO_fallback_fraction', 0.5),
                     ('Maltsev25_MCO_fallback_model', "B"),
                     ('Maltsev25_MCO_extrapolation_mode', "optimistic")]:
            changed = dict(model)
            changed[k] = v
            assert totest.get_SN_MODEL_NAME(changed) != n

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
