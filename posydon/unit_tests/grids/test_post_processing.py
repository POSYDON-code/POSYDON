"""Unit tests of posydon/grids/post_processing.py
"""

__authors__ = [
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>"
]

# import the module which will be tested
import posydon.grids.post_processing as totest
# aliases
np = totest.np

# import other needed code for the tests, which is not already imported in the
# module you like to test
from pytest import fixture, raises, warns
from inspect import isclass, isroutine
from posydon.utils.posydonwarning import InappropriateValueWarning


# define test classes collecting several test functions
class TestElements:
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = ['BinaryStar', 'CC_quantities',\
                    'CEE_parameters_from_core_abundance_thresholds',\
                    'Catch_POSYDON_Warnings',\
                    'DEFAULT_MARKERS_COLORS_LEGENDS', 'MODELS', 'Pwarn',\
                    'STAR_STATES_CC', 'SingleStar', 'StepSN',\
                    'TF1_POOL_STABLE', '__authors__', '__builtins__',\
                    '__cached__', '__credits__', '__doc__', '__file__',\
                    '__loader__', '__name__', '__package__', '__spec__',\
                    'add_post_processed_quantities',\
                    'assign_core_collapse_quantities_none',\
                    'calculate_Patton20_values_at_He_depl',\
                    'check_state_of_star', 'combine_TF12', 'copy', 'np',\
                    'post_process_grid', 'print_CC_quantities', 'tqdm']
        assert dir(totest) == elements, "There might be added or removed "\
                                        + "objects without an update on the "\
                                        + "unit test."

    def test_instance_CC_quantities(self):
        assert isinstance(totest.CC_quantities, (list)), "CC_quantities is "\
               + "of type: " + str(type(totest.CC_quantities))

    def test_instance_assign_core_collapse_quantities_none(self):
        assert isroutine(totest.assign_core_collapse_quantities_none)

    def test_instance_print_CC_quantities(self):
        assert isroutine(totest.print_CC_quantities)

    def test_instance_post_process_grid(self):
        assert isroutine(totest.post_process_grid)

    def test_instance_add_post_processed_quantities(self):
        assert isroutine(totest.add_post_processed_quantities)


class TestValues:
    # check that the values fit
    def test_value_CC_quantities(self):
        for q in ['state', 'SN_type', 'f_fb', 'mass', 'spin',\
                  'm_disk_accreted', 'm_disk_radiated',\
                  'CO_interpolation_class', 'M4', 'mu4', 'h1_mass_ej',\
                  'he4_mass_ej']:
            assert q in totest.CC_quantities


class TestFunctions:
    @fixture
    def star(self):
    # initialize a SingleStar instance, which is a required argument
        test_star = totest.SingleStar()
        test_star.state = "TestState"
        test_star.mass = 1.0
        test_star.spin = 0.0
        test_star.SN_type = "TestSNType"
        test_star.f_fb = 0.1
        test_star.m_disk_accreted = 0.2
        test_star.m_disk_radiated = 0.3
        test_star.M4 = 0.4
        test_star.mu4 = 0.5
        test_star.h1_mass_ej = np.nan
        test_star.he4_mass_ej = None
        return test_star

    # test functions
    def test_assign_core_collapse_quantities_none(self):
        # missing argument
        with raises(TypeError, match="missing 2 required positional "\
                                     +"arguments: 'EXTRA_COLUMNS' and "\
                                     +"'star_i'"):
            totest.assign_core_collapse_quantities_none()
        # bad input
        with raises(TypeError, match="'NoneType' object is not subscriptable"):
            totest.assign_core_collapse_quantities_none(None, 1)
        # bad input
        test_EXTRA_COLUMNS = {}
        with raises(ValueError, match="'star_i' should be 1 or 2."):
            totest.assign_core_collapse_quantities_none(test_EXTRA_COLUMNS, 0)
        # bad input
        with raises(TypeError, match="'MODEL_NAME' should be a string or "\
                                     +"None."):
            totest.assign_core_collapse_quantities_none(test_EXTRA_COLUMNS, 1,\
                                                        1)
        # bad input
        with raises(KeyError, match="'S1_MODEL01_state'"):
            totest.assign_core_collapse_quantities_none(test_EXTRA_COLUMNS, 1)
        # examples: all models in MODELS.py
        for s in [1, 2]:
            test_EXTRA_COLUMNS = {}
            for m in totest.MODELS.keys():
                for q in totest.CC_quantities:
                    test_EXTRA_COLUMNS[f'S{s}_{m}_{q}'] = []
            totest.assign_core_collapse_quantities_none(test_EXTRA_COLUMNS, s)
            for m in totest.MODELS.keys():
                for q in totest.CC_quantities:
                    assert test_EXTRA_COLUMNS[f'S{s}_{m}_{q}'] == [None]
        # examples: given model
        for s in [1, 2]:
            test_EXTRA_COLUMNS = {}
            for m in ["TESTMODEL"]:
                for q in totest.CC_quantities:
                    test_EXTRA_COLUMNS[f'S{s}_{m}_{q}'] = []
            totest.assign_core_collapse_quantities_none(test_EXTRA_COLUMNS, s,\
                                                        MODEL_NAME=m)
            for m in ["TESTMODEL"]:
                for q in totest.CC_quantities:
                    assert test_EXTRA_COLUMNS[f'S{s}_{m}_{q}'] == [None]

    def test_print_CC_quantities(self, star, capsys):
        # missing argument
        with raises(TypeError, match="missing 1 required positional "\
                                     +"argument: 'star'"):
            totest.print_CC_quantities()
        # bad input
        with warns(InappropriateValueWarning) as winfo:
            totest.print_CC_quantities(None)
        assert "Failed to print star values!\nWarning in preSN: 'NoneType' "\
               + "object has no attribute 'state'"\
               == winfo._list[0].message.args[0]
        output = capsys.readouterr().out.split("\n")
        assert len(output) == 4+1
        for h in ["mechanism", "state", "SN type", "f_fb", "mass [Msun]",\
                  "spin", "m_disk_accreted [Msun]", "m_disk_radiated [Msun]",\
                  "M4 [m/Msun]", "mu4 [(dm/Msun)/(dr/1000km)]",\
                  "h1_mass_ej [Msun]", "he4_mass_ej [Msun]"]:
            assert h in output[1]
        # bad input
        with warns(InappropriateValueWarning) as winfo:
            totest.print_CC_quantities(None, "TESTMODEL")
        assert "Failed to print star values!\nWarning in TESTMODEL: "\
               + "'NoneType' object has no attribute 'spin'"\
               == winfo._list[0].message.args[0]
        # examples: pre SN
        totest.print_CC_quantities(star)
        output = capsys.readouterr().out.split("\n")
        assert len(output) == 5+1
        for h in ["mechanism", "state", "SN type", "f_fb", "mass [Msun]",\
                  "spin", "m_disk_accreted [Msun]", "m_disk_radiated [Msun]",\
                  "M4 [m/Msun]", "mu4 [(dm/Msun)/(dr/1000km)]",\
                  "h1_mass_ej [Msun]", "he4_mass_ej [Msun]"]:
            assert h in output[1]
        for v in ["PRE SN STAR", star.state,\
                  "{:7.2f} {:12.2f}".format(star.mass, star.spin)]:
            assert v in output[3]
        # examples: SN MODEL
        totest.print_CC_quantities(star, "TESTMODEL")
        output = capsys.readouterr().out.split("\n")
        with capsys.disabled():
            print(output)
        assert len(output) == 1+1
        fmt = "{:1.2f} {:13.2f} {:12.2f} {:20.2f} {:20.2f} {:20.2f} {:20.2f} "\
              +"{:20.2f} {:20.2f}"
        for v in ["TESTMODEL", star.state, star.SN_type, fmt.format(star.f_fb,\
                   star.mass, star.spin, star.m_disk_accreted,\
                   star.m_disk_radiated, star.M4, star.mu4, star.h1_mass_ej,\
                   np.nan)]:
            assert v in output[0]

    def test_post_process_grid(self):
        # missing argument
        with raises(TypeError, match="missing 1 required positional "\
                                     +"argument: 'grid'"):
            totest.post_process_grid()
        # bad input
        with raises(AttributeError, match="'NoneType' object has no "\
                                          +"attribute 'MESA_dirs'"):
            totest.post_process_grid(None)
        # examples
        pass

    def test_add_post_processed_quantities(self):
        # missing argument
        with raises(TypeError, match="missing 3 required positional "\
                                     +"arguments: 'grid', "\
                                     +"'MESA_dirs_EXTRA_COLUMNS', and "\
                                     +"'EXTRA_COLUMNS'"):
            totest.add_post_processed_quantities()
        # bad input
        with raises(AttributeError, match="'NoneType' object has no "\
                                          +"attribute 'MESA_dirs'"):
            totest.add_post_processed_quantities(None, None, None)
        # examples
        pass
