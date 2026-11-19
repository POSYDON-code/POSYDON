"""Unit tests of posydon/binary_evol/MESA/step_mesa.py
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>"
]

# import other needed code for the tests, which is not already imported in the
# module you like to test
from pytest import fixture, raises

# import the module which will be tested
import posydon.binary_evol.MESA.step_mesa as totest
from posydon.binary_evol.simulationproperties import SimulationProperties
from posydon.binary_evol.singlestar import SingleStar

# aliases
np = totest.np
BinaryStar = totest.BinaryStar
GridError = totest.GridError


def make_RLO2_binary():
    """A binary with a compact object as star_1 and a donor as star_2."""
    binary = BinaryStar(star_1=SingleStar(mass=20.0, state='BH'),
                        star_2=SingleStar(mass=3.0,
                                          state='H-rich_Core_H_burning'),
                        state='RLO2', event='oRLO2', orbital_period=5.0,
                        time=1.0e7)
    binary.properties = SimulationProperties(max_simulation_time=13.8e9)
    return binary


class TestMesaGridStepFlipStars:
    """MesaGridStep.__call__ must leave star_1/star_2 in their original order,
    also when the step itself fails (issue #977)."""

    @fixture
    def step(self):
        # avoid loading any grids by bypassing __init__
        step = object.__new__(totest.MesaGridStep)
        step.flip_stars_before_step = True
        step.stop_method = 'stop_at_end'
        step.stop_var_name = None
        step.stop_value = None
        step.stop_star = 'star_1'
        step.stop_interpolate = True
        step.star_1_CO = False
        step.star_2_CO = False
        step.interpolation_method = 'linear3c_kNN'
        step.verbose = False
        step.get_final_MESA_step_time = lambda: 1.0e6
        step.stop_at_end = lambda binary, **kwargs: None
        return step

    def assert_original_order(self, binary):
        assert binary.star_1.state == 'BH'
        assert binary.star_2.state == 'H-rich_Core_H_burning'
        assert binary.star_1.mass_history[0] == 20.0
        assert binary.star_2.mass_history[0] == 3.0
        assert binary.state == 'RLO2'
        assert binary.event == 'oRLO2'

    def test_step_fails(self, step):
        binary = make_RLO2_binary()
        seen = {}

        def failing_step(binary, interp_method=None):
            # the step itself sees the flipped binary
            seen['star_1_state'] = binary.star_1.state
            raise GridError("interpolation failed")

        step.step = failing_step
        with raises(GridError, match="interpolation failed"):
            step(binary)
        assert seen['star_1_state'] == 'H-rich_Core_H_burning'
        self.assert_original_order(binary)

    def test_stop_method_fails(self, step):
        binary = make_RLO2_binary()
        step.step = lambda binary, interp_method=None: None

        def failing_stop(binary, **kwargs):
            raise ValueError("stop condition failed")

        step.stop_at_end = failing_stop
        with raises(ValueError, match="stop condition failed"):
            step(binary)
        self.assert_original_order(binary)

    def test_step_succeeds(self, step):
        binary = make_RLO2_binary()
        step.step = lambda binary, interp_method=None: None
        step(binary)
        self.assert_original_order(binary)

    def test_initial_RLOF(self, step):
        binary = make_RLO2_binary()
        step.get_final_MESA_step_time = lambda: np.nan
        step.step = lambda binary, interp_method=None: None
        step(binary)
        assert binary.state == 'initial_RLOF'
        assert binary.star_1.state == 'BH'
        assert binary.star_2.state == 'H-rich_Core_H_burning'
        assert binary.star_1.mass_history[0] == 20.0

    def test_no_flip(self, step):
        binary = make_RLO2_binary()
        step.flip_stars_before_step = False

        def failing_step(binary, interp_method=None):
            raise GridError("interpolation failed")

        step.step = failing_step
        with raises(GridError):
            step(binary)
        self.assert_original_order(binary)
