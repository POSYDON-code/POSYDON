"""Unit tests of posydon/binary_evol/simulationproperties.py
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
    "Anurag Sawarkar <anurag.sawarkar@unige.ch>"
]

# import other needed code for the tests, which is not already imported in the
# module you like to test
import os
from types import SimpleNamespace

from pytest import fixture, mark, raises

# import the module which will be tested
import posydon.binary_evol.simulationproperties as totest
from posydon.config import PATH_TO_POSYDON_DATA


class FakeGridStep:
    """Minimal step that, like a MesaGridStep, takes a grid and metallicity."""
    DEFAULT_KWARGS = {'metallicity': None, 'grid_path': None}

    def __init__(self, **kwargs):
        self.kwargs = kwargs


class TestGridPaths:
    @fixture
    def sim_prop(self):
        return totest.SimulationProperties()

    @mark.parametrize("step_name, path_name, grid_dir", [
        ("step_HMS_HMS", "HMS_HMS_path", "HMS-HMS"),
        ("step_HMS_HMS_RLO", "HMS_HMS_RLO_path", "HMS-HMS_RLO"),
        ("step_CO_HMS_RLO", "CO_HMS_RLO_path", "CO-HMS_RLO"),
        ("step_CO_HeMS", "CO_HeMS_path", "CO-HeMS"),
        ("step_CO_HeMS_RLO", "CO_HeMS_RLO_path", "CO-HeMS_RLO"),
    ])
    def test_default_grid_path(self, sim_prop, step_name, path_name,
                               grid_dir):
        expected = os.path.abspath(os.path.join(PATH_TO_POSYDON_DATA,
                                                grid_dir))
        assert getattr(sim_prop, path_name) == expected
        _, step_kwargs = sim_prop.check_step(0.1, None, step_name,
                                             (FakeGridStep, {}))
        assert step_kwargs['grid_path'] == expected
        assert step_kwargs['metallicity'] == 0.1

    def test_custom_HMS_HMS_RLO_path(self, tmp_path):
        sim_prop = totest.SimulationProperties(
            HMS_HMS_RLO_path=str(tmp_path))
        _, step_kwargs = sim_prop.check_step(1.0, None, "step_HMS_HMS_RLO",
                                             (FakeGridStep, {}))
        assert step_kwargs['grid_path'] == str(tmp_path)


class TestOutputTimes:
    """Test the output_times property read/write and the output_times_between method."""
    def test_output_times(self):
        assert len(totest.SimulationProperties().output_times) == 0
        sim_prop = totest.SimulationProperties(output_times=[3.0, 1.0, 3.0, 2.0])
        assert sim_prop.output_times.tolist() == [1.0, 2.0, 3.0]
        assert sim_prop.output_times_between(1.0, 3.0).tolist() == [2.0]
        assert sim_prop.output_times_between(3.0, 6.0).tolist() == []
        with raises(ValueError, match="non-finite"):
            sim_prop.output_times = [1.0, float('nan')]


class TestHooksOutputTimes:
    """Test the output_times timesteps' hooks for names and times"""
    @fixture
    def binary(self):
        return SimpleNamespace(event_history=['ZAMS', None, None, 'oRLO1'],
                               step_names=['initial_cond'],
                               step_times=[0.0])
    def test_step_names(self, binary):
        totest.StepNamesHooks().post_step(binary, 'step_HMS_HMS')
        assert binary.step_names == ['initial_cond'] + 3 * ['step_HMS_HMS']

    def test_step_times(self, binary):
        hooks = totest.TimingHooks()
        hooks.pre_step(binary, 'step_HMS_HMS')
        hooks.post_step(binary, 'step_HMS_HMS')
        assert binary.step_times[:3] == [0.0, None, None]
        assert binary.step_times[3] >= 0.0
