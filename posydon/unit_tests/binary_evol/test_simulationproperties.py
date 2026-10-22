"""Unit tests of posydon/binary_evol/simulationproperties.py
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>"
]

# import other needed code for the tests, which is not already imported in the
# module you like to test
import os

from pytest import fixture, mark

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
