"""Unit tests of posydon/binary_evol/DT/step_detached.py
"""

__authors__ = [
    "Anurag Sawarkar <anurag.sawarkar@unige.ch>"
]

from types import SimpleNamespace

from pytest import mark

import posydon.binary_evol.DT.step_detached as totest
from posydon.binary_evol.simulationproperties import SimulationProperties


class TestOutputTimes:
    """Test for 'output_times' history rows for detached_step
    """

    @mark.parametrize("dt, output_times, expected", [
        (None, None, [10.0]),
        (None, [1.0, 4.0, 10.0, 12.0], [4.0, 10.0]),
        (4.0, [6.0, 9.0], [6.0, 9.0, 10.0]),
    ])
    def test_get_time_after_evo(self, dt, output_times, expected):
        # Step from 2 to 10 years
        step = SimpleNamespace(dt=dt, n_o_steps_history=None,
                               res=SimpleNamespace(t=[2.0, 10.0]))
        binary = SimpleNamespace(
            time=2.0, properties=SimulationProperties(output_times=output_times))
        t = totest.detached_step.get_time_after_evo(step, binary)
        assert t.tolist() == expected
