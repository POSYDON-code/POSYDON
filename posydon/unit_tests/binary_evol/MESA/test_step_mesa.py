"""Unit tests of posydon/binary_evol/MESA/step_mesa.py
"""

__authors__ = [
    "Anurag Sawarkar <anurag.sawarkar@unige.ch>"
]

# import module which will be tested
import posydon.binary_evol.MESA.step_mesa as totest
np = totest.np

from types import SimpleNamespace

class TestOutputTimes:
    """ Tests for the helper functions relating to 'output_times' rows
    """

    def test_resample_histories(self):
        binary_history = np.array([(0.0, 1.0), (10.0, 2.0), (20.0, 4.0)],
                                  dtype=[('age', float), ('period_days', float)])
        history1 = np.array([(0.0,), (1.0,), (2.0,), (9.0,)],
                        dtype=[('center_h1', float)])
        run = SimpleNamespace(binary_history=binary_history,
                              history1=history1, history2=None)
        n_binary_history, star_histories = totest.resample_histories(
            run, [5.0, 15.0])
        assert n_binary_history['age'].tolist() == [5.0, 15.0]
        assert n_binary_history['period_days'].tolist() == [1.5, 3.0]
        assert star_histories[0]['center_h1'].tolist() == [0.5, 1.5]
        assert star_histories[1] is None

    def test_get_and_set_state(self):
        star_1 = SimpleNamespace(mass=10.0, mass_history=[10.0])
        star_2 = SimpleNamespace(mass=5.0, mass_history=[5.0])
        binary = SimpleNamespace(star_1=star_1, star_2=star_2, time=0.0,
                                 time_history=[0.0])
        state = totest._get_state(binary)
        # Change values
        binary.time = 3.0
        binary.time_history.append(1.0)
        star_1.mass = 8.0
        star_1.mass_history.append(9.0)
        star_2.new_attribute = 'new'
        totest._set_state(state)
        # Value restored, history kept
        assert binary.time == 0.0
        assert star_1.mass == 10.0
        assert not hasattr(star_2, 'new_attribute')
        assert binary.time_history == [0.0, 1.0]
        assert star_1.mass_history == [10.0, 9.0]