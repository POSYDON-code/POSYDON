import os
import unittest
from unittest.mock import patch

from posydon.config import PATH_TO_POSYDON
from posydon.grids.psygrid import PSyGrid

PATH_TO_GRID = os.path.join(
    PATH_TO_POSYDON,
    "posydon/tests/data/POSYDON-UNIT-TESTS/"
    "visualization/grid_unit_test_plot.h5"
)

if not os.path.exists(PATH_TO_GRID):
    raise ValueError("Test grid for unit testing was not found!")


class TestPlot1D(unittest.TestCase):
    def test_one_track_one_var_plotting(self):
        grid = PSyGrid(PATH_TO_GRID)
        with patch('matplotlib.pyplot.show') as show_patch:
            grid.plot(42,
                      "star_age",
                      "center_he4",
                      history="history1",
                      **{'show_fig': True})
            assert show_patch.called
        with patch('matplotlib.pyplot.show') as show_patch:
            grid.plot(42,
                      "age",
                      "star_1_mass",
                      history="binary_history",
                      **{'show_fig': True})
            assert show_patch.called

    def test_one_track_many_vars_plotting(self):
        grid = PSyGrid(PATH_TO_GRID)
        with patch('matplotlib.pyplot.show') as show_patch:
            grid.plot(42,
                      "star_age", ["center_he4", "log_LHe"],
                      history="history1",
                      **{'show_fig': True})
            assert show_patch.called
        with patch('matplotlib.pyplot.show') as show_patch:
            grid.plot(
                42,
                "age",
                ["star_1_mass", "binary_separation", "rl_relative_overflow_1"],
                history="binary_history",
                **{'show_fig': True})
            assert show_patch.called

    def test_many_tracks_one_var_plotting(self):
        grid = PSyGrid(PATH_TO_GRID)
        with patch('matplotlib.pyplot.show') as show_patch:
            grid.plot([42, 43, 44],
                      "star_age",
                      "center_he4",
                      history="history1",
                      **{'show_fig': True})
            assert show_patch.called

    def test_many_tracks_many_vars_plotting(self):
        grid = PSyGrid(PATH_TO_GRID)
        with patch('matplotlib.pyplot.show') as show_patch:
            grid.plot(
                [42, 43, 44],
                "age",
                ["star_1_mass", "binary_separation", "rl_relative_overflow_1"],
                history="binary_history",
                **{'show_fig': True})
            assert show_patch.called

    def test_one_track_one_var_extra_var_color_plotting(self):
        grid = PSyGrid(PATH_TO_GRID)
        with patch('matplotlib.pyplot.show') as show_patch:
            grid.plot(42,
                      "age",
                      "star_1_mass",
                      "period_days",
                      history="binary_history",
                      **{'show_fig': True})
            assert show_patch.called

    def test_many_tracks_one_var_extra_var_color_plotting(self):
        grid = PSyGrid(PATH_TO_GRID)
        with patch('matplotlib.pyplot.show') as show_patch:
            grid.plot([42, 43],
                      "age",
                      "star_1_mass",
                      "period_days",
                      history="binary_history",
                      **{'show_fig': True})
            assert show_patch.called

    def test_one_track_HR_plotting(self):
        grid = PSyGrid(PATH_TO_GRID)
        with patch('matplotlib.pyplot.show') as show_patch:
            grid.HR(42, history="history1", **{'show_fig': True})
            assert show_patch.called

    def test_many_tracks_HR_plotting(self):
        grid = PSyGrid(PATH_TO_GRID)
        with patch('matplotlib.pyplot.show') as show_patch:
            grid.HR([42, 43, 44], history="history1", **{'show_fig': True})
            assert show_patch.called


if __name__ == "__main__":
    unittest.main()
