import unittest
from unittest.mock import patch
from posydon.config import PATH_TO_POSYDON
import os
from posydon.grids.psygrid import PSyGrid

PATH_TO_GRID = os.path.join(
    PATH_TO_POSYDON,
    "posydon/tests/data/POSYDON-UNIT-TESTS/"
    "visualization/grid_unit_test_plot.h5"
)
if not os.path.exists(PATH_TO_GRID):
    raise ValueError("Test grid for unit testing was not found!")

# class TestPlot2D(unittest.TestCase):
#     def test_termination_flag_1_plotting(self):
#         grid = PSyGrid(PATH_TO_GRID)
#         with patch('matplotlib.pyplot.show') as show_patch:
#             grid.plot2D("star_1_mass",
#                         "period_days",
#                         "star_1_mass",
#                         termination_flag="termination_flag_1",
#                         grid_3D=True,
#                         slice_3D_var_str="star_2_mass",
#                         slice_3D_var_range=(2.5, 3.0),
#                         **{'show_fig': True})
#             assert show_patch.called
#         with patch('matplotlib.pyplot.show') as show_patch:
#             grid.plot2D("star_1_mass",
#                         "period_days",
#                         "c_core_mass",
#                         termination_flag="termination_flag_1",
#                         grid_3D=True,
#                         slice_3D_var_str="star_2_mass",
#                         slice_3D_var_range=(2.5, 3.0),
#                         **{'show_fig': True})
#             assert show_patch.called
#         with patch('matplotlib.pyplot.show') as show_patch:
#             grid.plot2D("star_1_mass",
#                         "period_days",
#                         "binary_separation",
#                         termination_flag="termination_flag_1",
#                         grid_3D=True,
#                         slice_3D_var_str="star_2_mass",
#                         slice_3D_var_range=(2.5, 3.0),
#                         **{'show_fig': True})
#             assert show_patch.called
#
#     def test_termination_flag_2_plotting(self):
#         grid = PSyGrid(PATH_TO_GRID)
#         with patch('matplotlib.pyplot.show') as show_patch:
#             grid.plot2D("star_1_mass",
#                         "period_days",
#                         None,
#                         termination_flag="termination_flag_2",
#                         grid_3D=True,
#                         slice_3D_var_str="star_2_mass",
#                         slice_3D_var_range=(2.5, 3.0),
#                         **{'show_fig': True})
#             assert show_patch.called
#
#     def test_termination_flag_3_plotting(self):
#         grid = PSyGrid(PATH_TO_GRID)
#         with patch('matplotlib.pyplot.show') as show_patch:
#             grid.plot2D("star_1_mass",
#                         "period_days",
#                         None,
#                         termination_flag="termination_flag_3",
#                         grid_3D=True,
#                         slice_3D_var_str="star_2_mass",
#                         slice_3D_var_range=(2.5, 3.0),
#                         **{'show_fig': True})
#             assert show_patch.called
#
#     def test_termination_flag_4_plotting(self):
#         grid = PSyGrid(PATH_TO_GRID)
#         with patch('matplotlib.pyplot.show') as show_patch:
#             grid.plot2D("star_1_mass",
#                         "period_days",
#                         None,
#                         termination_flag="termination_flag_4",
#                         grid_3D=True,
#                         slice_3D_var_str="star_2_mass",
#                         slice_3D_var_range=(2.5, 3.0),
#                         **{'show_fig': True})
#             assert show_patch.called
#
#     def test_all_termination_flags_plotting(self):
#         grid = PSyGrid(PATH_TO_GRID)
#         with patch('matplotlib.pyplot.show') as show_patch:
#             grid.plot2D("star_1_mass",
#                         "period_days",
#                         "binary_separation",
#                         termination_flag="all",
#                         grid_3D=True,
#                         slice_3D_var_str="star_2_mass",
#                         slice_3D_var_range=(2.5, 3.0),
#                         **{'show_fig': True})
#             assert show_patch.called
#
#     def test_RLO_plotting(self):
#         grid = PSyGrid(PATH_TO_GRID)
#         # with patch('matplotlib.pyplot.show') as show_patch:
#         #     grid.plot2D("star_1_mass",
#         #                 "period_days",
#         #                 "c_core_mass",
#         #                 termination_flag="termination_flag_1",
#         #                 grid_3D=True,
#         #                 slice_3D_var_str="star_2_mass",
#         #                 slice_3D_var_range=(2.5, 3.0),
#         #                 slice_at_RLO=True,
#         #                 **{
#         #                 'show_fig': True
#         #                 })
#         #     assert show_patch.called
#         # with patch('matplotlib.pyplot.show') as show_patch:
#         #     grid.plot2D("star_1_mass",
#         #                 "period_days",
#         #                 "star_1_mass",
#         #                 termination_flag="termination_flag_1",
#         #                 grid_3D=True,
#         #                 slice_3D_var_str="star_2_mass",
#         #                 slice_3D_var_range=(2.5, 3.0),
#         #                 slice_at_RLO=True,
#         #                 **{
#         #                 'show_fig': True
#         #                 })
#         #     assert show_patch.called
#
#     def test_extra_grid_plotting(self):
#         grid = PSyGrid(PATH_TO_GRID)
#         with patch('matplotlib.pyplot.show') as show_patch:
#             grid.plot2D("star_1_mass",
#                         "period_days",
#                         "star_1_mass",
#                         termination_flag="termination_flag_1",
#                         grid_3D=True,
#                         slice_3D_var_str="star_2_mass",
#                         slice_3D_var_range=(2.5, 3.0),
#                         extra_grid=grid,
#                         **{'show_fig': True})
#         with patch('matplotlib.pyplot.show') as show_patch:
#             grid.plot2D("star_1_mass",
#                         "period_days",
#                         "star_1_mass",
#                         termination_flag="termination_flag_1",
#                         grid_3D=True,
#                         slice_3D_var_str="star_2_mass",
#                         slice_3D_var_range=(2.5, 3.0),
#                         extra_grid=grid,
#                         **{'show_fig': True})
#             assert show_patch.called
#
#
# if __name__ == "__main__":
#     unittest.main()
