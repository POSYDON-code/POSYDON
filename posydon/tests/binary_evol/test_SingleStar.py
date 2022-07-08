import unittest
import os
import numpy as np
from posydon.utils.common_functions import PATH_TO_POSYDON
from posydon.binary_evol.singlestar import SingleStar
from posydon.grids.psygrid import PSyGrid

PATH_TO_GRID = os.path.join(
    PATH_TO_POSYDON,
    "posydon/tests/data/POSYDON-UNIT-TESTS/"
    "visualization/grid_unit_test_plot.h5"
)

if not os.path.exists(PATH_TO_GRID):
    raise ValueError("Test grid for unit testing was not found!")


class TestSingleStar(unittest.TestCase):
    def test_SingleStar_initialisation(self):
        # load an example grid: compact object + He-star
        grid = PSyGrid(PATH_TO_GRID)

        # initialise a star with the properties of run i=42
        i = 42

        # all STARPROPERTIES
        kwargs = {
            'state': 'stripped_He_Central_C_depletion',
            'metallicity': grid[i].initial_values['Z'],
            'mass': grid[i].history1['star_mass'][-1],
            'log_R': np.nan,
            'log_L': grid[i].history1['log_L'][-1],
            'lg_mdot': np.nan,
            'lg_system_mdot': np.nan,
            'lg_wind_mdot': np.nan,
            'he_core_mass': grid[i].history1['he_core_mass'][-1],
            'he_core_radius': np.nan,
            'c_core_radius': grid[i].history1['he_core_mass'][-1],
            'o_core_mass': np.nan,
            'o_core_radius': np.nan,
            'center_h1': grid[i].history1['center_h1'][-1],
            'center_he4': grid[i].history1['center_he4'][-1],
            'center_c12': grid[i].history1['center_c12'][-1],
            'center_n14': np.nan,
            'center_o16': np.nan,
            'surface_h1': grid[i].history1['surface_h1'][-1],
            'surface_he4': np.nan,
            'surface_c12': np.nan,
            'surface_n14': np.nan,
            'surface_o16': np.nan,
            'log_LH': grid[i].history1['log_LH'][-1],
            'log_LHe': grid[i].history1['log_LHe'][-1],
            'log_LZ': grid[i].history1['log_LZ'][-1],
            'log_Lnuc': grid[i].history1['log_Lnuc'][-1],
            'c12_c12': grid[i].history1['c12_c12'][-1],
            'surf_avg_omega_div_omega_crit': np.nan,
            'total_moment_of_inertia': np.nan,
            'log_total_angular_momentum': np.nan,
            'spin': np.nan,
            'profile': grid[i].final_profile1
        }

        star = SingleStar(**kwargs)

        # check that the above kwars have a history
        for item in kwargs.keys():
            self.assertIsInstance(getattr(star, item + '_history'), list)


if __name__ == '__main__':
    unittest.main()
