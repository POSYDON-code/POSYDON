import unittest
import os
from posydon.utils.common_functions import PATH_TO_POSYDON
from posydon.grids.psygrid import PSyGrid
import posydon.utils.constants as const
from posydon.binary_evol.singlestar import SingleStar
from posydon.binary_evol.SN.profile_collapse import do_core_collapse_BH, get_initial_BH_properties, compute_isco_properties

PATH_TO_GRID = os.path.join(
    PATH_TO_POSYDON, "posydon/tests/data/POSYDON-UNIT-TESTS/"
    "visualization/grid_unit_test_plot.h5")

if not os.path.isfile(PATH_TO_GRID):
    print(PATH_TO_GRID)
    raise ValueError("Test grid for unit testing was not found!")

# constants in CGS
G = const.standard_cgrav
c = const.clight
Mo = const.Msun


class TestProfileCollapse(unittest.TestCase):
    def test_r_isco(self):
        m_BH = 1. * Mo
        self.assertAlmostEqual(compute_isco_properties(0., m_BH)[0] /
                               (G * m_BH / c**2),
                               6.0,
                               places=5)
        self.assertAlmostEqual(compute_isco_properties(0.999, m_BH)[0] /
                               (G * m_BH / c**2),
                               1.1817646130335708,
                               places=5)

    def test_j_isco(self):
        m_BH = 1. * Mo
        self.assertAlmostEqual(compute_isco_properties(0, m_BH)[1] /
                               (G * m_BH / c),
                               3.464101615137754,
                               places=5)
        self.assertAlmostEqual(compute_isco_properties(0.999, m_BH)[1] /
                               (G * m_BH / c),
                               1.3418378380509774,
                               places=5)

    def test_radiation_efficiency(self):
        m_BH = 1. * Mo
        self.assertAlmostEqual((1 - compute_isco_properties(0., m_BH)[2]),
                               0.057190958417936644,
                               places=5)
        self.assertAlmostEqual((1 - compute_isco_properties(0.999, m_BH)[2]),
                               0.3397940734762088,
                               places=5)

    def test_low_spinning_He_star(self):
        grid = PSyGrid(PATH_TO_GRID)
        i = 42
        star = SingleStar(**{'profile': grid[i].final_profile1})
        m_rembar = grid[i].final_values['star_1_mass']
        mass_direct_collapse = 3.  # Msun
        delta_M = 0.5  # Msun
        results = do_core_collapse_BH(star, m_rembar, mass_direct_collapse,
                                      delta_M)
        self.assertAlmostEqual(results[0], 13.365071929231409, places=5)
        self.assertAlmostEqual(results[1], 8.98074719361575e-09, places=5)

    def test_midly_spinning_He_star(self):
        grid = PSyGrid(PATH_TO_GRID)
        i = 13
        star = SingleStar(**{'profile': grid[i].final_profile1})
        m_rembar = grid[i].final_values['star_1_mass']
        mass_direct_collapse = 3.  # Msun
        delta_M = 0.5  # Msun
        results = do_core_collapse_BH(star, m_rembar, mass_direct_collapse,
                                      delta_M)
        self.assertAlmostEqual(results[0], 5.60832288900688, places=5)
        self.assertAlmostEqual(results[1], 0.42583967572001924, places=5)

    def test_rapidly_spinning_He_star(self):
        grid = PSyGrid(PATH_TO_GRID)
        i = 6
        star = SingleStar(**{'profile': grid[i].final_profile1})
        m_rembar = grid[i].final_values['star_1_mass']
        mass_direct_collapse = 3.  # Msun
        delta_M = 0.5  # Msun
        results = do_core_collapse_BH(star, m_rembar, mass_direct_collapse,
                                      delta_M)
        self.assertAlmostEqual(results[0], 38.50844589130613, places=5)
        self.assertAlmostEqual(results[1], 0.9835226614001595, places=5)


if __name__ == "__main__":
    unittest.main()
