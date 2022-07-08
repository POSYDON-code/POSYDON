import unittest
import posydon.utils.constants as const
from posydon.utils.common_functions import (
    orbital_separation_from_period,
    orbital_period_from_separation,
)


class TestKeplerThirdLaw(unittest.TestCase):
    def test_sep_from_per_1(self):

        self.assertAlmostEqual(
            orbital_separation_from_period(10.0, 1.0, 1.0),
            24.60349357746796,
            places=5,
            msg="Should be 24.603493 Rsun.",
        )

    def test_sep_from_per_2(self):

        self.assertAlmostEqual(
            orbital_separation_from_period(5.0, 3.0, 1.0),
            19.527805793482745,
            places=5,
            msg="Should be 19.52780579 Rsun.",
        )

    def test_sep_from_per_3(self):

        self.assertAlmostEqual(
            orbital_separation_from_period(10.0, 0.0001, 1.0),
            19.52845669864617,
            places=5,
            msg="Should be 19.5284567 Rsun.",
        )

    def test_per_from_sep_1(self):

        self.assertAlmostEqual(
            orbital_period_from_separation(10.0, 1.0, 1.0),
            2.591606452663911,
            places=5,
            msg="Should be 2.5916 days.",
        )

    def test_per_from_sep_2(self):

        self.assertAlmostEqual(
            orbital_period_from_separation(5.0, 3.0, 1.0),
            0.6479016131659777,
            places=5,
            msg="Should be 0.6479 days.",
        )

    def test_per_from_sep_3(self):

        self.assertAlmostEqual(
            orbital_period_from_separation(10.0, 0.001, 1.0),
            3.6632538244566186,
            places=5,
            msg="Should be 3.6633 days.",
        )


if __name__ == "__main__":
    unittest.main()
