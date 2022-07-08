import unittest
from posydon.utils.common_functions import inspiral_timescale_from_separation


class TestMerger_sep(unittest.TestCase):
    def test_merger_sep_1(self):

        self.assertAlmostEqual(
            inspiral_timescale_from_separation(15.0, 30.0, 14.97, 0.0),
            372.668437,
            places=5,
            msg="Should be 372.6684Myr.",
        )

    def test_merger_sep_2(self):

        self.assertAlmostEqual(
            inspiral_timescale_from_separation(15.0, 30.0, 14.97, 0.2),
            321.23584297,
            places=5,
            msg="Should be 321.2358Myr.",
        )

    def test_merger_sep_3(self):

        self.assertAlmostEqual(
            inspiral_timescale_from_separation(15.0, 30.0, 14.97, 0.5),
            133.3401165,
            places=5,
            msg="Should be 133.3401Myr.",
        )

    def test_merger_sep_4(self):
        # Tests the continuity of the function approaching 0 ecc.
        self.assertAlmostEqual(
            inspiral_timescale_from_separation(15.0, 30.0, 14.97, 0.00001),
            372.668437,
            places=5,
            msg="Should be 372.668Myr.",
        )


if __name__ == "__main__":
    unittest.main()
