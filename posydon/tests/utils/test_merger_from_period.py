import unittest
from posydon.utils.common_functions import inspiral_timescale_from_orbital_period


class TestMerger_porb(unittest.TestCase):
    def test_merger_porb_1(self):

        self.assertAlmostEqual(
            inspiral_timescale_from_orbital_period(15.0, 30.0, 1.0, 0.0),
            372.10532488,
            places=5,
            msg="Should be 372.1053Myr.",
        )

    def test_merger_porb_2(self):

        self.assertAlmostEqual(
            inspiral_timescale_from_orbital_period(15.0, 30.0, 1.0, 0.2),
            320.75044687,
            places=5,
            msg="Should be 320.7504Myr.",
        )

    def test_merger_porb_3(self):

        self.assertAlmostEqual(
            inspiral_timescale_from_orbital_period(15.0, 30.0, 1.0, 0.5),
            133.138635978858,
            places=5,
            msg="Should be 133.1386Myr.",
        )


if __name__ == "__main__":
    unittest.main()
