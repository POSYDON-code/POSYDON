import unittest
from posydon.utils.common_functions import eddington_limit
from posydon.binary_evol.singlestar import SingleStar
from posydon.binary_evol.binarystar import BinaryStar


class TestEddingtonLimit(unittest.TestCase):
    def test_eddington_1(self):
        # One of the stars should be a compact object (NS or BH) and it should go in
        # the star1 argument position
        STAR1_PROP = {"mass": 2.831, "state": "NS"}
        STAR2_PROP = {"surface_h1": 0.715}

        star1 = SingleStar(**STAR1_PROP)
        star2 = SingleStar(**STAR2_PROP)
        binary = BinaryStar(star1, star2)

        self.assertAlmostEqual(
            eddington_limit(binary, idx=-1)[0],
            1.274160439238623e-07,
            places=5,
            msg="Should be 1.27416e-07.",
        )

    def test_eddington_2(self):
        # One of the stars should be a compact object (NS or BH) and it should go in
        # the star1 argument position
        STAR1_PROP = {"mass": 2.831, "state": "BH"}
        STAR2_PROP = {"surface_h1": 0.715}

        star1 = SingleStar(**STAR1_PROP)
        star2 = SingleStar(**STAR2_PROP)
        binary = BinaryStar(star1, star2)

        self.assertAlmostEqual(
            eddington_limit(binary, idx=-1)[1],
            0.0571909584179,
            places=5,
            msg="Should be 0.05719.",
        )

    def test_eddington_3(self):
        # One of the stars should be a compact object (NS or BH) and it should go in
        # the star1 argument position
        STAR1_PROP = {"mass": 1.5, "state": "NS"}
        STAR2_PROP = {"surface_h1": 0.85}

        star1 = SingleStar(**STAR1_PROP)
        star2 = SingleStar(**STAR2_PROP)
        binary = BinaryStar(star1, star2)

        self.assertAlmostEqual(
            eddington_limit(binary, idx=-1)[0],
            1.7768657384e-08,
            places=5,
            msg="Should be 1.77687e-08.",
        )


if __name__ == "__main__":
    unittest.main()
