import unittest
from posydon.utils.common_functions import beaming
from posydon.binary_evol.singlestar import SingleStar
from posydon.binary_evol.binarystar import BinaryStar

class TestEddingtonBeaming(unittest.TestCase):
    def test_beaming_1(self):
        # One of the stars should be a compact object (NS or BH) and it should go in
        STAR1_PROP = {"mass": 1.4, "state": "NS"}
        STAR2_PROP = {"surface_h1": 0.7}
        BINARY_PROP = {
            "lg_mtransfer_rate": -2
        }
        star1 = SingleStar(**STAR1_PROP)
        star2 = SingleStar(**STAR2_PROP)
        binary = BinaryStar(star1, star2, **BINARY_PROP)

        self.assertAlmostEqual(
            beaming(binary)[0],
            874.55367846,
            places=5,
            msg="Should be 874.5537."
        )

    def test_beaming_2(self):
        # One of the stars should be a compact object (NS or BH) and it should go in
        # the star1 argument position
        STAR1_PROP = {"mass": 1.4, "state": "NS"}
        STAR2_PROP = {"surface_h1": 0.7}
        BINARY_PROP = {
            "lg_mtransfer_rate": -2
        }
        star1 = SingleStar(**STAR1_PROP)
        star2 = SingleStar(**STAR2_PROP)
        binary = BinaryStar(star1, star2, **BINARY_PROP)

        self.assertAlmostEqual(
            beaming(binary)[1],
            3.524620449e-10,
            places=5,
            msg="Should be 3.524620e-10."
        )

    def test_beaming_3(self):
        # One of the stars should be a compact object (NS or BH) and it should go in
        STAR1_PROP = {"mass": 10.0, "state": "BH"}
        STAR2_PROP = {"surface_h1": 0.9}
        BINARY_PROP = {
            "lg_mtransfer_rate": -2
        }
        star1 = SingleStar(**STAR1_PROP)
        star2 = SingleStar(**STAR2_PROP)
        binary = BinaryStar(star1, star2, **BINARY_PROP)

        self.assertAlmostEqual(
            beaming(binary)[0],
            37.466319,
            places=5,
            msg="Should be 37.4663193."
        )

if __name__ == "__main__":
    unittest.main()
