import unittest
from posydon.utils.common_functions import is_number, rzams


class TestIsNumber(unittest.TestCase):
    def test_is_number_1(self):

        self.assertTrue(
            is_number(17.5),
            msg="Should be True.",
        )

    def test_is_number_2(self):

        self.assertFalse(
            is_number("not a float"),
            msg="Should be False.",
        )


class TestRZAMS(unittest.TestCase):
    def test_rzams_1(self):

        self.assertAlmostEqual(
            rzams(1.0),
            0.88824945,
            places=5,
            msg="Should be 0.8882 Rsun.",
        )

    def test_rzams_2(self):

        self.assertAlmostEqual(
            rzams(10.0),
            3.941905620,
            places=5,
            msg="Should be 3.9419 Rsun.",
        )

    def test_rzams_3(self):

        self.assertAlmostEqual(
            rzams(100.0),
            17.22611585,
            places=5,
            msg="Should be 17.2261 Rsun.",
        )


if __name__ == "__main__":
    unittest.main()
