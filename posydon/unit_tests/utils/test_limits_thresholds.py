"""Unit tests of posydon/utils/limits_thresholds.py
"""

__authors__ = [
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>"
]

# import the unittest module and the module which will be tested
import unittest
import posydon.utils.limits_thresholds as totest

# import other needed code for the tests
# not needed

# define test classes
class TestElements(unittest.TestCase):
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = ['LG_MTRANSFER_RATE_THRESHOLD', 'LOG10_BURNING_THRESHOLD',
                    'NEUTRINO_MASS_LOSS_UPPER_LIMIT',
                    'REL_LOG10_BURNING_THRESHOLD',
                    'RL_RELATIVE_OVERFLOW_THRESHOLD',
                    'STATE_NS_STARMASS_UPPER_LIMIT',
                    'THRESHOLD_CENTRAL_ABUNDANCE',
                    'THRESHOLD_CENTRAL_ABUNDANCE_LOOSE_C',
                    'THRESHOLD_HE_NAKED_ABUNDANCE',
                    'THRESHOLD_NUCLEAR_LUMINOSITY', '__authors__',
                    '__builtins__', '__cached__', '__doc__', '__file__',
                    '__loader__', '__name__', '__package__', '__spec__', 'np']
        self.assertListEqual(dir(totest), elements,
                             msg="There might be added or removed objects "
                                 "without an update on the unit test.")

    def test_instance_RL_RELATIVE_OVERFLOW_THRESHOLD(self):
        self.assertIsInstance(totest.RL_RELATIVE_OVERFLOW_THRESHOLD,
                              (float, int))

    def test_instance_LG_MTRANSFER_RATE_THRESHOLD(self):
        self.assertIsInstance(totest.LG_MTRANSFER_RATE_THRESHOLD,
                              (float, int))

    def test_instance_THRESHOLD_CENTRAL_ABUNDANCE(self):
        self.assertIsInstance(totest.THRESHOLD_CENTRAL_ABUNDANCE,
                              (float, int))

    def test_instance_THRESHOLD_CENTRAL_ABUNDANCE_LOOSE_C(self):
        self.assertIsInstance(totest.THRESHOLD_CENTRAL_ABUNDANCE_LOOSE_C,
                              (float, int))

    def test_instance_THRESHOLD_HE_NAKED_ABUNDANCE(self):
        self.assertIsInstance(totest.THRESHOLD_HE_NAKED_ABUNDANCE,
                              (float, int))

    def test_instance_THRESHOLD_NUCLEAR_LUMINOSITY(self):
        self.assertIsInstance(totest.THRESHOLD_NUCLEAR_LUMINOSITY,
                              (float, int))

    def test_instance_REL_LOG10_BURNING_THRESHOLD(self):
        self.assertIsInstance(totest.REL_LOG10_BURNING_THRESHOLD,
                              (float, int))

    def test_instance_LOG10_BURNING_THRESHOLD(self):
        self.assertIsInstance(totest.LOG10_BURNING_THRESHOLD,
                              (float, int))

    def test_instance_STATE_NS_STARMASS_UPPER_LIMIT(self):
        self.assertIsInstance(totest.STATE_NS_STARMASS_UPPER_LIMIT,
                              (float, int))

    def test_instance_NEUTRINO_MASS_LOSS_UPPER_LIMIT(self):
        self.assertIsInstance(totest.NEUTRINO_MASS_LOSS_UPPER_LIMIT,
                              (float, int))


class TestLimits(unittest.TestCase):
    # check for validity ranges
#    def test_limits_RL_RELATIVE_OVERFLOW_THRESHOLD(self):
        # has no limits

#    def test_limits_LG_MTRANSFER_RATE_THRESHOLD(self):
        # has no limits

    def test_limits_THRESHOLD_CENTRAL_ABUNDANCE(self):
        # an abundance should be in [0,1]
        self.assertGreaterEqual(totest.THRESHOLD_CENTRAL_ABUNDANCE, 0.0)
        self.assertLessEqual(totest.THRESHOLD_CENTRAL_ABUNDANCE, 1.0)

    def test_limits_THRESHOLD_CENTRAL_ABUNDANCE_LOOSE_C(self):
        # an abundance should be in [0,1]
        self.assertGreaterEqual(totest.THRESHOLD_CENTRAL_ABUNDANCE_LOOSE_C, 0.0)
        self.assertLessEqual(totest.THRESHOLD_CENTRAL_ABUNDANCE_LOOSE_C, 1.0)
        # it should be limited by THRESHOLD_CENTRAL_ABUNDANCE
        self.assertGreaterEqual(totest.THRESHOLD_CENTRAL_ABUNDANCE_LOOSE_C, totest.THRESHOLD_CENTRAL_ABUNDANCE)

    def test_limits_THRESHOLD_HE_NAKED_ABUNDANCE(self):
        # an abundance should be in [0,1]
        self.assertGreaterEqual(totest.THRESHOLD_HE_NAKED_ABUNDANCE, 0.0)
        self.assertLessEqual(totest.THRESHOLD_HE_NAKED_ABUNDANCE, 1.0)

    def test_limits_THRESHOLD_NUCLEAR_LUMINOSITY(self):
        # an fraction should be in [0,1]
        self.assertGreaterEqual(totest.THRESHOLD_NUCLEAR_LUMINOSITY, 0.0)
        self.assertLessEqual(totest.THRESHOLD_NUCLEAR_LUMINOSITY, 1.0)

    def test_limits_REL_LOG10_BURNING_THRESHOLD(self):
        # the log of a fraction should be <0
        self.assertLessEqual(totest.REL_LOG10_BURNING_THRESHOLD, 0.0)

#    def test_limits_LOG10_BURNING_THRESHOLD(self):
        # has no limits

    def test_limits_STATE_NS_STARMASS_UPPER_LIMIT(self):
        # a mass should be >0
        self.assertGreater(totest.STATE_NS_STARMASS_UPPER_LIMIT, 0.0)

    def test_limits_NEUTRINO_MASS_LOSS_UPPER_LIMIT(self):
        # a mass limit should be >=0
        self.assertGreaterEqual(totest.NEUTRINO_MASS_LOSS_UPPER_LIMIT, 0.0)


if __name__ == "__main__":
    unittest.main()
