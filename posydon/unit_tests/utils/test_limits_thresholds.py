"""Unit tests of posydon/utils/limits_thresholds.py
"""

__authors__ = [
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>"
]

# import the module which will be tested
import posydon.utils.limits_thresholds as totest

# import other needed code for the tests, which is not already imported in the
# module you like to test


# define test classes collecting several test functions
class TestElements:
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = ['LG_MTRANSFER_RATE_THRESHOLD', 'LOG10_BURNING_THRESHOLD',\
                    'NEUTRINO_MASS_LOSS_UPPER_LIMIT',\
                    'REL_LOG10_BURNING_THRESHOLD',\
                    'RL_RELATIVE_OVERFLOW_THRESHOLD',\
                    'STATE_NS_STARMASS_UPPER_LIMIT',\
                    'THRESHOLD_CENTRAL_ABUNDANCE',\
                    'THRESHOLD_CENTRAL_ABUNDANCE_LOOSE_C',\
                    'THRESHOLD_HE_NAKED_ABUNDANCE',\
                    'THRESHOLD_NUCLEAR_LUMINOSITY', '__authors__',\
                    '__builtins__', '__cached__', '__doc__', '__file__',\
                    '__loader__', '__name__', '__package__', '__spec__', 'np']
        assert dir(totest) == elements, "There might be added or removed "+\
               "objects without an update on the unit test."

    def test_instance_RL_RELATIVE_OVERFLOW_THRESHOLD(self):
        assert isinstance(totest.RL_RELATIVE_OVERFLOW_THRESHOLD, (float,int)),\
               "RL_RELATIVE_OVERFLOW_THRESHOLD is of type: "+\
               str(type(totest.RL_RELATIVE_OVERFLOW_THRESHOLD))

    def test_instance_LG_MTRANSFER_RATE_THRESHOLD(self):
        assert isinstance(totest.LG_MTRANSFER_RATE_THRESHOLD, (float,int)),\
               "LG_MTRANSFER_RATE_THRESHOLD is of type: "+\
               str(type(totest.LG_MTRANSFER_RATE_THRESHOLD))

    def test_instance_THRESHOLD_CENTRAL_ABUNDANCE(self):
        assert isinstance(totest.THRESHOLD_CENTRAL_ABUNDANCE, (float,int)),\
               "THRESHOLD_CENTRAL_ABUNDANCE is of type: "+\
               str(type(totest.THRESHOLD_CENTRAL_ABUNDANCE))

    def test_instance_THRESHOLD_CENTRAL_ABUNDANCE_LOOSE_C(self):
        assert isinstance(totest.THRESHOLD_CENTRAL_ABUNDANCE_LOOSE_C,\
               (float,int)), "THRESHOLD_CENTRAL_ABUNDANCE_LOOSE_C is of "\
               "type: "+str(type(totest.THRESHOLD_CENTRAL_ABUNDANCE_LOOSE_C))

    def test_instance_THRESHOLD_HE_NAKED_ABUNDANCE(self):
        assert isinstance(totest.THRESHOLD_HE_NAKED_ABUNDANCE, (float,int)),\
               "THRESHOLD_HE_NAKED_ABUNDANCE is of type: "+\
               str(type(totest.THRESHOLD_HE_NAKED_ABUNDANCE))

    def test_instance_THRESHOLD_NUCLEAR_LUMINOSITY(self):
        assert isinstance(totest.THRESHOLD_NUCLEAR_LUMINOSITY, (float,int)),\
               "THRESHOLD_NUCLEAR_LUMINOSITY is of type: "+\
               str(type(totest.THRESHOLD_NUCLEAR_LUMINOSITY))

    def test_instance_REL_LOG10_BURNING_THRESHOLD(self):
        assert isinstance(totest.REL_LOG10_BURNING_THRESHOLD, (float,int)),\
               "REL_LOG10_BURNING_THRESHOLD is of type: "+\
               str(type(totest.REL_LOG10_BURNING_THRESHOLD))

    def test_instance_LOG10_BURNING_THRESHOLD(self):
        assert isinstance(totest.LOG10_BURNING_THRESHOLD, (float,int)),\
               "LOG10_BURNING_THRESHOLD is of type: "+\
               str(type(totest.LOG10_BURNING_THRESHOLD))

    def test_instance_STATE_NS_STARMASS_UPPER_LIMIT(self):
        assert isinstance(totest.STATE_NS_STARMASS_UPPER_LIMIT, (float,int)),\
               "STATE_NS_STARMASS_UPPER_LIMIT is of type: "+\
               str(type(totest.STATE_NS_STARMASS_UPPER_LIMIT))

    def test_instance_NEUTRINO_MASS_LOSS_UPPER_LIMIT(self):
        assert isinstance(totest.NEUTRINO_MASS_LOSS_UPPER_LIMIT, (float,int)),\
               "NEUTRINO_MASS_LOSS_UPPER_LIMIT is of type: "+\
               str(type(totest.NEUTRINO_MASS_LOSS_UPPER_LIMIT))


class TestLimits:
    # check for validity ranges
#    def test_limits_RL_RELATIVE_OVERFLOW_THRESHOLD(self):
        # has no limits

#    def test_limits_LG_MTRANSFER_RATE_THRESHOLD(self):
        # has no limits

    def test_limits_THRESHOLD_CENTRAL_ABUNDANCE(self):
        # an abundance should be in [0,1]
        assert totest.THRESHOLD_CENTRAL_ABUNDANCE >= 0.0,\
               "THRESHOLD_CENTRAL_ABUNDANCE should be in the range [0,1], "+\
               "but it is below"
        assert totest.THRESHOLD_CENTRAL_ABUNDANCE <= 1.0,\
               "THRESHOLD_CENTRAL_ABUNDANCE should be in the range [0,1], "+\
               "but it is above"

    def test_limits_THRESHOLD_CENTRAL_ABUNDANCE_LOOSE_C(self):
        # an abundance should be in [0,1]
        assert totest.THRESHOLD_CENTRAL_ABUNDANCE_LOOSE_C >= 0.0,\
               "THRESHOLD_CENTRAL_ABUNDANCE_LOOSE_C should be in the range "+\
               "[0,1], but it is below"
        assert totest.THRESHOLD_CENTRAL_ABUNDANCE_LOOSE_C <= 1.0,\
               "THRESHOLD_CENTRAL_ABUNDANCE_LOOSE_C should be in the range "+\
               "[0,1], but it is above"
        # it should be limited by THRESHOLD_CENTRAL_ABUNDANCE
        assert totest.THRESHOLD_CENTRAL_ABUNDANCE_LOOSE_C >=\
               totest.THRESHOLD_CENTRAL_ABUNDANCE,\
               "the loose condition should be less strict"

    def test_limits_THRESHOLD_HE_NAKED_ABUNDANCE(self):
        # an abundance should be in [0,1]
        assert totest.THRESHOLD_HE_NAKED_ABUNDANCE >= 0.0,\
               "THRESHOLD_HE_NAKED_ABUNDANCE should be in the range [0,1], "+\
               "but it is below"
        assert totest.THRESHOLD_HE_NAKED_ABUNDANCE <= 1.0,\
               "THRESHOLD_HE_NAKED_ABUNDANCE should be in the range [0,1], "+\
               "but it is above"

    def test_limits_THRESHOLD_NUCLEAR_LUMINOSITY(self):
        # an fraction should be in [0,1]
        assert totest.THRESHOLD_NUCLEAR_LUMINOSITY >= 0.0,\
               "THRESHOLD_NUCLEAR_LUMINOSITY should be in the range [0,1], "+\
               "but it is below"
        assert totest.THRESHOLD_NUCLEAR_LUMINOSITY <= 1.0,\
               "THRESHOLD_NUCLEAR_LUMINOSITY should be in the range [0,1], "+\
               "but it is above"

    def test_limits_REL_LOG10_BURNING_THRESHOLD(self):
        # the log of a fraction should be <=0
        assert totest.REL_LOG10_BURNING_THRESHOLD <= 0.0, "a fraction should"+\
               " be in the range [0,1], hence the log of it can't be positive"

#    def test_limits_LOG10_BURNING_THRESHOLD(self):
        # has no limits

    def test_limits_STATE_NS_STARMASS_UPPER_LIMIT(self):
        # a mass should be >0
        assert totest.STATE_NS_STARMASS_UPPER_LIMIT > 0.0,\
               "a mass has to be positve"

    def test_limits_NEUTRINO_MASS_LOSS_UPPER_LIMIT(self):
        # a mass limit should be >=0
        assert totest.NEUTRINO_MASS_LOSS_UPPER_LIMIT >= 0.0,\
               "there shouldn't be a mass gain"
