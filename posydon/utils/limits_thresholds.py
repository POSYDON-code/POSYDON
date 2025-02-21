"""Constant limits and thresholds.
"""


__authors__ = [
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>"
]

import numpy as np


# THRESHOLD USED FOR INFERRING MASS-TRANSFER, AND RLO

# This needs to be aligned with run_binary_extras.f
# old: RL_RELATIVE_OVERFLOW_THRESHOLD = -0.05
RL_RELATIVE_OVERFLOW_THRESHOLD = 0.00  # relative overflow threshold
# This needs to be aligned with run_binary_extras.f
# old: LG_MTRANSFER_RATE_THRESHOLD = -12
LG_MTRANSFER_RATE_THRESHOLD = -10      # mass-transfer rate threshold (in Lsol)
# Numerical limit from the number of datapoints the automatically detected
# boundary for initial RLO: the value should be ≳3 (geometrical arguments) and
# <<1000 (well below grid size).
MIN_COUNT_INITIAL_RLO_BOUNDARY = 20


# CONSTANTS RELATED TO INFERRING STAR STATES

THRESHOLD_CENTRAL_ABUNDANCE = 0.01   # central abundance for flagging depletion
THRESHOLD_CENTRAL_ABUNDANCE_LOOSE_C = 0.1  # central C abundance allow cut
THRESHOLD_HE_NAKED_ABUNDANCE = 0.01  # for surface abundance for stripped_He
THRESHOLD_NUCLEAR_LUMINOSITY = 0.97  # for element fraction in nuclear burning
# relative burning threshold with respect to nuclear luminosity
REL_LOG10_BURNING_THRESHOLD = np.log10(1.0 - THRESHOLD_NUCLEAR_LUMINOSITY)
LOG10_BURNING_THRESHOLD = -10.0      # burning luminosity threshold (in Lsol)


# COMPACT OBJECT LIMITS

# Roting WDs can have a mass above the classical Chandrasekhar limit of 1.37 to
# 1.4 Msol depending on the composition
STATE_WD_STARMASS_UPPER_LIMIT = 1.48 # maximum mass of a white dwarf (in Msol)
# The limits for neutron stars depend on the underlying equation of state
STATE_NS_STARMASS_LOWER_LIMIT = 1.0  # minimum mass of a neutron star (in Msol)
STATE_NS_STARMASS_UPPER_LIMIT = 2.5  # maximum mass of a neutron star (in Msol)
# During a collapse escaping neutrinos reduce the mass ending up in the compact
# object
NEUTRINO_MASS_LOSS_UPPER_LIMIT = 0.5 # maximum loss in neutrinos (in Msol)


