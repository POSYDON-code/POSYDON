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


# CONSTANTS RELATED TO INFERRING STAR STATES

THRESHOLD_CENTRAL_ABUNDANCE = 0.01   # central abundance for flagging depletion
THRESHOLD_CENTRAL_ABUNDANCE_LOOSE_C = 0.1  # central C abundance allow cut
THRESHOLD_HE_NAKED_ABUNDANCE = 0.01  # for surface abundance for stripped_He
THRESHOLD_NUCLEAR_LUMINOSITY = 0.97  # for element fraction in nuclear burning
# relative burning threshold with respect to nuclear luminosity
REL_LOG10_BURNING_THRESHOLD = np.log10(1.0 - THRESHOLD_NUCLEAR_LUMINOSITY)
LOG10_BURNING_THRESHOLD = -10.0      # burning luminosity threshold (in Lsol)


# COMPACT OBJECT LIMITS

STATE_NS_STARMASS_UPPER_LIMIT = 2.5   # maximum mass of a neutron star (in Msol)
NEUTRINO_MASS_LOSS_UPPER_LIMIT = 0.5  # maximum loss in neutrinos (in Msol)


