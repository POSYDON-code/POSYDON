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

# Theoretical WD mass limits of spherical WDs depend on composition
## classical Chandrasekhar limit, see Eq. (63) in Chandrasekhar (1935):
## M3 = 5.728 / µ^2 Msol = 1.432 Msol for µ = 2
## from Eq. (10.3) in Pols (2009) or Eq. (43) in Bressan & Shepherd (2024)
##  following Kippenhahn & Weigert (1990): MCh = 1.459 * (2/µe)^2 Msol
## from Eq. (57) in Rueda & Ruffini (2013): 1.44 Msol
## from Table 1 in Carvalho, Marinho & Malheiro (2017): M(Newt.) = 1.4546 Msol,
##  M(SR) = 1.4358 Msol, M(GR) = 1.4154 Msol, M(non-rel. Newt.) = 1.4564 Msol
## review from Saumon et al. (2022): canonical Chandrasekhar mass of 1.4 Msol
# Roting and/or magnetic WDs can have a mass above the Chandrasekhar limit
## from Das & Mukhopadhyay (2012) about magnetic WDs: up to 2.3-2.6 Msol
# There are direct and indirect observational constraints
## most massive WD from observations in Vennes et al. (1997): 1.41+-0.04 Msol
## massive WDs in solar neighborhood in Kilic et al. (2021): 1.351+-0.006 Msol
## overluminous Type Ia supernova require WDs up to 2.8 Msol, see Das &
##  Mukhopadhyay (2013)
# We use a conservative value for non-rotating WDs of 1.46 Msol
STATE_WD_STARMASS_UPPER_LIMIT = 1.46 # maximum mass of a white dwarf (in Msol)

# The limits for neutron stars depend on the underlying equation of state
## secure upper bound on the neutron star mass in Kalogera & Baym (1996):
##  mmax = 2.9 Msol
## maximum non-rotating neutron star masses in Fryer et al. (2015):
##  mmax = 2.3-2.4 Msol
## from Table 4 in Alsing, Silva & Berti (2017): mmax = 1.96-2.79 Msol
# Most extreme masses of neutron stars observed in NS+WD binaries:
## largest mass from Antoniadis et al. (2013): MNS = 2.01+-0.04 Msol
## lowest mass from Fonseca et al. (2016): MNS = 1.18+0.10-0.09 Msol
# Most extreme masses of neutron stars observed from black widows:
## largest millisecond pulsar from Freire et al. (2008): MNS = 2.74+-0.21 Msol
## from van Kerkwijk, Breton & Kulkarni (2011): MNS = 2.40+-0.12 Msol
## from Romani et al. (2012): MNS = 2.68+-0.14 Msol
## from Romani et al. (2015): MNS = 2.63+0.3-0.2 Msol
# Most massive NS in X-ray binaries:
## from Falanga et al. (2015): MNS = 2.12+-0.16 Msol
# Combining different EM observations from Alsing, Silva & Berti (2017):
#  MNS < 2.6 Msol
# Most massive NS in GW mergers:
## GW190814 from Abbott et al. (2021): m2 = 2.59+0.08-0.09 Msol
## GW200210_092254 from Abbott et al. (2023): m2 = 2.83+0.47-0.42 Msol
# We use a conservative value of >1.1 Msol and the canonical value of <2.5 Msol
STATE_NS_STARMASS_LOWER_LIMIT = 1.1  # minimum mass of a neutron star (in Msol)
STATE_NS_STARMASS_UPPER_LIMIT = 2.5  # maximum mass of a neutron star (in Msol)

# During a collapse escaping neutrinos reduce the mass ending up in the compact
# object
NEUTRINO_MASS_LOSS_UPPER_LIMIT = 0.5 # maximum loss in neutrinos (in Msol)

