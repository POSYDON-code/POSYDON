"""Mathematical and physical constants (in cgs).

Using the 2010 CODATA recommended values of the physical constants:

    https://physics.nist.gov/cuu/Constants/Preprints/lsa2010.pdf

"""


__authors__ = [
    "Nam Tran <tranhn03@gmail.com>",
    "Emmanouil Zapartas <ezapartas@gmail.com>",
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Konstantinos Kovlakas <Konstantinos.Kovlakas@unige.ch>",
]


# MATHEMATICAL CONSTANTS

pi = 3.1415926535897932384626433832795028841971693993751
a2rad = pi / 180.0                  # angle to radians
rad2a = 180.0 / pi                  # radians to angle


# PHYSICAL CONSTANTS

standard_cgrav = 6.67428e-8         # gravitational constant (g^-1 cm^3 s^-2)
planck_h = 6.62606896e-27           # Planck's constant (erg s)
hbar = planck_h / (2*pi)
qe = 4.80320440e-10                 # electron charge (esu = g cm^3 s^-2)^1/2)
avo = 6.02214129e23                 # Avogadro's constant (mole^-1)
clight = 2.99792458e10              # speed of light in vacuum (cm s^-1)
kerg = 1.3806504e-16                # Boltzmann's constant (erg K^-1)
boltzm = kerg
cgas = boltzm*avo                   # ideal gas constant; erg/K
kev = 8.617385e-5                   # converts temp to ev (ev K^-1)
amu = 1./avo                        # atomic mass unit (g)
mn = 1.6749286e-24                  # neutron mass (g)
mp = 1.6726231e-24                  # proton mass (g)
me = 9.1093826e-28                  # electron mass (g)
rbohr = hbar*hbar/(me * qe * qe)    # Bohr radius (cm)
fine = qe*qe/(hbar*clight)          # fine structure constant
hion = 13.605698140                 # hydrogen ionization energy (eV)
ev2erg = 1.602176565e-12            # electron volt to erg
inversecm2erg = planck_h * clight   # cm^-1 to erg
mev_to_ergs = 1e6 * ev2erg
mev_amu = mev_to_ergs / amu
Qconv = mev_to_ergs*avo
boltz_sigma = 5.670373e-5   # boltzmann's sigma = a*c/4 (erg cm^-2 K^-4 s^-1)
# radiation density constant, a ~ 7.5657e-15 erg cm^-3 K^-4; Prad=crad*T^4 / 3
crad = boltz_sigma * 4 / clight
ssol = boltz_sigma
asol = crad
weinlam = planck_h * clight / (kerg * 4.965114232)
weinfre = 2.821439372*kerg/planck_h
rhonuc = 2.342e14                   # density of nucleus (g cm^3)


# ASTRONOMICAL CONSTANTS.
# Solar age, L, and R values from Bahcall et al, ApJ 618 (2005) 1049-1056.

Zsun = 0.0142                   # Asplund+09
# Y abundances @ ZAMS given Aplund+09 scaling, w/ keys as Z/Zsun
zams_table = {2.: 2.915e-01,
              1.: 2.703e-01,
              0.45: 2.586e-01,
              0.2: 2.533e-01,
              0.1: 2.511e-01,
              0.01: 2.492e-01,
              0.001: 2.49e-01,
              0.0001: 2.49e-01}
msol = 1.9892e33                # solar mass, gravitational not baryonic (g)
rsol = 6.9598e10                # solar radius (cm)
lsol = 3.8418e33                # solar luminosity (erg s^-1)
agesol = 4.57e9                 # solar age (years)
Msun = msol
Rsun = rsol
Lsun = lsol
Msun33 = msol*1e-33
Rsun11 = rsol*1e-11
Lsun33 = lsol*1e-33
ly = 9.460528e17                # light year (cm)
pc = 3.261633 * ly              # parsec (cm)
secyer = 3.1558149984e7         # seconds per year
dayyer = 365.25                 # days per year
age_of_universe = 13.8e9        # (year)

Teffsol = 5777.0
loggsol = 4.4378893534131256    # w. default MESA msol, rsol and standard_cgrav

mbolsun = 4.746                 # Bolometric magnitude of the Sun
mbolsol = mbolsun

m_earth = 5.9764e27             # earth mass; 3.004424e-6 Msun (g)
r_earth = 6.37e8                # earth radius (cm)
au = 1.495978921e13             # astronomical unit (cm)
aursun = 214.95

m_jupiter = 1.8986e30           # jupiter mass, 0.954454e-3 Msun (g)
r_jupiter = 6.9911e9            # jupiter mean radius (cm)
semimajor_axis_jupiter = 7.7857e13  # jupiter semimajor axis (cm)


# CONVERSIONS (not in MESA)

km2cm = 1.0e5                   # Km (in cm)
day2sec = 24.0*60.0*60.0        # day (in seconds)
H_weight = 1.00782503207        # atomic mass; see TABLE II
He_weight = 4.002603254131      # atomic mass; see TABLE III
SNcheck_ERR = 1e-10             # machine precision
