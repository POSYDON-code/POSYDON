"""Common functions to be used while running populations."""


__authors__ = [
    "Konstantinos Kovlakas <Konstantinos.Kovlakas@unige.ch>",
    "Devina Misra <devina.misra@unige.ch>",
    "Emmanouil Zapartas <ezapartas@gmail.com>",
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Nam Tran <tranhn03@gmail.com>",
    "Ying Qin <<yingqin2013@hotmail.com>",
    "Jeffrey Andrews <jeffrey.andrews@northwestern.edu>",
    "Tassos Fragos <Anastasios.Fragkos@unige.ch>",
    "Scott Coughlin <scottcoughlin2014@u.northwestern.edu>",
    "Kyle Akira Rocha <kylerocha2024@u.northwestern.edu>",
]


import os
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import newton
from scipy.integrate import quad
from posydon.utils import constants as const
import copy
import warnings
from scipy.interpolate import PchipInterpolator


PATH_TO_POSYDON = os.environ.get("PATH_TO_POSYDON")


# Constants related to inferring star states
THRESHOLD_CENTRAL_ABUNDANCE = 0.01   # central abundance for flagging depletion
THRESHOLD_HE_NAKED_ABUNDANCE = 0.01  # for surface abundance for stripped_He
THRESHOLD_NUCLEAR_LUMINOSITY = 0.97  # for element fraction in nuclear burning
# relative burning threshold with respect to nuclear luminosity
REL_LOG10_BURNING_THRESHOLD = np.log10(1.0 - THRESHOLD_NUCLEAR_LUMINOSITY)
LOG10_BURNING_THRESHOLD = -10.0      # burning luminosity threshold (in Lsol)
STATE_NS_STARMASS_UPPER_LIMIT = 2.5
STATE_UNDETERMINED = "undetermined_evolutionary_state"

# ALL POSSIBLE STAR STATES
BURNING_STATES = ["Core_H_burning", "Core_He_burning",
                  "Shell_H_burning", "Central_He_depleted",
                  "Central_C_depletion"]
RICHNESS_STATES = ["H-rich", "stripped_He"]
COMPACT_OBJECTS = ["WD", "NS", "BH","massless_remnant"]

ALL_STAR_STATES = COMPACT_OBJECTS + [STATE_UNDETERMINED]
ALL_STAR_STATES.extend(["{}_{}".format(rich_in, burning)
                        for rich_in in RICHNESS_STATES
                        for burning in BURNING_STATES])

# `ALL_STAR_STATES` includes the following strings:
# 'WD', 'NS', 'BH', 'undetermined_evolutionary_state',
# 'H-rich_Core_H_burning', 'H-rich_Core_He_burning',
# 'H-rich_Shell_H_burning', 'H-rich_Central_He_depleted',
# 'stripped_He_Core_He_burning', 'stripped_Central_He_depleted',
# 'H-rich_Central_C_depletion', 'stripped_He_Central_C_depletion'

# Mass-transfer cases in form of integer flags
MT_CASE_NO_RLO = 0
MT_CASE_A = 1
MT_CASE_B = 2
MT_CASE_C = 3
MT_CASE_BA = 4
MT_CASE_BB = 5
MT_CASE_BC = 6
MT_CASE_NONBURNING = 8
MT_CASE_UNDETERMINED = 9

# All cases meaning RLO is happening
ALL_RLO_CASES = set([MT_CASE_A, MT_CASE_B, MT_CASE_C,
                     MT_CASE_BA, MT_CASE_BB, MT_CASE_BC,
                     MT_CASE_NONBURNING])

# Conversion of integer mass-transfer flags to strings
MT_CASE_TO_STR = {
    MT_CASE_NO_RLO: "no_RLO",
    MT_CASE_A: "A",
    MT_CASE_B: "B",
    MT_CASE_C: "C",
    MT_CASE_BA: "BA",
    MT_CASE_BB: "BB",
    MT_CASE_BC: "BC",
    MT_CASE_NONBURNING: "nonburning",
    MT_CASE_UNDETERMINED: "undetermined_MT"
}

# Conversion of strings to integer mass-transfer flags
MT_STR_TO_CASE = {string: integer for integer, string
                  in MT_CASE_TO_STR.items()}

# Threshold used for inferring mass-transfer, and RLO
RL_RELATIVE_OVERFLOW_THRESHOLD = -0.05
LG_MTRANSFER_RATE_THRESHOLD = -12


DEFAULT_CE_OPTION_FOR_LAMBDA = \
    "lambda_from_profile_gravitational_plus_internal_minus_recombination"


def is_number(s):
    """Check if the input can be converted to a float."""
    try:
        float(s)
        return True
    except ValueError:
        return False


def stefan_boltzmann_law(L, R):
    """Compute the effective temperature give the luminosity and radius."""
    return (L * const.Lsun / (4.0 * np.pi * (R*const.Rsun) ** 2.0)
            / const.boltz_sigma) ** (1.0 / 4.0)


def rzams(m, z=0.02, Zsun=0.02):
    """Evaluate the zero age main sequence radius [1]_.

    Parameters
    ----------
    m : array_like
        The masses of the stars in Msun.
    z : float
        The metallicity of the star.

    Returns
    -------
    ndarray
        Array of same size as `m` containing the ZAMS radii of the stars (Rsun)

    References
    ----------
    .. [1] Tout C. A., Pols O. R., Eggleton P. P., Han Z., 1996, MNRAS, 281,
        257


    """
    m = np.asanyarray(m)
    xz = [
        0.0, 3.970417e-01, -3.2913574e-01, 3.4776688e-01, 3.7470851e-01,
        9.011915e-02, 8.527626e+00, -2.441225973e+01, 5.643597107e+01,
        3.706152575e+01, 5.4562406e+00, 2.5546e-04, -1.23461e-03, -2.3246e-04,
        4.5519e-04, 1.6176e-04, 5.432889e+00, -8.62157806e+00, 1.344202049e+01,
        1.451584135e+01, 3.39793084e+00, 5.563579e+00, -1.032345224e+01,
        1.944322980e+01, 1.897361347e+01, 4.16903097e+00, 7.8866060e-01,
        -2.90870942e+00, 6.54713531e+00, 4.05606657e+00, 5.3287322e-01,
        5.86685e-03, -1.704237e-02, 3.872348e-02, 2.570041e-02, 3.83376e-03,
        1.715359e+00, 6.2246212e-01, -9.2557761e-01, -1.16996966e+00,
        -3.0631491e-01, 6.597788e+00, -4.2450044e-01, -1.213339427e+01,
        -1.073509484e+01, -2.51487077e+00, 1.008855000e+01, -7.11727086e+00,
        -3.167119479e+01, -2.424848322e+01, -5.33608972e+00, 1.012495e+00,
        3.2699690e-01, -9.23418e-03, -3.876858e-02, -4.12750e-03, 7.490166e-02,
        2.410413e-02, 7.233664e-02, 3.040467e-02, 1.97741e-03, 1.077422e-02,
        3.082234e+00, 9.447205e-01, -2.15200882e+00, -2.49219496e+00,
        -6.3848738e-01, 1.784778e+01, -7.4534569e+00, -4.896066856e+01,
        -4.005386135e+01, -9.09331816e+00, 2.2582e-04, -1.86899e-03,
        3.88783e-03, 1.42402e-03, -7.671e-05
    ]
    lzs = np.log10(z / Zsun)

    msp = np.zeros(17)
    msp[0] = 0.0
    msp[1] = xz[1] + lzs * (xz[2] + lzs * (xz[3] + lzs
                                           * (xz[4] + lzs * xz[5])))
    msp[2] = xz[6] + lzs * (xz[7] + lzs * (xz[8] + lzs
                                           * (xz[9] + lzs * xz[10])))
    msp[3] = xz[11] + lzs * (xz[12] + lzs * (xz[13] + lzs
                                             * (xz[14] + lzs * xz[15])))
    msp[4] = xz[16] + lzs * (xz[17] + lzs * (xz[18] + lzs
                                             * (xz[19] + lzs * xz[20])))
    msp[5] = xz[21] + lzs * (xz[22] + lzs * (xz[23] + lzs
                                             * (xz[24] + lzs * xz[25])))
    msp[6] = xz[26] + lzs * (xz[27] + lzs * (xz[28] + lzs
                                             * (xz[29] + lzs * xz[30])))
    msp[7] = xz[31] + lzs * (xz[32] + lzs * (xz[33] + lzs
                                             * (xz[34] + lzs * xz[35])))
    msp[8] = xz[36] + lzs * (xz[37] + lzs * (xz[38] + lzs
                                             * (xz[39] + lzs * xz[40])))
    msp[9] = xz[41] + lzs * (xz[42] + lzs * (xz[43] + lzs
                                             * (xz[44] + lzs * xz[45])))
    msp[10] = xz[46] + lzs * (xz[47] + lzs * (xz[48] + lzs
                                              * (xz[49] + lzs * xz[50])))
    msp[11] = xz[51] + lzs * (xz[52] + lzs * (xz[53] + lzs
                                              * (xz[54] + lzs * xz[55])))
    msp[12] = xz[56] + lzs * (xz[57] + lzs * (xz[58] + lzs
                                              * (xz[59] + lzs * xz[60])))
    msp[13] = xz[61]
    msp[14] = xz[62] + lzs * (xz[63] + lzs * (xz[64] + lzs
                                              * (xz[65] + lzs * xz[66])))
    msp[15] = xz[67] + lzs * (xz[68] + lzs * (xz[69] + lzs
                                              * (xz[70] + lzs * xz[71])))
    msp[16] = xz[72] + lzs * (xz[73] + lzs * (xz[74] + lzs
                                              * (xz[75] + lzs * xz[76])))
    mx = np.sqrt(m)
    r = ((msp[8] * m**2 + msp[9] * m**6) * mx + msp[10] * m**11
         + (msp[11] + msp[12] * mx) * m**19) / (
             msp[13] + msp[14] * m**2
             + (msp[15] * m**8 + m**18 + msp[16] * m**19) * mx)

    return r


'''


Receives:
q ->
a_orb ->

Returns:
RL -> Roche lobe radius in similar units as a_orb
'''


def roche_lobe_radius(q, a_orb=1):
    """Approximate the Roche lobe radius from [1]_.

    Parameters
    ----------
    q : float
        Dimensionless mass ratio = MRL/Mcomp, where
        MRL is the mass of the star we calculate the RL and
        Mcomp is the mass of its companion star.
    a_orb : float
        Orbital separation. The return value will have the same unit.

    Returns
    -------
    float
        Roche lobe radius in similar units as a_orb
    References
    ----------
    .. [1] Eggleton, P. P. 1983, ApJ, 268, 368

    """
    RL = a_orb * (0.49 * q**(2. / 3.)) / (
        0.6 * q**(2. / 3.) + np.log(1 + q**(1. / 3))
    )
    return RL


def orbital_separation_from_period(period_days, m1_solar, m2_solar):
    """Apply the Third Kepler law.

    Parameters
    ----------
    period_days : float
        Orbital period in days.
    m1_solar : float
        Mass of the one of the stars, in solar units.
    m2_solar : float
        Mass of the other star, in solar units.

    Returns
    -------
    float
        The separation of the binary in solar radii.

    """
    # cast to float64 to avoid overflow
    m1_solar = np.float64(m1_solar)
    m2_solar = np.float64(m2_solar)
    period_days = np.float64(period_days)
    
    separation_cm = (const.standard_cgrav
                     * (m1_solar * const.Msun + m2_solar * const.Msun)
                     / (4.0 * const.pi**2.0)
                     * (period_days * const.day2sec) ** 2.0)**(1.0/3.0)
    separation_solar = separation_cm / const.Rsun
    return separation_solar


def orbital_period_from_separation(separation, m1, m2):
    """Apply the Third Kepler law.

    Parameters
    ----------
    separation : float
        Orbital separation (semi-major axis) in Rsun.
    m1 : float
        Mass of one of the stars in solar units.
    m2 : type
        Mass of the other star in solar units.

    Returns
    -------
    float
        The orbital period in days.

    """
    return const.dayyer * ((separation / const.aursun)**3.0 / (m1 + m2)) ** 0.5


def BSE_to_POSYDON(ktype):
    """Convert BSE numerical type to POSYDON string.

    Parameters
    ----------
    ktype : int
        The BSE numerical type.
    core_mass : float
        Core mass of the star.

    Returns
    -------
    str
        The corresponding POSYDON string.

    """
    if ktype in (0, 1):
        state = 'H-rich_Core_H_burning'
    elif ktype in (2, 3):   # Hertzprung Gap and Giant Branch star
        state = 'H-rich_Shell_H_burning'
    elif ktype == 4:        # Core Helium Burning
        state = 'H-rich_Core_He_burning'
    elif ktype in (5, 6):   # Early AGB and TPAGB
        state = 'H-rich_Central_He_depleted'
    elif ktype == 7:        # Naked Helium Star MS
        state = 'stripped_He_Core_He_burning'
    elif ktype in (8, 9):   # Naked Helium star (HG and GB)
        state = 'stripped_He_Central_He_depleted'
    elif ktype in (10, 11, 12):
        state = 'WD'
    elif ktype == 13:
        state = 'NS'
    elif ktype == 14:
        state = 'BH'
    elif ktype == 15:
        state = 'Massless remnant'
    else:
        raise ValueError('Conversion of the ktype {} is not supported.'.
                         format(ktype))
    return state


def POSYDON_to_BSE(star):
    """Convert POSYDON state to BSE numerical type.

    Parameters
    ----------
    star : SingleStar
        The star which state requires conversion.

    Returns
    -------
    int
        The corresponding BSE numerical type.
    """
    if star.state == 'H-rich_Core_H_burning':
        if star.mass < 0.7:
            ktype = 0
        else:
            ktype = 1
    elif star.state == 'NS':
        ktype = 13
    elif star.state == 'BH':
        ktype = 14
    else:
        raise ValueError(
            'Conversion of the state {} is currently not supported.'.format(
                star.state))
    return ktype


def eddington_limit(binary, idx=-1):
    """Calculate the Eddington limit & radtiative efficiency of compact object.

    Parameters
    ----------
    binary : BinaryStar
        The binary object.


    Returns
    -------
    list
        The Eddington accretion limit and radiative efficiency in solar units.

    """
    if binary.star_1.state in ['NS', 'BH', 'WD']:
        accretor = binary.star_1
        donor = binary.star_2
    elif binary.star_2.state in ['NS', 'BH', 'WD']:
        accretor = binary.star_2
        donor = binary.star_1
    else:
        raise Exception("Eddington-limit to be calculated for non-CO?")

    state_acc = np.atleast_1d(
        np.asanyarray([*accretor.state_history, accretor.state])[idx])
    m_acc = np.atleast_1d(
        np.asanyarray([*accretor.mass_history, accretor.mass],
                      dtype=float)[idx]) * const.msol
    surface_h1 = np.atleast_1d(
        np.asanyarray([*donor.surface_h1_history, donor.surface_h1],
                      dtype=float)[idx])

    for i in range(len(state_acc)):
        if state_acc[i] is None:
            if any(j == 'NS' for j in accretor.state_history):
                state_acc[i] = 'NS'
            elif any(j == 'BH' for j in accretor.state_history):
                state_acc[i] = 'BH'
            elif any(j == 'WD' for j in accretor.state_history):
                state_acc[i] = 'WD'
            else:
                raise ValueError('COtype must be "BH", "NS", or "WD"')

        if surface_h1[i] is None:
            surface_h1[i] = 0.7155
        if state_acc[i] == "BH":
            r_isco = 6
            # m_ini is the accretor mass at zero spin
            m_ini = m_acc[i] * np.sqrt(r_isco / 6)
            eta = 1 - np.sqrt(1 - (min(m_acc[i],
                                       np.sqrt(6) * m_ini) / (3 * m_ini))**2)
        else:
            # 1.25 * 10**6 cm
            acc_radius = CO_radius(m_acc[i], state_acc[i]) * const.Rsun
            eta = const.standard_cgrav*m_acc[i]/(acc_radius*const.clight**2)

        # This is the mass accretion rate corresponding to the
        # Eddington luminosity, L_edd = eta * mdot_edd * clightˆ2 (cgs units)
        mdot_edd = (4 * np.pi * const.standard_cgrav * m_acc[i]
                    / (0.2 * (1 + surface_h1[i]) * eta * const.clight))
    return mdot_edd / (const.msol / const.secyer), eta


def beaming(binary):
    """Calculate the geometrical beaming of a super-Eddington accreting source
    [1]_, [2]_.

    Compute the super-Eddington isotropic-equivalent accretion rate and the
    beaming factor of a star. This does not change the intrinsic accretion onto
    the accretor and is an observational effect due to the inflated structure
    of the accretion disc that beams the outgoing X-ray emission. This is
    important for observing super-Eddington X-ray sources
    (e.g. ultraluminous X-ray sources). In case of a BH we are assuming that it
    has zero spin which is not a good approximation for high accretion rates.

    Parameters
    ----------
    binary : BinaryStar
        The binary object.

    Returns
    -------
    list
        The super-Eddington isotropic-equivalent accretion rate and beaming
        factor respcetively in solar units.

    References
    ----------
    .. [1] Shakura, N. I. & Sunyaev, R. A. 1973, A&A, 24, 337
    .. [2] King A. R., 2008, MNRAS, 385, L113

    """
    mdot_edd = eddington_limit(binary, idx=-1)[0]

    rlo_mdot = 10**binary.lg_mtransfer_rate

    if rlo_mdot >= mdot_edd:
        if rlo_mdot > 8.5 * mdot_edd:
            # eq. 8 in King A. R., 2009, MNRAS, 393, L41-L44
            b = 73 / (rlo_mdot / mdot_edd)**2
        else:
            b = 1
        # Shakura, N. I. & Sunyaev, R. A. 1973, A&A, 24, 337
        mdot_beam = mdot_edd * (1 + np.log(rlo_mdot / mdot_edd)) / b
    else:
        b = 1
        mdot_beam = 10**binary.lg_mtransfer_rate

    return mdot_beam, b


def bondi_hoyle(binary, accretor, donor, idx=-1, wind_disk_criteria=True,
                scheme='Hurley+2002'):
    """Calculate the Bondi-Hoyle accretion rate of a binary [1]_.

    Parameters
    ----------
    binary : BinaryStar
        The binary which accretion rate is required.
    accretor : SingleStar
        The accretor in the binary.
    donor : SingleStar
        The donor in the binary.
    idx : int
        default: -1
    wind_disk_criteria : bool
        default: True, see [5]_
    scheme : str
        There are different options:

        - 'Hurley+2002' : following [3]_
        - 'Kudritzki+2000' : following [7]_

    Returns
    -------
    float
        The Bondi-Hoyle accretion rate in solar units.

    Notes
    -----
    An approximation is used for the accretion rate [2]_ and the wind velocity
    of the donor is moddeled as in [3]_, [6]_. Also see [4]_.

    References
    ----------
    .. [1] Bondi, H., & Hoyle, F. 1944, MNRAS, 104, 273
    .. [2] Boffin, H. M. J., & Jorissen, A. 1988, A&A, 205, 155
    .. [3] Hurley, J. R., Tout, C. A., & Pols, O. R. 2002, MNRAS, 329, 897
    .. [4] Belczynski, K., Kalogera, V., Rasio, F. A., et al. 2008, ApJS, 174,
        223
    .. [5] Sen, K. ,Xu, X. -T., Langer, N., El Mellah, I. , Schurmann, C., &
        Quast, M., 2021, A&A
    .. [6] Sander A. A. C., Vink J. S., 2020, MNRAS, 499, 873
    .. [7] Kudritzki, R.-P., & Puls, J. 2000, ARA&A, 38, 613

    """
    alpha = 1.5
    G = const.standard_cgrav * 1e-3     # 6.67428e-11 m3 kg-1 s-2
    Msun = const.Msun * 1e-3            # 1.988547e30  kg
    Rsun = const.Rsun * 1e-2            # 6.9566e8 m

    sep = np.atleast_1d(
        np.asanyarray([*binary.separation_history, binary.separation],
                      dtype=float)[idx])
    ecc = np.atleast_1d(
        np.asanyarray([*binary.eccentricity_history, binary.eccentricity],
                      dtype=float)[idx])
    m_acc = np.atleast_1d(
        np.asanyarray([*accretor.mass_history, accretor.mass],
                      dtype=float)[idx])
    m = np.atleast_1d(
        np.asanyarray([*donor.mass_history, donor.mass], dtype=float)[idx])
    lg_mdot = np.atleast_1d(
        np.asanyarray([*donor.lg_wind_mdot_history, donor.lg_wind_mdot],
                      dtype=float)[idx])
    he_core_mass = np.atleast_1d(
        np.asanyarray([*donor.he_core_radius_history, donor.he_core_radius],
                      dtype=float)[idx])
    log_R = np.atleast_1d(
        np.asanyarray([*donor.log_R_history, donor.log_R], dtype=float)[idx])
    radius = 10**log_R  # Rsun
    surface_h1 = np.atleast_1d(
        np.asanyarray([*donor.surface_h1_history, donor.surface_h1],
                      dtype=float)[idx])
    L = np.atleast_1d(
        np.asanyarray([*donor.log_L_history, donor.log_L], dtype=float)[idx])
    Teff = stefan_boltzmann_law(10**L, radius)

    beta = np.empty_like(sep)

    # Hurley, J. R., Tout, C. A., & Pols, O. R. 2002, MNRAS, 329, 897
    if scheme == 'Hurley+2002':
        # For H-rich stars
        beta[np.logical_and(he_core_mass, radius > 900.0)] = 0.125
        beta[m > 120.0] = 7.0
        beta[m < 1.4] = 0.5
        cond = np.logical_and(m >= 1.4, m <= 120.0)
        beta[cond] = 0.5 + (m[cond] - 1.4) / (120.0 - 1.4) * (6.5)

        # For He-rich stars
        beta[np.logical_and(surface_h1 <= 0.01, m > 120.0)] = 7.0
        beta[np.logical_and(surface_h1 <= 0.01, m < 10.0)] = 0.125
        mass_cond = np.logical_and(m >= 10.0, m <= 120.0)
        he_cond = np.logical_and(mass_cond, surface_h1 <= 0.01)
        beta[he_cond] = 0.125 + (m[he_cond] - 10.0) / (120.0 - 10.0) * (6.875)

        f_m = np.sqrt(beta)

    # Kudritzki, R.-P., & Puls, J. 2000, ARA&A, 38, 613
    elif scheme == 'Kudritzki+2000':
        for i in range(len(m)):
            if Teff[i] >= 21000:
                f_m = 2.65
            elif Teff[i] <= 10000:
                f_m = 1.0
            else:
                f_m = 1.4

    v_esc = np.sqrt(2 * G * m * Msun / (radius * Rsun))     # m/s
    v_wind = v_esc * f_m                                    # m/s

    # Sander A. A. C., Vink J. S., 2020, MNRAS, 499, 873
    for i in range(len(m)):
        if (surface_h1[i] < 0.4 and Teff[i] > 1.0e4):
            if lg_mdot[i] >= -5.25:
                slope = (3.7 - 3.25) / (-2.5 + 5.25)
            else:
                slope = (3.25 - 3.75) / (-5.25 + 7.25)
            v_wind[i] = 10 ** (slope * lg_mdot[i] + 3.25 + 5.25 * slope) * 1000
        else:
            pass

    n = np.sqrt((G * (m_acc + m) * Msun) / ((radius * Rsun)**3))
    t0 = np.random.rand(len(sep)) * 2 * np.pi / n
    E = newton(lambda x: x - ecc * np.sin(x) - n * t0,
               np.ones_like(sep) * np.pi / 2,
               maxiter=100)

    b = sep * Rsun * np.sqrt(1 - ecc**2)
    r_vec = np.array([sep * Rsun * (np.cos(E) - ecc), b * np.sin(E)])
    v_dir = np.array([-sep * Rsun * np.sin(E), b * np.cos(E)])
    r = np.linalg.norm(r_vec, axis=0)
    v_dir_norm = np.linalg.norm(v_dir, axis=0)

    k = np.einsum('ij,ij->j', r_vec, v_dir) / (r * v_dir_norm)  # cos(angle)
    v = np.sqrt(G * (m + m_acc) * Msun * ((2 / r) - (1 / (sep * Rsun))))  # m/s
    v_rel = np.sqrt(v**2 + v_wind**2 + 2 * v * v_wind * k)                # m/s

    # Bondi, H., & Hoyle, F. 1944, MNRAS, 104, 273
    mdot_acc = alpha * ((G * m_acc * Msun)**2
                        / (2 * v_rel**3 * v_wind * r**2)) * 10**lg_mdot

    # eq. 10 in Sen, K. ,Xu, X. -T., Langer, N., El Mellah, I. , Schurmann, C.,
    # Quast, M., 2021, A&A
    if wind_disk_criteria:      # check if accretion disk will form
        eta = 1.0/3.0           # wind accretion efficiency between 1 and 1/3
        gamma = 1.0             # for non-spinning BH
        q = m / m_acc
        rdisk_div_risco = (
            (2/3) * (eta / (1 + q)) ** 2
            * (v / (const.clight * 0.01)) ** (-2)
            * (1 + (v_wind / v) ** 2) ** (-4) * gamma ** (-1))
        for i in range(len(rdisk_div_risco)):
            if rdisk_div_risco[i] <= 1:         # No disk formed
                mdot_acc[i] = 10**-99.0

    # make it Eddington-limited
    mdot_edd = eddington_limit(binary, idx=idx)[0]
    mdot_acc = np.minimum(mdot_acc, mdot_edd)

    return np.squeeze(mdot_acc)


def rejection_sampler(x=None, y=None, size=1, x_lim=None, pdf=None):
    """Generate a sample from a 1d PDF using the acceptance-rejection method.

    Parameters
    ----------
    x : array_like
        The x-values of the PDF.
    y : array_like
        The y-values of the PDF.
    size : int
        The number of random numbers to generate.
    x_lim : array float
        The boundary where the pdf function is defined if passed as pdf.
    pdf : func
        The pdf function

    Returns
    -------
    ndarray
        An array with `size` random numbers generated from the PDF.

    """
    if pdf is None:
        assert np.all(y >= 0.0)
        try:
            pdf = interp1d(x, y)
        except ValueError:
            idxs = np.argsort(x)
            pdf = interp1d(x.take(idxs), y.take(idxs))

        x_rand = np.random.uniform(x.min(), x.max(), size)
        y_rand = np.random.uniform(0, y.max(), size)
        values = x_rand[y_rand <= pdf(x_rand)]
        while values.shape[0] < size:
            n = size - values.shape[0]
            x_rand = np.random.uniform(x.min(), x.max(), n)
            y_rand = np.random.uniform(0, y.max(), n)
            values = np.hstack([values, x_rand[y_rand <= pdf(x_rand)]])
    else:
        x_rand = np.random.uniform(x_lim[0], x_lim[1], size)
        pdf_max = max(pdf(np.random.uniform(x_lim[0], x_lim[1], 50000)))
        y_rand = np.random.uniform(0, pdf_max, size)
        values = x_rand[y_rand <= pdf(x_rand)]
        while values.shape[0] < size:
            n = size - values.shape[0]
            x_rand = np.random.uniform(x_lim[0], x_lim[1], n)
            y_rand = np.random.uniform(0, pdf_max, n)
            values = np.hstack([values, x_rand[y_rand <= pdf(x_rand)]])

    return values


def inverse_sampler(x, y, size=1):
    """Sample from a PDF using the inverse-sampling method.

    Parameters
    ----------
    x : array-like
        The x-axis coordinates of the points where the PDF is defined.
    y : array-like
        The probablity density at `x` (or a scaled version of it).
    size : int
        Number of samples to generate.

    Returns
    -------
    array
        The sample drawn from the PDF.

    """
    # x should be sorted
    assert np.all(np.diff(x) > 0)
    # y should be above 0
    assert np.all(y >= 0.0)

    # compute the area of each trapezoid
    segment_areas = 0.5 * (y[1:]+y[:-1]) * (x[1:]-x[:-1])
    # their cumulative sum denotes the scaled CDF at each x value
    cumsum_areas = np.cumsum(segment_areas)
    total_area = cumsum_areas[-1]

    # start the inverse sampling
    u = np.random.uniform(size=size)
    # index of "bin" where each sampled value corresponds too
    u_indices = np.searchsorted(cumsum_areas, u * total_area)
    # the area that should be covered from the end of the previous bin
    delta_y = total_area * u - cumsum_areas[u_indices-1]
    delta_y[u_indices == 0] = total_area * u[u_indices == 0]

    # the width and height (of the cap) of each sample's bin
    dx_bins = x[u_indices+1] - x[u_indices]
    dy_bins = y[u_indices+1] - y[u_indices]

    sample = x[u_indices] + dx_bins * (
        (y[u_indices]**2.0
         + 2.0 * delta_y * dy_bins/dx_bins)**0.5 - y[u_indices]) / dy_bins

    # if nan values found, then flat CDF for which the inverse is undefined...
    where_nan = np.where(~np.isfinite(sample))
    n_where_nan = len(where_nan)
    # ... in that case, simply sample randomly from each flat bin!
    if n_where_nan:
        assert np.all(dy_bins[where_nan] == 0)
        sample[where_nan] = x[u_indices][where_nan] + (
            dx_bins[where_nan] * np.random.uniform(size=n_where_nan))

    # make sure that everything worked as expected
    assert np.all(np.isfinite(sample))
    return sample


def histogram_sampler(x_edges, y, size=1):
    """Sample from an empirical PDF represented by a histogram.

    Parameters
    ----------
    x_edges : array-like
        The edges of the bins of the histrogram.
    y : array-like
        The counts (or a scaled version of them) of the histogram.
    size : int
        Number of random values to produce.

    Returns
    -------
    array
        The random sample.

    """
    assert np.all(y >= 0.0)

    # make sure that the lengths of the input arrays are correct
    n_bins = len(y)
    assert n_bins > 0 and len(x_edges) == n_bins + 1
    # first decide which will be the bin of each element in the sample
    bin_sample = np.random.choice(n_bins, replace=True, p=y/sum(y), size=size)

    sample = np.ones(size) * np.nan

    # select each bin and based on its uniform distribution, decide the sample
    bins_found = set(bin_sample)
    for bin_index in bins_found:
        in_this_bin = np.argwhere(bin_sample == bin_index)[:, 0]
        sample[in_this_bin] = np.random.uniform(
            x_edges[bin_index], x_edges[bin_index+1], size=len(in_this_bin))

    assert(np.all(np.isfinite(sample)))
    return np.squeeze(sample)


def read_histogram_from_file(path):
    """Read a histogram from a CSV file.

    The expected format is:

    # comment line
    x[0], x[1], ..., x[n], x[n+1]
    y[0], y[1], ..., y[n]

    where # denotes a comment line (also empty lines are ignored), and `n` is
    the number of bins (notice that the first line contains n+1 elements.)

    Usage: bin_edges, bin_counts = read_histogram_from_file("a_histogram.csv").


    Parameters
    ----------
    path : str
        The path of the CSV file containing the histogram information.

    Returns
    -------
    list of arrays
        The bin edges and bin counts of the histogram.

    """
    with open(path, "r") as f:
        arrays = []
        for line in f:
            line = line.strip()
            if len(line) == 0 or line.startswith("#"):
                continue
            arrays.append(np.fromstring(line.strip(), dtype=float, sep=","))
            if len(arrays) > 2:
                raise Exception("More than two lines found in the document.")

    return arrays


def inspiral_timescale_from_separation(star1_mass, star2_mass,
                                       separation, eccentricity):
    """Compute the timescale of GW inspiral using the orbital separation.

    Based on [1]_:
        https://journals.aps.org/pr/abstract/10.1103/PhysRev.136.B1224

    Parameters
    ----------
    star1_mass : float
        Mass of star 1 in Msun.
    star2_mass : float
        Mass of star 2 in Msun.
    separation : type
        Binary separation in Rsun.
    eccentricity : type
        Eccentricity of the binary orbit. Must be 0 <= ecc <1.

    Returns
    -------
    float
        The inspiral time scale of the two black holes.

    References
    ----------
    .. [1] Peters 1964 Phys. Rev. 136, B1224

    """
    # NOTE we should check that this matches up with Maggiori equations
    G = const.standard_cgrav
    c = const.clight
    Msun = const.Msun
    Rsun = const.Rsun
    secyer = const.secyer

    m1 = star1_mass * Msun
    m2 = star2_mass * Msun
    a = separation * Rsun
    ecc = eccentricity

    if m1 <= 0:
        raise ValueError(
            "Mass of star 1 is <= 0, which is not a physical value."
        )
    if m2 <= 0:
        raise ValueError(
            "Mass of star 2 is <= 0, which is not a physical value."
        )
    if a <= 0:
        raise ValueError(
            "Separation is <= 0, which is not a physical value.")

    if ecc < 0:
        raise ValueError(
            "Eccentricity is < 0, which is not a physical value.")
    if ecc >= 1:
        raise ValueError(
            "Eccentricity is >= 1, which is not a physical value."
        )

    # Eq. (5.9) in Peters 1964 Phys. Rev. 136, B1224
    beta = (64.0 / 5) * (G**3) * m1 * m2 * (m1 + m2) / (c**5)

    if ecc == 0:
        # Eq. (5.10) in Peters 1964 Phys. Rev. 136, B1224
        T_merger = a**4 / (4 * beta)
    else:
        # Eq. (5.11) at e_0, i.e. a_0 = a(e_0), solved for c_0 in
        # Peters 1964 Phys. Rev. 136, B1224
        def c_0(a0, e0):
            return ((a0 * (1 - e0**2)) * (e0**(-12.0 / 19))
                    * (1 + (121.0 / 304) * e0**2)**(-870.0 / 2299))

        c0 = c_0(a, ecc)

        # Eq. (5.14)
        def integrand(e):
            return (e**(29. / 19) * (1 + (121. / 304) * e**2)**(1181.0 / 2299)
                    / (1 - e**2)**(3.0 / 2))

        # assume binary circularizes by the time it merges
        T_merger = (12.0 / 19) * (
            (c0**4) / beta) * quad(integrand, 0.0, ecc)[0]

    # return T_merge in Myr
    return T_merger / (secyer * 1e6)


def inspiral_timescale_from_orbital_period(star1_mass, star2_mass,
                                           orbital_period, eccentricity):
    """Compute the timescale of GW inspiral using the orbital period.

    Based on [1]_:
        https://journals.aps.org/pr/abstract/10.1103/PhysRev.136.B1224

    Parameters
    ----------
    star1_mass : float
        Mass of star 1 in Msun.
    star2_mass : float
        Mass of star 2 in Msun.
    orbital_separation : type
        Binary separation in Rsun.
    eccentricity : type
        Eccentricity of the binary orbit. Must be 0 <= ecc <1.

    Returns
    -------
    float
        The inspiral time scale of the two black holes in Myr.

    References
    ----------
    .. [1] Peters 1964 Phys. Rev. 136, B1224

    """
    # NOTE we should check that this matches up with Maggiori equations
    separation = orbital_separation_from_period(orbital_period, star1_mass,
                                                star2_mass)
    T_merge = inspiral_timescale_from_separation(star1_mass, star2_mass,
                                                 separation, eccentricity)
    return T_merge


def spin_stable_mass_transfer(spin_i, star_mass_preMT, star_mass_postMT):
    """Calculate the spin of an accreting BH under stable mass transfer.

    Based on Thorne 1974 eq. 2a.

    """
    if star_mass_preMT is None or star_mass_postMT is None:
        return None
    z1 = 1+(1-spin_i**2)**(1/3)*((1+spin_i)**(1/3)+(1-spin_i)**(1/3))
    z2 = (3*spin_i**2+z1**2)**0.5
    r_isco = 3 + z2 - ((3-z1)*(3+z1+2*z2))**0.5
    if (1 <= star_mass_postMT / star_mass_preMT
            and star_mass_postMT / star_mass_preMT <= r_isco**0.5):
        spin = r_isco**(0.5) / 3 * (star_mass_preMT / star_mass_postMT) * (
            4 - (3 * r_isco * star_mass_preMT**2 / star_mass_postMT**2 - 2)**(0.5))
    elif star_mass_postMT / star_mass_preMT > r_isco**(0.5):
        spin = 1.
    else:
        spin = np.nan

    return spin


def check_state_of_star(star, i=None, star_CO=False):
    """Get the state of a SingleStar.

    Arguments
    ----------
    star: SingleStar
        The star for which the state will be computed.
    i : integer
        Index of the model for which we want to calculate the state of the
        SingleStar object. Default = -1 (the final).
    star_CO: bool
        True if we want to assume a compact object (WD/NS/BH).
        False if the state will be calculated by its attributes.

    Returns
    -------
    state : str
        state at of the star object at index i of the history.

    """
    if star_CO:
        star_mass = star.mass_history[i] if i is not None else star.mass
        if 'WD' == star.state_history[-1] or star.state == 'WD':
            return 'WD'
        else:
            return infer_star_state(star_mass=star_mass, star_CO=True)

    if i is None:
        star_mass = star.mass
        surface_h1 = star.surface_h1
        center_h1 = star.center_h1
        center_he4 = star.center_he4
        center_c12 = star.center_c12
        log_LH = star.log_LH
        log_LHe = star.log_LHe
        log_Lnuc = star.log_Lnuc
    else:
        star_mass = star.mass_history[i]
        surface_h1 = star.surface_h1_history[i]
        center_h1 = star.center_h1_history[i]
        center_he4 = star.center_he4_history[i]
        center_c12 = star.center_c12_history[i]
        log_LH = star.log_LH_history[i]
        log_LHe = star.log_LHe_history[i]
        log_Lnuc = star.log_Lnuc_history[i]

    return infer_star_state(star_mass=star_mass,
                            surface_h1=surface_h1,
                            center_h1=center_h1,
                            center_he4=center_he4,
                            center_c12=center_c12,
                            log_LH=log_LH,
                            log_LHe=log_LHe,
                            log_Lnuc=log_Lnuc,
                            star_CO=False)


def check_state_of_star_history_array(star, N=1, star_CO=False):
    """Calculate the evolutionary states with an array of history data.

    Parameters
    ----------
    star : SingleStar
        The star for which the state will be computed.
    N : int
        Number of history steps (from the end), to calculate the state.
    star_CO : bool
        If `True`, assume it's a compact object.

    Returns
    -------
    list of str
        The state(s) in a list.

    """
    return [check_state_of_star(star, i, star_CO) for i in range(-N, 0)]


def get_binary_state_and_event_and_mt_case(binary, interpolation_class=None,
                                           i=None, verbose=False):
    """Infer the binary state, event and mass transfer case.

    Parameters
    ----------
    binary : BinaryStar
        The POSYDON binary star.
    i : int
        The index of the history in which we want to calculate the state of the
        binary object. Default (None) is the current state.

    Returns
    -------
    binary_state : str
        One of 'detached', 'contact', 'RLO1' and 'RLO2'.

    binary_event : str
        Options are 'oRLO1' or 'oRLO2' (onset of RLO, the start of RLO),
        'oCE1', 'oCE2', 'oDoubleCE1', 'oDoubleCE2', 'CC1', 'CC2'.

    mass transfer case : str
        'caseA', 'caseB', etc.

    Examples
    --------
    If 'detached' then returns ['detached', None, None].

    If 'contact' then returns ['contact', None] or ['oCE1', 'oDoubleCE1'] or
    ['oDoubleCE2',  None].

    If RLO then returns either ['RLO1',  None, 'caseXX'] or
    ['RLO2',  None, 'caseXX'] or maybe ['RLO2',  'oRLO2', 'caseXX'].

    """
    # initializing: ['binary_state','binary_event','MT_case']
    if binary is None:
        return [None, None, 'None']

    if interpolation_class == 'not_converged':
        return [None, None, 'None']
    elif interpolation_class == 'initial_MT':
        return ['initial_RLOF', None, 'None']

    if i is None:
        lg_mtransfer = binary.lg_mtransfer_rate
        rl_overflow1 = binary.rl_relative_overflow_1
        rl_overflow2 = binary.rl_relative_overflow_2
        state1, state2 = binary.star_1.state, binary.star_2.state
        gamma1, gamma2 = binary.star_1.center_gamma, binary.star_2.center_gamma
    else:
        lg_mtransfer = binary.lg_mtransfer_rate_history[i]
        rl_overflow1 = binary.rl_relative_overflow_1_history[i]
        rl_overflow2 = binary.rl_relative_overflow_2_history[i]
        state1 = binary.star_1.state_history[i]
        state2 = binary.star_2.state_history[i]
        try:
            gamma1 = binary.star_1.center_gamma_history[i]
        except IndexError:  # this happens if compact object
            gamma1 = None
        try:
            gamma2 = binary.star_2.center_gamma_history[i]
        except IndexError:  # this happens if compact object
            gamma2 = None
        # prev_state currently unused, causing errors
        # prev_state = binary.state_history[i-1] if i != 0 else MT_CASE_NO_RLO

    # get numerical MT cases
    mt_flag_1 = infer_mass_transfer_case(rl_overflow1, lg_mtransfer, state1,
                                         verbose=verbose)
    mt_flag_2 = infer_mass_transfer_case(rl_overflow2, lg_mtransfer, state2,
                                         verbose=verbose)
    # convert to strings
    mt_flag_1_str = cumulative_mass_transfer_string([mt_flag_1])
    mt_flag_2_str = cumulative_mass_transfer_string([mt_flag_2])

    rlof1 = mt_flag_1 in ALL_RLO_CASES
    rlof2 = mt_flag_2 in ALL_RLO_CASES
    no_rlof = (mt_flag_1 == MT_CASE_NO_RLO) and (mt_flag_2 == MT_CASE_NO_RLO)

    if rlof1 and rlof2:                             # contact condition
        result = ['contact', None, 'None']
        if interpolation_class == 'unstable_MT':
            result = ['contact', 'oCE1', 'None']
    elif no_rlof:                                   # no MT in any star
        result = ['detached', None, 'None']
    elif rlof1 and not rlof2:                       # only in star 1
        result = ['RLO1', None, mt_flag_1_str]
        if interpolation_class == 'unstable_MT':
            return ['RLO1', 'oCE1', mt_flag_1_str]
        # if prev_state not in ALL_RLO_CASES:
        #    return ['RLO1', 'oRLO1', mt_flag_1_str]
    elif rlof2 and not rlof1:                       # only in star 2
        result = ['RLO2', None, mt_flag_2_str]
        if interpolation_class == 'unstable_MT':
            return ['RLO2', 'oCE2', mt_flag_2_str]
    else:                                           # undetermined in any star
        result = ["undefined", None, 'None']

    if result[1] == "oCE1":
        # Check for double CE
        comp_star = binary.star_2
        if comp_star.state not in [
                "H-rich_Core_H_burning",
                "stripped_He_Core_He_burning", "WD", "NS", "BH"]:
            result[1] = "oDoubleCE1"
    elif result[1] == "oCE2":
        # Check for double CE
        comp_star = binary.star_1
        if comp_star.state not in [
                "H-rich_Core_H_burning",
                "stripped_He_Core_He_burning", "WD", "NS", "BH"]:
            result[1] = "oDoubleCE2"

    if ("Central_C_depletion" in state1
            or "Central_He_depleted" in state1
            or (gamma1 is not None and gamma1 >= 10.0)):    # WD formation
        result[1] = "CC1"
    elif ("Central_C_depletion" in state2
          or "Central_He_depleted" in state2
          or (gamma2 is not None and gamma2 >= 10.0)):      # WD formation
        result[1] = "CC2"

    return result


def get_binary_state_and_event_and_mt_case_array(binary, N=None,
                                                 verbose=False):
    """Calculate the evolutionary states with an array of history data.

    Arguments
    ----------
    binary: POSYDON BinaryStar object
    N : array
        index of the history in which we want to calculate the state of the
        SingleStar object

    Returns
    -------
    binary_state : str(see in check_state_of_star)
    binary_event :
    MT_case :

    """
    if N is None:               # focus on current evolutonary state
        result = get_binary_state_and_event_and_mt_case(binary, i=None,
                                                        verbose=verbose)
        binary_state = result[0]
        binary_event = result[1]
        MT_case = result[2]
    else:
        binary_state = []
        binary_event = []
        MT_case = []
        for index in range(-N, 0):
            result = get_binary_state_and_event_and_mt_case(binary, i=index,
                                                            verbose=verbose)
            binary_state.append(result[0])
            binary_event.append(result[1])
            MT_case.append(result[2])

    return binary_state, binary_event, MT_case


def CO_radius(M, COtype):
    """Calculate the radius of a compact object based on its type and mass.

    Parameters
    ----------
    M : float
        CO mass in Msun
    COtype : str
        Tyep of compact object. Accepted values: "BH", "NS", "WD"

    Returns
    -------
    float
        Compact object radius in solar radii

    """
    if M <= 0.:
        raise ValueError('Compact object mass must be a positive value')

    if COtype == "BH":
        R = Schwarzschild_Radius(M)
    elif COtype == "NS":
        # Some references:
        # Most, E. R., Weih, L. R., Rezzolla, L., & Schaffner-Bielich, J. 2018,
        #   PhRvL, 120, 261103
        # Abbott, B. P., Abbott, R., Abbott, T. D., et al. 2020, ApJL, 892, L3
        # Landry, P., Essick, R., & Chatziioannou, K. 2020, PhRvD, 101, 123007
        R = 12.5e5/const.Rsun
    elif COtype == "WD":
        # Hansen C. J., Kawaler S. D., Trimble V., 2004,
        #   Stellar Interiors. Springer New York
        R = 2.9e8*(M)**(-1./3.)/const.Rsun
    else:
        raise ValueError(
            'COtype not in the list of valid options: "BH", "NS", "WD"')

    return R


def He_MS_lifetime(mass):
    """Calculate the lifetime of He burning in a star.

    Parameters
    ----------
    mass : type
        Mass of star in solar masses

    Returns
    -------
    float
        He MS time duration in yr.

    """
    if mass < 2.0:
        he_t_ms = 10 ** 8
    elif mass >= 2.0 and mass < 10.0:
        he_t_ms = 10**(-2.6094 * np.log10(mass) + 8.7855)
    elif mass >= 10.0 and mass < 100.0:
        he_t_ms = 10**(-0.69897 * np.log10(mass) + 6.875)
    elif mass >= 100.0:
        he_t_ms = 3 * 10 ** 5
    return he_t_ms


def Schwarzschild_Radius(M):
    """Calculate the Schwarzschild Radius of BH with mass M.

    Parameters
    ----------
    M : float
        BH mass in Msun

    Returns
    -------
    float
        Schwarzschild Radius in solar radii

    """
    G = const.standard_cgrav
    c = const.clight
    Rsun = const.Rsun
    Msun = const.Msun

    # Kutner, M. L. 2003, Astronomy: A physical perspective,
    #   Cambridge University Press
    return (2 * G * M * Msun / (c**2 * Rsun))


def flip_stars(binary):
    """Short summary.

    Parameters
    ----------
    binary : type
        Description of parameter `binary`.

    Returns
    -------
    type
        Description of returned object.

    """
    star_1 = copy.copy(binary.star_1)
    star_2 = copy.copy(binary.star_2)
    setattr(binary, 'star_1', star_2)
    setattr(binary, 'star_2', star_1)

    state = getattr(binary, 'state')
    if state == 'RLO1':
        setattr(binary, 'state', 'RLO2')
    elif state == 'RLO2':
        setattr(binary, 'state', 'RLO1')
    event = getattr(binary, 'event')
    if event == 'oRLO1':
        setattr(binary, 'event', 'oRLO2')
    elif event == 'oRLO2':
        setattr(binary, 'event', 'oRLO1')
    if event == 'oCE1':
        setattr(binary, 'event', 'oCE2')
    elif event == 'oCE2':
        setattr(binary, 'event', 'oCE1')
    if event == 'CC1':
        setattr(binary, 'event', 'CC2')
    elif event == 'CC2':
        setattr(binary, 'event', 'CC1')

    state_history = np.array(getattr(binary, 'state_history'))
    cond_RLO2 = state_history == 'RLO1'
    cond_RLO1 = state_history == 'RLO2'
    state_history[cond_RLO2] = 'RLO2'
    state_history[cond_RLO1] = 'RLO1'
    setattr(binary, 'state_history', state_history.tolist())

    event_history = np.array(getattr(binary, 'event_history'))
    cond_CC2 = event_history == 'CC1'
    cond_CC1 = event_history == 'CC2'
    event_history[cond_CC2] = 'CC2'
    event_history[cond_CC1] = 'CC1'
    cond_oRLO2 = event_history == 'oRLO1'
    cond_oRLO1 = event_history == 'oRLO2'
    event_history[cond_oRLO2] = 'oRLO2'
    event_history[cond_oRLO1] = 'oRLO1'
    cond_oCE2 = event_history == 'oCE1'
    cond_oCE1 = event_history == 'oCE2'
    event_history[cond_oCE2] = 'oCE2'
    event_history[cond_oCE1] = 'oCE1'
    setattr(binary, 'event_history', event_history.tolist())

    for i in ['t_sync_rad_', 't_sync_conv_', 'rl_relative_overflow_']:

        value1 = getattr(binary, i+'1')
        value2 = getattr(binary, i+'2')
        value1_history = getattr(binary, i+'1_history')
        value2_history = getattr(binary, i+'2_history')

        setattr(binary, i+'1', value2)
        setattr(binary, i+'2', value1)
        setattr(binary, i+'1_history', value2_history)
        setattr(binary, i+'2_history', value1_history)


def infer_star_state(star_mass=None, surface_h1=None,
                     center_h1=None, center_he4=None, center_c12=None,
                     log_LH=None, log_LHe=None, log_Lnuc=None, star_CO=False):
    """Infer the star state (corresponding to termination flags 2 and 3)."""
    if star_CO:
        return "NS" if star_mass <= STATE_NS_STARMASS_UPPER_LIMIT else "BH"

    if surface_h1 is None:
        return STATE_UNDETERMINED

    rich_in = ("H-rich" if surface_h1 > THRESHOLD_HE_NAKED_ABUNDANCE
               else "stripped_He")
    burning_H = (log_LH > LOG10_BURNING_THRESHOLD
                 and log_LH - log_Lnuc > REL_LOG10_BURNING_THRESHOLD)
    burning_He = (log_LHe > LOG10_BURNING_THRESHOLD
                  and log_LHe - log_Lnuc > REL_LOG10_BURNING_THRESHOLD)

    H_in_core = center_h1 > THRESHOLD_CENTRAL_ABUNDANCE
    He_in_core = center_he4 > THRESHOLD_CENTRAL_ABUNDANCE
    C_in_core = center_c12 > THRESHOLD_CENTRAL_ABUNDANCE

    if not (H_in_core or He_in_core):   # H and He are depleted
        if not C_in_core:
            burning = "Central_C_depletion"
        else:
            burning = "Central_He_depleted"
        # from now on, either H or He in core
    elif H_in_core:                     # no matter what the He abundance is
        if burning_H:
            burning = "Core_H_burning"
        else:
            burning = "non_burning"
    else:                               # from now on: only He, not H in core
        if burning_He:
            burning = "Core_He_burning"
        elif burning_H:
            burning = "Shell_H_burning"
        else:
            burning = "non_burning"

    return "{}_{}".format(rich_in, burning)


def infer_mass_transfer_case(rl_relative_overflow,
                             lg_mtransfer_rate,
                             donor_state,
                             verbose=False):
    """Infer the mass-transfer case of a given star.

    Parameters
    ----------
    rl_relative_overflow : float
    lg_mtransfer_rate : float
    donor_state : str
        Values of star parameters at a specific step.

    Returns
    -------
    int
        The mass-transfer case integer flag.

    """
    if rl_relative_overflow is None or lg_mtransfer_rate is None:
        return MT_CASE_NO_RLO

    if (rl_relative_overflow <= RL_RELATIVE_OVERFLOW_THRESHOLD
            or (lg_mtransfer_rate <= LG_MTRANSFER_RATE_THRESHOLD
                and rl_relative_overflow < 0.0)):
        if verbose:
            print("checking rl_relative_overflow / lg_mtransfer_rate,",
                  rl_relative_overflow, lg_mtransfer_rate)
        return MT_CASE_NO_RLO

    if "non_burning" in donor_state:
        return MT_CASE_NONBURNING
    elif "H-rich" in donor_state:
        if "Core_H_burning" in donor_state:
            return MT_CASE_A
        if ("Core_He_burning" in donor_state
                or "Shell_H_burning" in donor_state):
            return MT_CASE_B
        if ("Central_He_depleted" in donor_state
                or "Central_C_depletion" in donor_state):
            return MT_CASE_C
    elif "stripped_He" in donor_state:
        if "Core_He_burning" in donor_state:
            return MT_CASE_BA
        if ("Central_He_depleted" in donor_state
                or "Central_C_depletion" in donor_state):
            return MT_CASE_BB
    return MT_CASE_UNDETERMINED


def cumulative_mass_transfer_numeric(MT_cases):
    """Summarize the history of MT cases in a short list of integers.

    Parameters
    ----------
    MT_cases : array-like
        A list of the integer MT flags at sequential history steps.

    Returns
    -------
    list of int
        A shorter list of integer MT flags, following these rules:
        i.   If undetermined MT at any step, it is indicated in the beginning
             (as a warning) and ignored later.
        ii.  If no RLO anywhere, then this is the only element in the returned
             list (or the second if undetermined MT at any step). Intermediate
             no-RLO phases (between RLO cases) are ignored.
        iii. Consequent same cases are reported only once (e.g., A, A, A -> A.)

    The end result is a list of integers indicating whether undetermined MT
    anywhere, if no RLO everywhere, or "changes" of MT cases.

    """
    if len(MT_cases) == 0:
        return [MT_CASE_UNDETERMINED]

    if isinstance(MT_cases[0], str):
        cases = [MT_STR_TO_CASE[case] for case in MT_cases]
    else:
        cases = MT_cases.copy()

    result = []

    # if undetermined MT anywhere, report it at the beginning and forget it
    if MT_CASE_UNDETERMINED in cases:
        result.append(MT_CASE_UNDETERMINED)
        cases = [case for case in cases if case != MT_CASE_UNDETERMINED]
        if len(cases) == 0:
            return result

    # if no RLO at all steps, report this, otherwise forget about these phases
    if MT_CASE_NO_RLO in cases:
        cases = [case for case in cases if case != MT_CASE_NO_RLO]
        if len(cases) == 0:
            result.append(MT_CASE_NO_RLO)
            return result

    # from now on... undetermined, and no_RLO are not in the list...
    curr_case = cases[0]
    result.append(curr_case)
    prev_case = curr_case
    for curr_case in cases[1:]:
        if curr_case != prev_case:
            result.append(curr_case)
            prev_case = curr_case

    return result


def cumulative_mass_transfer_string(cumulative_integers):
    """Convert a cumulative MT sequence to a string.

    Parameters
    ----------
    cumulative_integers : list of int
        Typically, the output of `cumulative_mass_transfer_numeric`.

    Returns
    -------
    str
        A summarization of the mass-tranfer cases in the form of a string, as
        opposed to the output of `cumulative_mass_transfer_numeric` (see this
        function to understand the rules of summarization). The character `?`
        in the beginning of the string indicates undetermined MT case at some
        point of the evolution. `no_RLO` indicates no RLO at any evolutionary
        step. `caseA`, `caseB`, etc., denote the corresponding MT cases. When
        multiple cases are found, they are indicated only when the begin, and
        are combined using `/`. For example:

        ?no_RLO     : undetermined MT at few steps, otherwise no RLO.
        caseA       : case A MT only (possible no_RLO at some points).
        ?caseA/B    : undetermined MT somewhere, then case A, then case B.
        caseA/B/A   : case A, then B, and A again (although unphysical).

    """
    assert len(cumulative_integers) != 0
    result = ""
    added_case_word = False
    for integer in cumulative_integers:
        if integer == MT_CASE_UNDETERMINED:
            result += "?"
        elif integer == MT_CASE_NO_RLO:
            result += "no_RLO"
        else:
            if not added_case_word:
                result += "case_"
                added_case_word = True
            else:
                result += "/"
            if integer in MT_CASE_TO_STR:
                result += MT_CASE_TO_STR[integer] + '1' # from star 1
            else:
                result += MT_CASE_TO_STR[integer-10] + '2' # from star 2
    return result


def cumulative_mass_transfer_flag(MT_cases):
    """Get the cumulative MT string from a list of integer MT casses."""
    return cumulative_mass_transfer_string(
        cumulative_mass_transfer_numeric(MT_cases)
    )


def calculate_Patton20_values_at_He_depl(star):
    """Calculate the carbon core mass and abundance very close to ignition.

    This is important for using the Patton+2020 SN prescription

    Arguments
    ----------
    star: SingleStar object holding the history of its attributes.

    Returns
    -------
    None

    It updates the following values in the star object
    co_core_mass_at_He_depletion: float
        co_core_mass at He core depletion
        (almost at the same time as carbon core ignition)
    avg_c_in_c_core_at_He_depletion : float
        avg carbon abundance inside CO_core_mass at He core depletion
        (almost at the same time as carbon core ignition)

    """
    co_core_mass_at_He_depletion = None
    avg_c_in_c_core_at_He_depletion = None
    if star.state_history is not None:
        if ("H-rich_Central_He_depleted" in star.state_history):
            i_He_depl = np.argmax(
                np.array(star.state_history) == "H-rich_Central_He_depleted")
            co_core_mass_at_He_depletion = star.co_core_mass_history[i_He_depl]
            avg_c_in_c_core_at_He_depletion = star.avg_c_in_c_core_history[
                i_He_depl]
        elif ("stripped_He_Central_He_depleted" in star.state_history):
            i_He_depl = np.argmax(np.array(star.state_history)
                                  == "stripped_He_Central_He_depleted")
            co_core_mass_at_He_depletion = star.co_core_mass_history[i_He_depl]
            avg_c_in_c_core_at_He_depletion = star.avg_c_in_c_core_history[
                i_He_depl]
    else:
        co_core_mass_at_He_depletion = None
        avg_c_in_c_core_at_He_depletion = None

    # return co_core_mass_at_He_depletion, avg_c_in_c_core_at_He_depletion
    star.co_core_mass_at_He_depletion = co_core_mass_at_He_depletion
    star.avg_c_in_c_core_at_He_depletion = avg_c_in_c_core_at_He_depletion


def CEE_parameters_from_core_abundance_thresholds(star, verbose=False):
    """Find the envelope mass for different core boundary abundance thresholds.

    The results are meant to be used in collabration with the respective
    `lambda_CE_*cent`, `lambda_CE_pure_He_star_10cent`.

    Arguments
    ----------
    star: SingleStar object holding the history of its attributes.

    Returns
    -------
    None

    It updates the following values in the star object
    co_core_mass_at_He_depletion: float
        co_core_mass at He core depletion
        (almost at the same time as carbon core ignition).
    avg_c_in_c_core_at_He_depletion : float
        avg carbon abundance inside CO_core_mass at He core depletion
        (almost at the same time as carbon core ignition).

    """
    mass = star.mass
    radius = 10.**star.log_R
    m_core_CE_1cent = 0.0
    m_core_CE_10cent = 0.0
    m_core_CE_30cent = 0.0
    m_core_CE_pure_He_star_10cent = 0.0
    r_core_CE_1cent = 0.0
    r_core_CE_10cent = 0.0
    r_core_CE_30cent = 0.0
    r_core_CE_pure_He_star_10cent = 0.0
    profile = star.profile              # final profile of a star in a MESA run

    if profile is not None and isinstance(profile, np.ndarray):
        mass_prof = profile["mass"]
        star_state = star.state

        m_core = 0.0
        r_core = 0.0

        if "H-rich" in star_state:
            for element_frac in [0.01, 0.1, 0.3]:
                ind_core = calculate_core_boundary(
                    mass_prof, star_state, profile,
                    core_element_fraction_definition=element_frac)
                lambda_CE, m_core, r_core = calculate_lambda_from_profile(
                    profile=profile, donor_star_state=star_state,
                    m1_i=mass, radius1=radius,
                    core_element_fraction_definition=element_frac,
                    ind_core=ind_core, verbose=verbose)

                if element_frac == 0.01:
                    m_core_CE_1cent = m_core
                    r_core_CE_1cent = r_core
                    lambda_CE_1cent = lambda_CE
                elif element_frac == 0.1:
                    m_core_CE_10cent = m_core
                    r_core_CE_10cent = r_core
                    lambda_CE_10cent = lambda_CE
                if element_frac == 0.3:
                    m_core_CE_30cent = m_core
                    r_core_CE_30cent = r_core
                    lambda_CE_30cent = lambda_CE

            # calculate also potential CO core, for consistency
            for element_frac in [0.1]:
                ind_core = calculate_core_boundary(
                    mass_prof, star_state, profile,
                    core_element_fraction_definition=element_frac,
                    CO_core_in_Hrich_star=True)
                lambda_CE, m_core, r_core = calculate_lambda_from_profile(
                    profile=profile, donor_star_state=star_state,
                    m1_i=mass, radius1=radius,
                    core_element_fraction_definition=element_frac,
                    ind_core=ind_core, CO_core_in_Hrich_star=True,
                    verbose=verbose)
                m_core_CE_pure_He_star_10cent = m_core
                r_core_CE_pure_He_star_10cent = r_core
                lambda_CE_pure_He_star_10cent = lambda_CE
        elif "stripped_He" in star_state:
            for element_frac in [0.1]:
                ind_core = calculate_core_boundary(
                    mass_prof, star_state, profile,
                    core_element_fraction_definition=element_frac)
                lambda_CE, m_core, r_core = calculate_lambda_from_profile(
                    profile=profile, donor_star_state=star_state,
                    m1_i=mass, radius1=radius,
                    core_element_fraction_definition=element_frac,
                    ind_core=ind_core,
                    verbose=verbose)
                m_core_CE_pure_He_star_10cent = m_core
                m_core_CE_1cent = mass
                m_core_CE_10cent = mass
                m_core_CE_30cent = mass
                r_core_CE_pure_He_star_10cent = r_core
                r_core_CE_1cent = radius
                r_core_CE_10cent = radius
                r_core_CE_30cent = radius
                lambda_CE_pure_He_star_10cent = lambda_CE
                lambda_CE_1cent = 1e99
                lambda_CE_10cent = 1e99
                # TODO: decide whether this should be None
                #       for interpolation reasons
                lambda_CE_30cent = 1e99
        else:   # CO-object or undetermined_evolutionary_state?
            m_core_CE_pure_He_star_10cent = np.nan
            m_core_CE_1cent = np.nan
            m_core_CE_10cent = np.nan
            m_core_CE_30cent = np.nan
            r_core_CE_pure_He_star_10cent = np.nan
            r_core_CE_1cent = np.nan
            r_core_CE_10cent = np.nan
            r_core_CE_30cent = np.nan
            lambda_CE_1cent = np.nan
            lambda_CE_10cent = np.nan
            lambda_CE_30cent = np.nan
            lambda_CE_pure_He_star_10cent = np.nan
            if verbose:
                print('star state {} is not what expected during the '
                      'CEE_parameters_from_core_abundance_thresholds.'.
                      format(star_state))

    else:                                                       # no profile
        m_core_CE_pure_He_star_10cent = np.nan
        m_core_CE_1cent = np.nan
        m_core_CE_10cent = np.nan
        m_core_CE_30cent = np.nan
        r_core_CE_pure_He_star_10cent = np.nan
        r_core_CE_1cent = np.nan
        r_core_CE_10cent = np.nan
        r_core_CE_30cent = np.nan
        lambda_CE_1cent = np.nan
        lambda_CE_10cent = np.nan
        lambda_CE_30cent = np.nan
        lambda_CE_pure_He_star_10cent = np.nan

    if verbose:
        print("star_state", star_state)
        print("m_core_CE_1cent,m_core_CE_10cent,m_core_CE_30cent,"
              "m_core_CE_pure_He_star_10cent",
              m_core_CE_1cent, m_core_CE_10cent, m_core_CE_30cent,
              m_core_CE_pure_He_star_10cent)
        print("r_core_CE_1cent,r_core_CE_10cent,r_core_CE_30cent,"
              "r_core_CE_pure_He_star_10cent",
              r_core_CE_1cent, r_core_CE_10cent, r_core_CE_30cent,
              r_core_CE_pure_He_star_10cent)
        print("lambda_CE_1cent,lambda_CE_10cent,lambda_CE_30cent,"
              "lambda_CE_pure_He_star_10cent",
              lambda_CE_1cent, lambda_CE_10cent, lambda_CE_30cent,
              lambda_CE_pure_He_star_10cent)

    star.m_core_CE_1cent = m_core_CE_1cent
    star.m_core_CE_10cent = m_core_CE_10cent
    star.m_core_CE_30cent = m_core_CE_30cent
    star.m_core_CE_pure_He_star_10cent = m_core_CE_pure_He_star_10cent
    star.r_core_CE_1cent = r_core_CE_1cent
    star.r_core_CE_10cent = r_core_CE_10cent
    star.r_core_CE_30cent = r_core_CE_30cent
    star.r_core_CE_pure_He_star_10cent = r_core_CE_pure_He_star_10cent
    star.lambda_CE_1cent = lambda_CE_1cent
    star.lambda_CE_10cent = lambda_CE_10cent
    star.lambda_CE_30cent = lambda_CE_30cent
    star.lambda_CE_pure_He_star_10cent = lambda_CE_pure_He_star_10cent


def initialize_empty_array(arr):
    """Initialize an empty record array with NaNs and empty strings."""
    res = arr.copy()
    for colname in res.dtype.names:
        if np.issubsctype(res[colname], float):
            res[colname] = np.nan
        if np.issubsctype(res[colname], str):
            res[colname] = np.nan
    return res


def calculate_core_boundary(donor_mass,
                            donor_star_state,
                            profile,
                            mc1_i=None,
                            core_element_fraction_definition=None,
                            CO_core_in_Hrich_star=False):
    """Calculate the shell where the core is - envelope boundary.

    Parameters
    ----------
    donor_mass : array
        Profile column of enclosed mass of the star.
    donor_star_state : string
        The POSYDON evolutionary state of the donor star
    profile : numpy.array
        Donor's star profile from MESA
    mc1_i : float
        core mass and total mass of the donor star.
    core_element_fraction_definition : float
        The mass fraction of the envelope abundant chemical element at
        the core-envelope boundary to derive the donor's core mass from
        the profile
    CO_core_in_Hrich_star: Bool
        This should be true if we want to calculate the boundary of CO core in
        a H-rich star (and not of the helium core).

    Returns
    -------
    ind_core : int
        The value of the cell position of the core - envelope boundary, at the
        donor's MESA  profile (for a profile that starts from the surface).
        More specifically, it returns the index of the first cell (starting
        from the surface), that the elements conditions for the core are met.
        - Returns 0 for a star that is all core
        - Returns -1 for a star that is all envelope.

    """
    # the threshold of the elements that need to be high in the core
    core_element_high_fraction_definition = 0.1
    # ENHANCEMENT: this list needs to be imported from e.g. flow_chart.py
    STAR_STATES_H_RICH = [
        "H-rich_Core_H_burning",
        "H-rich_Shell_H_burning",
        "H-rich_Core_He_burning",
        "H-rich_Central_He_depleted",
        "H-rich_Core_C_burning",
        "H-rich_Central_C_depletion",
        "H-rich_non_burning",
    ]
    # ENHANCEMENT: this list needs to be imported from e.g. flow_chart.py
    STAR_STATE_He = [
        'stripped_He_Core_He_burning',
        'stripped_He_Central_He_depleted',
        'stripped_He_Central_C_depletion',
        'stripped_He_non_burning'
    ]

    if core_element_fraction_definition is not None:
        if ((donor_star_state in STAR_STATES_H_RICH)
                and ('x_mass_fraction_H' in profile.dtype.names)
                and ('y_mass_fraction_He' in profile.dtype.names)
                and ('z_mass_fraction_metals' in profile.dtype.names)):
            if not CO_core_in_Hrich_star:
                element = profile['x_mass_fraction_H']
                element_core = np.add(profile['y_mass_fraction_He'],
                                      profile['z_mass_fraction_metals'])
            else:
                element = np.add(profile['x_mass_fraction_H'],
                                 profile['y_mass_fraction_He'])
                element_core = profile['z_mass_fraction_metals']
        elif (donor_star_state in STAR_STATE_He
              and 'x_mass_fraction_H' in profile.dtype.names
              and 'y_mass_fraction_He' in profile.dtype.names
              and 'z_mass_fraction_metals' in profile.dtype.names):
            # Recalculate the core from a chemical mass fraction threshold.
            # Here we assume Xelement=0.1 is the default MESA core definition.
            # element = profile['y_mass_fraction_He']
            element = np.add(profile['x_mass_fraction_H'],
                             profile['y_mass_fraction_He'])
            element_core = profile['z_mass_fraction_metals']
        # ind_core=np.argmax(element[::-1]>=core_element_fraction_definition)
        else:
            ind_core = -1
            warnings.warn("Profile columns not enough to calculate the core "
                          "boundaries for CE, all star considered an envelope")
            return ind_core

        # starting from the surface, both conditions become True when element
        # (of which the envelope is rich) decreases and the element_core which
        # is in the core increases.
        both_conditions = (
            element <= core_element_fraction_definition).__and__(
                element_core >= core_element_high_fraction_definition)
        if not np.any(both_conditions):
            # the whole star is an envelope, from surface towards the core
            # the "both_conditions" never becomes True
            ind_core = -1
        else:
            # starting from the surface we find the first time that we get True
            ind_core = np.argmax(both_conditions)
            # This includes the case that the whole star is an "core", as it
            # will find only "True" in both_conditions and will return
            # ind_core=0 (first, surface cell for a MESA profile)
    elif (mc1_i is not None):
        # calculate the cell position of the core boundary. In principle you
        # calculate index of the first time the inequality becomes False (your
        # lowest value) for your MESA profile, so starting from the surface.
        ind_core = np.argmin(donor_mass >= mc1_i)
    else:
        raise ValueError(
            "Not possible to calculate the core boundary of the donor in CE")
    return ind_core


def period_evol_wind_loss(M_current, M_init, Mcomp, P_init):
    """Calculate analytically the period widening due to wind mass loss [1]_.

    Parameters
    ----------
    M_current : float
        Current mass of the mass losing star (Msun)
    M_init : float
        Initial mass of the mass losing star when the calculation starts (Msun)
    Mcomp : float
        (Constant) mass of the companion star (Msun)
    P_init : float
        Initial binary period when the calculation starts (days)

    Returns
    -------
     float
        Current binary period  in days
    References
    ----------
    .. [1] Tauris, T. M., & van den Heuvel, E. 2006, Compact stellar X-ray
        sources, 1, 623

    """
    log10P = (-2.*np.log10(M_current+Mcomp)
              + 2.*np.log10(M_init+Mcomp) + np.log10(P_init))
    return 10.0**log10P


def separation_evol_wind_loss(M_current, M_init, Mcomp, A_init):
    """Calculate analytically the separation widening due to wind mass loss [1]_.

    Parameters
    ----------
    M_current : float
        Current mass of the mass losing star (Msun)
    M_init : float
        Initial mass of the mass losing star when the calculation starts (Msun)
    Mcomp : float
        (Constant) mass of the companion star (Msun)
    A_init : float
        Initial binary separation when the calculation starts (Rsun)

    Returns
    -------
     float
        Current binary separation  in Rsun
    References
    ----------
    .. [1] Tauris, T. M., & van den Heuvel, E. 2006, Compact stellar X-ray
        sources, 1, 623.

    """
    log10A = (-np.log10(M_current+Mcomp)
              + np.log10(M_init+Mcomp) + np.log10(A_init))
    return 10.0**log10A


def period_change_stabe_MT(period_i, Mdon_i, Mdon_f, Macc_i,
                           alpha=0.0, beta=0.0):
    """Change the binary period after a semi-detahed stable MT phase.

    Calculated in Sorensen, Fragos et al.  2017A&A...597A..12S.
    Note that MT efficiencies are assumed constant (i.e., not time-dependent)
    throughout the MT phase.

    Parameters
    ----------
    period_i : float
        Initial period
    Mdon_i: float
        Initial donor mass
    Mdon_f: float
        Final donor mass (should be in the same unit's of Mdon_i)
    Mdon_i: float
        Initial accretor's mass (should be in the same unit's of Mdon_i)
    alpha : float [0-1]
        Fraction of DM_don (= Mdon_i - Mdon_f) from the donor,
        lost from the donor's vicinity
    beta : float [0-1]
        Fraction of Mdot from the L1 point (= (1-alpha)*DM_don),
        lost from the accretor's vicinity.
        The final accreted rate is (1-beta)(1-alpha)*DM_don

    Returns
    -------
    period_f : float
        final period at the end of stable MT, in the same units as period_i

    """
    DM_don = Mdon_i - Mdon_f    # mass lost from donor (>0)
    Macc_f = Macc_i + (1.-beta)*(1.-alpha)*DM_don
    if alpha < 0.0 or beta < 0.0 or alpha > 1.0 or beta > 1.0:
        raise ValueError("In period_change_stabe_MT, mass transfer "
                         "efficiencies, alpha, beta {}{} are not in the [0-1] "
                         "range.".format(alpha, beta))
    if beta != 1.0:      # Eq. 7 of Sorensen+Fragos et al. 2017
        period_f = (period_i * (Mdon_f/Mdon_i)**(3.*(alpha-1.))
                    * (Macc_f/Macc_i)**(3./(beta-1.))
                    * ((Mdon_i + Macc_i)/(Mdon_f + Macc_f))**(2.))
    else:
        # fully non-conservative MT. Eq. 8 of Sorensen+Fragos et al. 2017,
        # were we already assumed beta=1
        period_f = (period_i * (Mdon_f/Mdon_i)**(3.*(alpha-1.))
                    * np.exp(3.*(1.-alpha)*(Mdon_f - Mdon_i)/Macc_f)
                    * ((Mdon_i + Macc_i)/(Mdon_f + Macc_f))**(2.))

    return period_f


def linear_interpolation_between_two_cells(array_y, array_x, x_target,
                                           top=None, bot=None, verbose=False):
    """Interpolate quantities between two star profile shells."""
    if ((np.isnan(top) or top is None) and (np.isnan(bot) or bot is None)):
        top = np.argmax(array_x >= x_target)
        bot = top - 1
    elif np.isnan(bot) or bot is None:
        bot = top - 1
    elif np.isnan(top) or top is None:
        top = bot + 1

    if top > len(array_x):
        y_target = array_y[top]
    if bot < 0:
        bot = 0

    if top == bot:
        y_target = array_y[top]
        warnings.warn("linear interpolation between the same point")
        if verbose:
            print("linear interpolation, but at the edge")
            print("x_target,top, bot, len(array_x), y_target",
                  x_target, top, bot, len(array_x), y_target)
    else:
        x_top = array_x[top]
        x_bot = array_x[bot]

        y_top = array_y[top]
        y_bot = array_y[bot]

        slope = (y_top - y_bot) / (x_top - x_bot)
        const = (y_top*x_bot - y_bot*x_top) / (x_top - x_bot)
        y_target = slope * x_target - const

        if verbose:
            print("linear interpolation")
            print("x_target,top, bot, len(array_x)",
                  x_target, top, bot, len(array_x))
            print("x_top, x_bot, y_top, y_bot, y_target",
                  x_top, x_bot, y_top, y_bot, y_target)

    return y_target


def calculate_lambda_from_profile(
        profile, donor_star_state,  m1_i=np.nan, radius1=np.nan,
        common_envelope_option_for_lambda=DEFAULT_CE_OPTION_FOR_LAMBDA,
        core_element_fraction_definition=0.1, ind_core=None,
        common_envelope_alpha_thermal=1.0, tolerance=0.001,
        CO_core_in_Hrich_star=False, verbose=False):
    """Calculate common-enevelope lambda from the profile of a star.

     We also pass a more accurate calculation of the donor core mass for the
     purposes of common-envelope evolution.

    Parameters
    ----------
    profile : numpy.array
        Donor's star profile from MESA
    donor_star_state : string
        The POSYDON evolutionary state of the donor star
    common_envelope_option_for_lambda : str
        Available options:
        * 'default_lambda': using for lambda the constant value of
        common_envelope_lambda_default parameter
        * 'lambda_from_profile_gravitational': calculating the lambda
        parameter from the donor's profile by using the gravitational
        binding energy from the surface to the core
        (needing "mass", and "radius" as columns in the profile)
        * 'lambda_from_profile_gravitational_plus_internal': as above
         but taking into account a factor of common_envelope_alpha_thermal *
         internal energy too in the binding energy (needing also "energy" as
         column in the profile)
        * 'lambda_from_profile_gravitational_plus_internal_minus_recombination'
        as above but not taking into account the recombination energy in the
        internal energy (needing also "y_mass_fraction_He", "x_mass_fraction_H"
        "neutral_fraction_H", "neutral_fraction_He", and "avg_charge_He" as
        column in the profile)
    core_element_fraction_definition : float
        The mass fraction of the envelope abundant chemical element at
        the core-envelope boundary to derive the donor's core mass from
        the profile.
    ind_core : int
        The value of the cell position of the core - envelope boundary, at the
        donor's MESA  profile (for a profile that starts from the surface).
        More specifically, it returns the index of the first cell (starting
        from the surface), that the elements conditions for the core are met.
        It is 0 for a star that is all core and -1 for a star that is all
        envelope. If it is not given (None), it will be calculated inside the
        function.
    common_envelope_alpha_thermal : float
        Used and explained depending on the common_envelope_option_for_lambda
        option above.
    tolerance : float
        The tolerance of numerical difference in two floats when comparing
        and testing results.
    CO_core_in_Hrich_star: Bool
        This should be true if we want to calculate the boundary of CO core in
        a H-rich star (and not of the helium core).
    verbose : bool
        In case we want information about the CEE.

    Returns
    -------
    lambda_CE: float
        lambda_CE for the envelope of the donor in CEE,
        calculated from profile
    mc1_i: float
        More accurate calculation of the donor core mass for the purposes
        of CEE.
    rc1_i: float
        More accurate calculation of the donor core radius for the purposes
        of CEE.

    """
    # get mass and radius and dm from profile
    donor_mass, donor_radius, donor_dm = get_mass_radius_dm_from_profile(
        profile, m1_i, radius1, tolerance)
    # if np.isnan(m1_i) or m1_i is None or np.isnan(radius1) or radius1 is None
    m1_i = donor_mass[0]
    radius1 = donor_radius[0]
    specific_internal_energy = get_internal_energy_from_profile(
        common_envelope_option_for_lambda, profile, tolerance)

    if ind_core is None:
        # To be used in a MESA profile (so one that starts from the surface)
        ind_core = calculate_core_boundary(donor_mass, donor_star_state,
                                           profile, None,
                                           core_element_fraction_definition,
                                           CO_core_in_Hrich_star)
    if ind_core == 0:       # all star is a core, immediate successful ejection
        Ebind_i = 0.0
        # TODO: decide whether this should be None for interpolation reasons
        lambda_CE = 1.0e99
        mc1_i = m1_i
        rc1_i = radius1
    else:
        if ind_core == -1:  # all star is an envelope, calculate for whole star
            Ebind_i = calculate_binding_energy(
                donor_mass, donor_radius, donor_dm, specific_internal_energy,
                len(donor_mass), common_envelope_alpha_thermal, verbose)
            mc1_i = 0.0
            rc1_i = 0.0
        else:
            if "H-rich" in donor_star_state:
                if not CO_core_in_Hrich_star:
                    elem_prof = profile["x_mass_fraction_H"]
                else:
                    elem_prof = profile["y_mass_fraction_He"]
            elif "stripped_He" in donor_star_state:
                elem_prof = profile["y_mass_fraction_He"]
            mc1_i = linear_interpolation_between_two_cells(
                donor_mass, elem_prof, core_element_fraction_definition,
                ind_core, ind_core-1, verbose)
            rc1_i = linear_interpolation_between_two_cells(
                donor_radius, elem_prof, core_element_fraction_definition,
                ind_core, ind_core-1, verbose)

            # linear interpolation
            Ebind_i_top = calculate_binding_energy(
                donor_mass, donor_radius, donor_dm, specific_internal_energy,
                ind_core, common_envelope_alpha_thermal, verbose)
            Ebind_i_bot = calculate_binding_energy(
                donor_mass, donor_radius, donor_dm, specific_internal_energy,
                ind_core-1, common_envelope_alpha_thermal, verbose)
            if verbose:
                lambda_CE_top = (-const.standard_cgrav * m1_i * const.Msun
                                 * (m1_i - donor_mass[ind_core]) * const.Msun
                                 / (Ebind_i_top*radius1*const.Rsun))
                lambda_CE_bot = (-const.standard_cgrav * m1_i * const.Msun
                                 * (m1_i - donor_mass[ind_core-1]) * const.Msun
                                 / (Ebind_i_bot*radius1*const.Rsun))
                print("lambda_CE_top, lambda_CE_bot",
                      lambda_CE_top, lambda_CE_bot)
            weight = ((core_element_fraction_definition-elem_prof[ind_core-1])
                      / (elem_prof[ind_core] - elem_prof[ind_core-1]))
            Ebind_i = Ebind_i_bot + weight*(Ebind_i_top - Ebind_i_bot)
        # lambda of the donor is calculated from the profile
        lambda_CE = (-const.standard_cgrav * m1_i*const.Msun
                     * (m1_i - mc1_i)*const.Msun/(Ebind_i*radius1*const.Rsun))
    if verbose:
        print("m1_i, radius1, len(profile) vs ind_core, mc1_i, rc1_i",
              m1_i, radius1, len(donor_mass), " vs ", ind_core, mc1_i, rc1_i)
        print("Ebind_i from profile ", Ebind_i)
        print("lambda_CE ", lambda_CE)
    if not (lambda_CE > -tolerance):
        raise Exception("CEE problem, lamda_CE has negative value.")
    return lambda_CE, mc1_i, rc1_i


def get_mass_radius_dm_from_profile(profile, m1_i=0.0,
                                    radius1=0.0, tolerance=0.001):
    """TODO: add summary.

    Reads and returns the profile columms of enclosed mass radius and
    mass per shell of the donor star from a MESA profile.

    Parameters
    ----------
    m1_i : float
        m1_i is the value passed from the singlestar object,
        for testing purposes only
    radius1 : float
        radius1 is the value passed from the singlestar object,
        for testing purposes only
    profile : numpy.array
        Donor's star profile from MESA
    tolerance : float
        The tolerance of numerical difference in two floats when comparing and
        testing results.

    Returns
    -------
    donor_mass : array
        Profile column of enclosed mass of the star.
    donor_radius : array
        Profile column of the radius the star.
    donor_dm : array
        Profile mass per shell of the star.

    """
    if (("mass" in profile.dtype.names)
            and (("radius" in profile.dtype.names)
                 or ("log_R" in profile.dtype.names))):
        donor_mass = profile["mass"]

        if ("radius" in profile.dtype.names):
            donor_radius = profile["radius"]
        elif ("log_R" in profile.dtype.names):
            donor_radius = 10**profile["log_R"]

        # checking if mass of profile agrees with the mass of the binary object
        if np.abs(donor_mass[0] - m1_i) > tolerance:
            warnings.warn("Donor mass from the binary class object "
                          "and the profile do not agree")
            print("mass profile/object:", (donor_mass[0]), (m1_i))
        # checking if radius of profile agrees with the radius of the binary
        if np.abs(donor_radius[0] - radius1) > tolerance:
            warnings.warn("Donor radius from the binary class object "
                          "and the profile do not agree")
            print("radius profile/object:", (donor_radius[0]), (radius1))

        # MANOS: if dm exists as a column, else calculate it from mass column
        if "dm" in profile.dtype.names:
            donor_dm = profile["dm"]
            # dm in MESA is in cgs, not in MSsun units, so we transform it
            donor_dm = donor_dm / const.Msun
        else:
            donor_dm = np.concatenate((-1 * np.diff(donor_mass),
                                       [donor_mass[-1]]))
    else:
        raise ValueError("One or many of the mass and/or radius needed columns"
                         " in the profile is not provided for the CEE")
    return donor_mass, donor_radius, donor_dm


def get_internal_energy_from_profile(common_envelope_option_for_lambda,
                                     profile, tolerance=0.001):
    """Calculate the specific internal energy per shell of the donor.

    Parameters
    ----------
    common_envelope_option_for_lambda : str
        Available options:
        * 'default_lambda': using for lambda the constant value of
        common_envelope_lambda_default parameter
        * 'lambda_from_profile_gravitational': calculating the lambda
        parameter from the donor's profile by using the gravitational
        binding energy from the surface to the core
        (needing "mass", and "radius" as columns in the profile)
        * 'lambda_from_profile_gravitational_plus_internal': as above
         but taking into account a factor of common_envelope_alpha_thermal *
         internal energy too in the binding energy (needing also "energy" as
         column in the profile)
        * 'lambda_from_profile_gravitational_plus_internal_minus_recombination'
        as above but not taking into account the recombination energy in the
        internal energy (needing also "y_mass_fraction_He", "x_mass_fraction_H"
        "neutral_fraction_H", "neutral_fraction_He", and "avg_charge_He" as
        column in the profile)
    profile : numpy.array
        Donor's star profile from MESA
    tolerance : float
        The tolerance of numerical difference in two floats when comparing
        and testing results.

    Returns
    -------
    specific_donor_internal_energy : array
        Value of the specific internal energy per shell of the donor

    """
    if (common_envelope_option_for_lambda
            == "lambda_from_profile_gravitational"):
        # initiate specific internal energy as 0
        specific_donor_internal_energy = profile["radius"] * 0.0
    elif ((common_envelope_option_for_lambda
           != "lambda_from_profile_gravitational")
            and (not("energy" in profile.dtype.names))):
        warnings.warn("Profile does not include internal energy -- "
                      "Proceeding with 'lambda_from_profile_gravitational'")
        # initiate specific internal energy as 0
        specific_donor_internal_energy = profile["radius"] * 0.0
    elif ((common_envelope_option_for_lambda
           == "lambda_from_profile_gravitational_plus_internal")
            and ("energy" in profile.dtype.names)):
        # specific internal energy - if we would have used the "total_energy"
        # it would include (internal+potential+kinetic+rotation)
        specific_donor_internal_energy = profile["energy"]
        # if (np.any(specific_donor_internal_energy < -tolerance)):#if negative
        # raise ValueError("CEE problem calculating internal energy, "
        #                  "giving negative values")
        if not (np.any(specific_donor_internal_energy > -tolerance)):
            raise Exception("CEE problem calculating internal energy, "
                            "giving negative values.")

        specific_donor_H2recomb_energy = calculate_H2recombination_energy(
            profile, tolerance)
        # we still need to subtract the H2 recombination energy which is
        # included in the "energy" column of the profile
        # internal_energy - H2 recombination energy per shell
        specific_donor_internal_energy = (
            specific_donor_internal_energy - specific_donor_H2recomb_energy)
        if not (np.any(specific_donor_internal_energy > -tolerance)):
            raise Exception(
                "CEE problem calculating recombination (and H2 recombination) "
                "energy, remaining internal energy giving negative values.")
    elif ((common_envelope_option_for_lambda == "lambda_from_profile_"
           "gravitational_plus_internal_minus_recombination")
            and ("energy" in profile.dtype.names)):
        # specific internal energy. If we have used the "total_energy" it would
        # include (internal+potential+kinetic+rotation)
        specific_donor_internal_energy = profile["energy"]
        if not (np.any(specific_donor_internal_energy > -tolerance)):
            raise Exception("CEE problem calculating internal energy, "
                            "giving negative values.")

        # we still need to subtract the H2 recombination energy which is
        # included in the "energy" column of the profile
        specific_donor_H2recomb_energy = calculate_H2recombination_energy(
            profile, tolerance)
        specific_donor_recomb_energy = calculate_recombination_energy(
            profile, tolerance)
        # internal_energy - recombination energy - H2 recombination energy per
        # shell (so I think it is left with thermal + radiation)
        specific_donor_internal_energy = (
            specific_donor_internal_energy
            - specific_donor_recomb_energy
            - specific_donor_H2recomb_energy)
        # if (np.any(specific_donor_internal_energy < -tolerance)):#if negative
        #     raise ValueError("CEE problem calculating recombination (and H2 "
        #     "recombination) energy, remaining internal energy "
        #     "giving negative values")
        if not (np.any(specific_donor_internal_energy > -tolerance)):
            raise Exception(
                "CEE problem calculating recombination (and H2 recombination) "
                "energy, remaining internal energy giving negative values.")
    return specific_donor_internal_energy


def calculate_H2recombination_energy(profile, tolerance=0.001):
    """Compute the recombination energy of H2 per shell in erg.

    Parameters
    ----------
    profile : array
        Donor's star profile from MESA
    tolerance : float
        The tolerance of numerical difference in two floats when comparing
        and testing results.

    Returns
    -------
    specific_donor_H2recomb_energy : array
        recombination energy of H2 per shell in ergs

    """
    if "x_mass_fraction_H" not in profile.dtype.names:
        warnings.warn("Profile does not include Hydrogen mass fraction "
                      "calculate H2 recombination energy -- "
                      "H2 recombination energy is assumed 0")
        specific_donor_H2recomb_energy = profile["radius"] * 0.0
    else:
        # Dissociation energy [cm^1] from Cheng+2018:
        # https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.121.013001
        specific_donor_H2recomb_energy = (
            35999.582894 * const.inversecm2erg / (2.0 * const.H_weight)
            * profile['x_mass_fraction_H'] * const.avo)
        # http://www.nat.vu.nl/~griessen/STofHinM/ChapIIHatomMoleculeGas.pdf
        # if np.any(specific_donor_H2recomb_energy < -tolerance): # if negative
        #     raise ValueError("CEE problem calculating H2 recombination "
        #     "energy, giving negative values")
        if not (np.any(specific_donor_H2recomb_energy > -tolerance)):
            raise Exception("CEE problem calculating H2 recombination energy, "
                            "giving negative values")
    # return specific_donor_H2recomb_energy * const.ev2erg
    return specific_donor_H2recomb_energy


def calculate_recombination_energy(profile, tolerance=0.001):
    """Compute the recombination energy per shell in erg.

    Parameters
    ----------
    profile : numpy.array
        Donor's star profile from MESA
    tolerance : float
        The tolerance of numerical difference in two floats when comparing
        and testing results.

    Returns
    -------
    specific_donor_recomb_energy : array
        recombination energy per shell in ergs

    """
    if (not(("y_mass_fraction_He" in profile.dtype.names)
            and ("x_mass_fraction_H" in profile.dtype.names)
            and ("neutral_fraction_H" in profile.dtype.names)
            and ("neutral_fraction_He" in profile.dtype.names)
            and ("avg_charge_He" in profile.dtype.names))):
        warnings.warn("Profile does not include mass fractions and ionizations"
                      " of elements to calculate recombination energy "
                      "-- recombination energy is assumed 0")
        specific_donor_recomb_energy = profile["radius"] * 0.0
    else:
        # from MESA/binary/private/binary_ce.f90
        frac_HI = profile["neutral_fraction_H"]
        for i in range(len(frac_HI)):
            # TODO: For some reason this is a bit > 1 in the outer envelope.
            #       Maybe not important but better to solve it.
            frac_HI[i] = min(1., frac_HI[i])
        frac_HII = 1.0 - frac_HI

        frac_HeI = profile["neutral_fraction_He"]
        avg_charge_He = profile["avg_charge_He"]
        for i in range(len(frac_HI)):
            frac_HeI[i] = min(1., frac_HeI[i])
            # knowing the frac_HeI and the avg_charge_He,
            # we can solve for frac_HeII and frac_HeIII
            avg_charge_He[i] = max(avg_charge_He[i], 0.0)
        frac_HeII = 2.0 - 2.0 * frac_HeI - avg_charge_He
        frac_HeIII = 1.0 - frac_HeII - frac_HeI

        specific_donor_recomb_energy = profile_recomb_energy(
            profile['x_mass_fraction_H'], profile['y_mass_fraction_He'],
            frac_HII, frac_HeII, frac_HeIII)
        # if (np.any(specific_donor_recomb_energy < -tolerance)): # if negative
        #     print(specific_donor_recomb_energy)
        #     raise ValueError("CEE problem calculating recombination energy,"
        #     " giving negative values")
        if not (np.any(specific_donor_recomb_energy > -tolerance)):
            raise Exception("CEE problem calculating recombination energy, "
                            "giving negative values.")
    return specific_donor_recomb_energy


def profile_recomb_energy(x_mass_fraction_H, y_mass_fraction_He,
                          frac_HII, frac_HeII, frac_HeIII):
    """Calculate the recombination energy per shell.

    Parameters
    ----------
    x_mass_fraction_H : array
        Mass fraction of H per shell
    y_mass_fraction_He : array
        Mass fraction of He per shell.
    frac_HII : array
        Fraction of single ionized H per shell.
    frac_HeII : array
        Fraction of single ionized He per shell.
    frac_HeIII : array
        Fraction of double ionized He per shell.

    Returns
    -------
    recomb_energy : 1D numpy.arrays
        recombination energy per shell in ergs
        (same dimension as x_mass_fraction_H and y_mass_fraction_He)
    """
    recomb_energy = (
        54.4177650 * frac_HeIII / const.He_weight
        * y_mass_fraction_He * const.avo
        + 24.587388 * (frac_HeII + frac_HeIII) / const.He_weight
        * y_mass_fraction_He * const.avo
        + 13.598434 * frac_HII / const.H_weight
        * x_mass_fraction_H * const.avo)
    return recomb_energy * const.ev2erg


def calculate_binding_energy(donor_mass, donor_radius, donor_dm,
                             specific_internal_energy,
                             ind_core,
                             factor_internal_energy,
                             verbose, tolerance=0.001
                             ):
    """Calculate the total binding energy of the envelope of the star.

    Parameters
    ----------
    donor_mass : array
        Enclosed mass coordinate of the donor star from its profile.
    donor_radius : array
        Radius coordinate of the donor star from its profile.
    donor_dm : array
        Mass enclosed per shell in the donor star's profile.
    specific_internal_energy : array
        Specific internal energy of the donor star.
    ind_core : int
        The value of the shell position of the core - envelope boundary,
        at the donor's MESA profile
    factor_internal_energy : float
        The factor to multiply with internal energy to be taken into
        account when we calculate the  binding energy of the enevelope
    verbose : bool
        In case we want information about the CEE  (the default is False).

    Returns
    -------
    Ebind_i : float
        The total binding energy of the envelope of the star

    """
    # Sum of gravitational energy from surface to core boundary
    Grav_energy = 0.0
    # Sum of internal energy from surface to core boundary. This is 0 if
    # 'lambda_from_profile_gravitational' or (thermal+radiation+recombination)
    # for "lambda_from_profile_gravitational_plus_internal" or
    # (thermal+radiation) for
    # "lambda_from_profile_gravitational_plus_internal_minus_recombination"
    U_i = 0.0
    # sum from surface to the core. Your core boundary is in element [ind_core]
    # in a normal MESA (and POSYDON) profile
    for i in range(ind_core):
        Grav_energy_of_cell = (-const.standard_cgrav * donor_mass[i]
                               * const.Msun * donor_dm[i]*const.Msun
                               / (donor_radius[i]*const.Rsun))
        # integral of gravitational energy as we go deeper into the star
        Grav_energy = Grav_energy + Grav_energy_of_cell
        U_i = U_i + specific_internal_energy[i]*donor_dm[i]*const.Msun
    if Grav_energy > 0.0:
        print("Grav_energy, donor_mass, donor_dm, donor_radius",
              Grav_energy, donor_mass, donor_dm, donor_radius)
        if not (Grav_energy < tolerance):
            raise Exception("CEE problem calculating gravitational energy, "
                            "giving positive values.")
    # binding energy of the enevelope equals its gravitational energy +
    # an a_th fraction of its internal energy
    Ebind_i = Grav_energy + factor_internal_energy * U_i
    if Ebind_i > -tolerance:
        warnings.warn("Ebind_i of the envelope is found positive")
    if verbose:
        print("integration of gravitational energy surface to core "
              "[Grav_energy], integration of internal energy surface to "
              "core [U_i] (0 if not taken into account) ", Grav_energy, U_i)
        print("Ebind = Grav_energy + factor_internal_energy*U_i  :  ", Ebind_i)
    return Ebind_i

def calculate_Mejected_for_integrated_binding_energy(profile, Ebind_threshold,
                             mc1_i, rc1_i,
                             m1_i = 0.0, radius1 = 0.0,
                             factor_internal_energy=1.0,tolerance=0.001
                             ):
    """Calculate the mass lost from the envelope for an energy budget of Ebind_threshold

    Parameters
    ----------
    profile : numpy.array
        Donor's star profile from MESA
    Ebind_threshold : float
        Orbital energy used from the spiral in to partial unbind the envelope. Positive
        We integrate from surface to calcualte the partial loss of mass during CE that merges.
    factor_internal_energy : float
        The factor to multiply with internal energy to be taken into
        account when we calculate the  binding energy of the enevelope
    verbose : bool
        In case we want information about the CEE  (the default is False).

    Returns
    -------
    Ebind_i : float
        The total binding energy of the envelope of the star

    """

    donor_mass, donor_radius, donor_dm = get_mass_radius_dm_from_profile(
        profile, m1_i, radius1, tolerance)
    specific_internal_energy = get_internal_energy_from_profile(
        common_envelope_option_for_lambda = "lambda_from_profile_gravitational_plus_internal_minus_recombination",
        profile = profile, tolerance = tolerance)

    # Sum of gravitational energy from surface towards inside
    Grav_energy = 0.0
    # Sum of internal energy from surface towards.
    U_energy = 0.0
    # sum from surface to the core. Your threshold is in element [ind_threshold]
    # in a normal MESA (and POSYDON) profile
    i = 0
    Ebind_so_far = 0.0 # the integration from surface going inwards of the binding energy (negative in principle)

    while (abs(Ebind_so_far) < Ebind_threshold) and (i<len(donor_mass)):
        Grav_energy_of_cell = (-const.standard_cgrav * donor_mass[i]
                               * const.Msun * donor_dm[i]*const.Msun
                               / (donor_radius[i]*const.Rsun))
        # integral of gravitational energy as we go deeper into the star
        Grav_energy = Grav_energy + Grav_energy_of_cell
        U_energy = U_energy + specific_internal_energy[i]*donor_dm[i]*const.Msun
        Ebind_so_far = Grav_energy + factor_internal_energy * U_energy
        i=i+1
    ind_threshold = i-1

    if donor_mass[ind_threshold]< mc1_i or  donor_radius[ind_threshold]<rc1_i:
        warnings.warn("partial mass ejected found more than the envelope mass")
        print("M_ejected, M_envelope = ", donor_mass[0] - donor_mass[ind_threshold], donor_mass[0] - mc1_i)
        donor_mass[ind_threshold] = mc1_i

    M_ejected = donor_mass[0] - donor_mass[ind_threshold]

    return M_ejected


class PchipInterpolator2:
    """Interpolation class."""

    def __init__(self, *args, positive=False, **kwargs):
        """Initialize the interpolator."""
        self.interpolator = PchipInterpolator(*args, **kwargs)
        self.positive = positive

    def __call__(self, *args, **kwargs):
        """Use the interpolator."""
        result = self.interpolator(*args, **kwargs)
        if self.positive:
            result = np.maximum(result, 0.0)
        return result

def convert_metallicity_to_string(Z):
    """Check if metallicity is supported by POSYDON v2."""
    # check supported metallicity
    valid_Z = [2e+00,1e+00,4.5e-01,2e-01,1e-01,1e-02,1e-03,1e-04]
    if not Z in valid_Z:
        raise ValueError(f'Metallicity {Z} not supported! Available metallicities in POSYDON v2 are {valid_Z}.')
    return f'{Z:1.1e}'.replace('.0','')
