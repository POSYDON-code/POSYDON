"""Tides."""

__authors__ = [
    "Zepei Xing <Zepei.Xing@unige.ch>",
    "Jeffrey Andrews <jeffrey.andrews@northwestern.edu>",
]


import numpy as np

from posydon.utils.common_functions import orbital_period_from_separation
import posydon.utils.constants as const


def calculate_kT_convective(P_orb, x_star, verbose=False):
    """Calculate the tidal k/T term for convective stars.

    Timescale is calculated from eq. (31) of Hurley et al. 2002, generalized
    for convective layers. k/T is calculated from eq. (30) of Hurley et al.
    2002

    Parameters
    ----------
    P_orb : float
        Orbital period
    x_star : tuple
        Tuple containing quantities related to the present star.
    verbose : bool
        Flag indicating whether verbose outputs are turned on

    Returns
    -------
    kT_conv : float
        k/T term for convective stars

    """
    # Unpack variables
    M, R, M_env, DR_env, Renv_middle, L, Omega, I, conv_mx1_top_r, \
        conv_mx1_bot_r, surface_h1 = x_star

    # Check envelope mass
    if M_env <= 0.0 or np.isnan(M_env) or np.isinf(M_env):
        if verbose:
            print("Problematic envelope mass (M_env):", M_env)
        return 1e99

    # Check DR_env
    if DR_env <= 0.0 or np.isnan(DR_env) or np.isinf(DR_env):
        if verbose:
            print("Problematic DR envelope (DR_env):", DR_env)
        return 1e99

    # Check Renv_middle
    if Renv_middle <= 0.0 or np.isnan(Renv_middle) or np.isinf(Renv_middle):
        if verbose:
            print("Problematic middle radius of the envelope (Renv_middle):",
                  Renv_middle)
        return 1e99

    # If all checks are complete, calculate convective timescale
    tau_conv = 0.431 * ((M_env * DR_env * Renv_middle / (3 * L)) ** (1/3))

    # Calculate periods
    P_spin = 2 * np.pi / Omega
    P_tid = np.abs(1 / (1 / P_orb - 1 / P_spin))
    f_conv = np.min([1, (P_tid / (2 * tau_conv)) ** 2])

    # Calculate the convective k/T term
    F_tid = 1.  # not 50 as before
    kT_conv = (2. / 21) * (f_conv / tau_conv) * (M_env / M)

    # Run a check on the calculated kT_conv
    if kT_conv is None or not np.isfinite(kT_conv):
        kT_conv = 0.0

    # this is the 1/timescale of all d/dt calculted below in yr^-1
    if verbose and verbose != 1:
        print("Equilibrium tides in deep convective envelope",
              M_env, DR_env, Renv_middle, M)
        print("convective tiimescales and efficiencies",
              tau_conv, P_orb, P_spin, P_tid, f_conv, F_tid, kT_conv)

    return kT_conv


def calculate_kT_radiative(a, P_orb, q, x_star, verbose=False):
    """Calculate the tidal k/T term for radiative stars.

    Timescale is calculated from eq. (31) of Hurley et al. 2002, generalized
    for convective layers. k/T is calculated from eq. (30) of Hurley et al.
    2002

    Parameters
    ----------
    a : float
        Orbital separation
    P_orb : float
        Orbital period
    q : float
        Mass ratio
    x_star : tuple
        Tuple containing quantities related to the present star.
    verbose : bool
        Flag indicating whether verbose outputs are turned on

    Returns
    -------
    kT_rad : float
        k/T term for radiative stars

    """

    def calculate_E2(a, q, M, R, conv_mx1_top_r, conv_mx1_bot_r, surface_h1,
                     verbose=False):
        """Calculate the tidal constant E2.

        Timescale is calculated from eq. 9 in Qin et al. 2018, 616, A28

        Parameters
        ----------
        a : float
            orbital separation
        q : float
            Mass ratio
        M : float
            Stellar mass
        R : float
            Stellar radius
        conv_mx1_bot_r : float
            radius of the bottom of the outermost convective zone
        conv_mx1_top_r : float
            radius of the top of the outermost convective zone
        surface_h1 : float
            The h1 fraction of the surface
        verbose : bool
            Verbose outputs flag (default: false)

        Returns
        -------
        E2 : float
            E2 tidal constant, used for radiative stars

        """
        R_conv = conv_mx1_top_r - conv_mx1_bot_r

        if (R_conv > R or R_conv <= 0.0 or conv_mx1_bot_r / R > 0.1):
            # we switch to Zahn+1975 calculation of E2
            E2 = 1.592e-9 * M ** (2.84)
            if verbose and verbose != 1:
                print("R_conv of the convective core is not behaving well or "
                      "we are not calculating the convective core, we switch "
                      "to Zahn+1975 calculation of E2",
                      R_conv, R, conv_mx1_top_r, conv_mx1_bot_r, E2)
        else:
            if R <= 0:
                E2 = 0
            elif surface_h1 > 0.4:
                # eq. (9) of Qin et al. 2018, 616, A28
                E2 = 10.0 ** (-0.42) * (R_conv / R) ** 7.5
            elif surface_h1 <= 0.4:  # "HeStar":
                # eq. (9) of Qin et al. 2018, 616, A28
                E2 = 10.0 ** (-0.93) * (R_conv / R) ** 6.7
            else:  # in principle we should not go here
                # eq. (43) of Hurley et al. 2002 from Zahn+1975 Depreciated
                E2 = 1.592e-9 * M ** 2.84
        # kT = 1.9782e4 * np.sqrt(M * R**2 / a**5) * (1 + q)**(5. / 6) * E2
        # eq. (42) of Hurley et al. 2002. Depreciated

        return E2

    # Unpack variables
    M, R, M_env, DR_env, Renv_middle, L, Omega, I, conv_mx1_top_r, \
        conv_mx1_bot_r, surface_h1 = x_star

    # First, calculate E2
    E2 = calculate_E2(a, q, M, R, conv_mx1_top_r, conv_mx1_bot_r, surface_h1)

    # Now, calculate the k/T term
    kT_rad = np.sqrt(const.standard_cgrav * (M * const.msol) *
        (R * const.rsol)**2 / (a * const.rsol)**5.0 ) * (1 + q) ** (5.0 / 6)
    kT_rad = kT_rad * E2 * const.secyer

    return kT_rad


def calculate_tidal_derivatives(x_orb, x_star, kT, F_tide, q, n,
                                f1, f2, f3, f4, f5):
    """Calculate the tidal derivatives.

    Tidal effects are calculated using the equilibrium tidal prescription
    from Hut (1981).

    Parameters
    ----------
    x_orb : tuple
        Tuple containing quantities related to orbit.
    x_star : tuple
        Tuple containing quantities related to the present star.
    kT : float
        k/T term of the star
    F_tide : float
        Strength of tides
    q : float
        Mass ratio
    n : float
        Mean motion resonance
    f1, f2, f3, f4, f5 : float
        Terms in eccentricity expansion of derivatives

    Returns
    -------
    da_tides : float
        The derivative of the orbital separation (in Rsun/yr) due to tidal
        interaction with a companion.
    de_tides : float
        The derivative of the eccentricity due to tidal interaction with a
        companion.
    dOmega_tides : float
        The rotational angular velocity of the present star (in 1/yr).

    """
    # Unpack variables
    a, e = x_orb
    M, R, M_env, DR_env, Renv_middle, L, Omega, I, conv_mx1_top_r, \
        conv_mx1_bot_r, surface_h1 = x_star

    # eq. (9) of Hut 1981, 99, 126
    da_tides = -6 * F_tide * kT * q * (1 + q) * (R / a) ** 8 \
        * (a / (1 - e ** 2) ** (15 / 2)) \
        * (f1 - (1 - e ** 2) ** (3 / 2) * f2 * Omega / n)

    # eq. (10) of Hut 1981, 99, 126
    de_tides = -27 * F_tide * kT * q * (1 + q) * (R / a) ** 8\
        * (e / (1 - e ** 2) ** (13 / 2))\
        * (f3 - (11 / 18) * (1 - e ** 2) ** (3 / 2) * f4 * Omega/n)

    # eq. (11) of Hut 1981, 99, 126
    dOmega_tides = 3 * F_tide * kT * q ** 2 * M * R ** 2 / I\
        * (R / a) ** 6 * n / (1 - e ** 2) ** 6\
        * (f2 - (1 - e ** 2) ** (3 / 2) * f5 * Omega / n)

    return da_tides, de_tides, dOmega_tides


def calculate_tides(x_orb, x_pri, x_sec, F_tide=1, verbose=False):
    """Calculate the orbital effect of tides.

    Tidal effects are calculated using the equilibrium tidal prescription
    from Hut (1981).

    Parameters
    ----------
    x_orb : tuple
        Tuple containing quantities related to orbit.
    x_pri : tuple
        Tuple containing quantities related to primary star.
    x_sec : tuple
        Tuple containing quantities related to secondary star.
    F_tide : float
        Strength factor of the tidal force (default: 1)
    verbose : bool
        Turn on verbosity (default: False)

    Returns
    -------
    da_tides_sec : float
        The derivative of the orbital separation (in Rsun/yr) due to tidal
        interaction with the secondary.
    de_tides_sec : float
        The derivative of the eccentricity due to tidal interaction with the
        secondary.
    dOmega_tides_sec : float
        The rotational angular velocity of the primary star (in 1/yr).
    da_tides_pri : float
        The derivative of the orbital separation (in Rsun/yr) due to tidal
        interaction with the primary.
    de_tides_pri : float
        The derivative of the eccentricity due to tidal interaction with the
        primary.
    dOmega_tides_pri : float
        The rotational angular velocity of the primary star (in 1/yr).

    """
    # Unpack variables
    a, e = x_orb
    M_pri, R_pri, M_env_pri, DR_env_pri, Renv_middle_pri, L_pri, Omega_pri, \
        I_pri, conv_mx1_top_r_pri, conv_mx1_bot_r_pri, surface_h1_pri \
        = x_pri
    M_sec, R_sec, M_env_sec, DR_env_sec, Renv_middle_sec, L_sec, Omega_sec, \
        I_sec, conv_mx1_top_r_sec, conv_mx1_bot_r_sec, surface_h1_sec \
        = x_sec

    q1 = M_pri / M_sec
    q2 = M_sec / M_pri
    P_orb = orbital_period_from_separation(a, M_pri, M_sec)
    n = 2.0 * const.pi / P_orb  # mean orbital ang. vel. in rad/year

    # Eccentricity expansion terms
    f1 = 1 + (31 / 2) * e ** 2 + (255 / 8) * e ** 4 + (185 / 16) * e ** 6 \
        + (25 / 64) * e ** 8
    f2 = 1 + (15 / 2) * e ** 2 + (45 / 8) * e ** 4 + (5 / 16) * e ** 6
    f3 = 1 + (15 / 4) * e ** 2 + (15 / 8) * e ** 4 + (5 / 64) * e ** 6
    f4 = 1 + (3 / 2) * e ** 2 + (1 / 8) * e ** 4
    f5 = 1 + 3 * e ** 2 + (3 / 8) * e ** 4

    # Calculate the convective k/T terms
    kT_conv_sec = calculate_kT_convective(P_orb, x_sec, verbose)
    kT_conv_pri = calculate_kT_convective(P_orb, x_pri, verbose)

    # Calculate the radiative k/T terms
    kT_rad_sec = calculate_kT_radiative(a, P_orb, q1, x_sec, verbose)
    kT_rad_pri = calculate_kT_radiative(a, P_orb, q2, x_pri, verbose)

    # Use the maximum of the two k/T terms
    kT_sec = max(kT_conv_sec, kT_rad_sec)
    kT_pri = max(kT_conv_pri, kT_rad_pri)

    if verbose and verbose != 1:
        print("kT_conv/rad of tides is ", kT_conv_sec, kT_rad_sec,
              kT_conv_pri, kT_rad_pri, "in 1/yr, and we picked the ",
              kT_sec, kT_pri)

    da_tides_sec, de_tides_sec, dOmega_tides_sec = \
        calculate_tidal_derivatives(x_orb, x_sec, kT_sec, F_tide, q1, n,
                                    f1, f2, f3, f4, f5)

    da_tides_pri, de_tides_pri, dOmega_tides_pri = \
        calculate_tidal_derivatives(x_orb, x_pri, kT_pri, F_tide, q2, n,
                                    f1, f2, f3, f4, f5)

    return da_tides_sec, de_tides_sec, dOmega_tides_sec, da_tides_pri, \
        de_tides_pri, dOmega_tides_pri
