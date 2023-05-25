"""Magnetic Braking."""


__authors__ = [
    "Seth Gossage <seth.gossage@northwestern.edu>"
    "Jeffrey Andrews <jeffrey.andrews@ufl.edu>",
]


import numpy as np

import posydon.utils.constants as const


def calculate_magnetic_braking(p_pri, p_sec, magnetic_braking_mode):
    """Calculate the impact of magnetic braking on a binary.

    domega_mb / dt = torque_mb / I is calculated, i.e., the amount of change
    in Omega over 1 year.

    Parameters
    ----------
    p_pri, p_sec : two tuples
        The properties of the primary and secondary stars, respectively. Each
        tuple contains, in order, M_pri (the star's mass in Msun), R_pri
        (the star's radius in Rsun), Omega_pri (the star's angular velocity in
        1/yr), I_pri (the star's moment of inertia in Msun Rsun^2),
        tau_conv_pri (the convective timescale of the primary),
        and Mdot_pri (the mass loss rate of the star in Msun/yr)
    magnetic_braking_mode : string
        The magnetic braking prescription used.

    Returns
    -------
    dOmega_mb_sec, dOmega_mb_pri : float
        The spin angular velocity derivative for the secondary and primary
        stars, respectively (units of yr^-2)

    """
    if magnetic_braking_mode == "RVJ83":
        return calculate_magnetic_braking_RVJ83(p_pri, p_sec)
    elif magnetic_braking_mode == "M15":
        return calculate_magnetic_braking_M15(p_pri, p_sec)
    elif magnetic_braking_mode == "G18":
        return calculate_magnetic_braking_G18(p_pri, p_sec)
    elif magnetic_braking_mode == "CARB":
        return calculate_magnetic_braking_CARB(p_pri, p_sec)
    else:
        raise Exception("WARNING: Magnetic braking is not being calculated in"
                        "the detached step. The given magnetic_braking_mode"
                        "string \"", magnetic_braking_mode, "\" does not match"
                        "the available built-in cases. To enable magnetic"
                        "braking, please set magnetic_braking_mode to one of"
                        "the following strings:\n"
                        "\"RVJ83\" for Rappaport, Verbunt, & Joss 1983\n"
                        "\"G18\" for Garraffo et al. 2018\n"
                        "\"M15\" for Matt et al. 2015\n"
                        "\"CARB\" for Van & Ivanova 2019\n")


def calculate_magnetic_braking_RVJ83(p_pri, p_sec):
    """Use the Rappaport, Verbunt, and Joss (1983) prescription.

    Calculate the impact of magnetic braking on a binary using the
    Rappaport, Verbunt, and Joss 1983, ApJ, 275, 713 prescription. The torque
    is taken from eq.36 of Rapport+1983, with γ = 4. Torque units
    # converted from cgs units to [Msol], [Rsol], [yr] as all stellar
    # parameters are given in units of [Msol], [Rsol], [yr] and so that
    # dOmega_mb/dt is in units of [yr^-2].

    Parameters
    ----------
    p_pri, p_sec : two tuples
        The properties of the primary and secondary stars, respectively. Each
        tuple contains, in order, M_pri (the star's mass in Msun), R_pri
        (the star's radius in Rsun), Omega_pri (the star's angular velocity in
        1/yr), I_pri (the star's moment of inertia in Msun Rsun^2),
        tau_conv_pri (the convective timescale of the primary),
        and Mdot_pri (the mass loss rate of the star in Msun/yr)

    Returns
    -------
    dOmega_mb_sec, dOmega_mb_pri : float
        The spin angular velocity derivative for the secondary and primary
        stars, respectively [units of yr^-2]

    """
    M_pri, R_pri, Omega_pri, I_pri, tau_conv_pri, Mdot_pri = p_pri
    M_sec, R_sec, Omega_sec, I_sec, tau_conv_sec, Mdot_sec = p_sec

    # Converting units:
    # The constant 3.8e-30 from Rappaport+1983 has units of [cm^-2 s]
    # which need to be converted...
    #
    # -3.8e-30 [cm^-2 s] * (const.rsol**2/const.secyer) -> [Rsol^-2 yr]
    # * M [Msol]
    # * R ** 4 [Rsol^4]
    # * Omega ** 3 [yr^-3]
    # / I [Msol Rsol^2 ]
    #
    # Thus, dOmega/dt comes out to [yr^-2]

    dOmega_mb_sec = (
        -3.8e30 * (const.rsol**2 / const.secyer)
        * M_sec
        * R_sec**4
        * Omega_sec**3
        / I_sec
        * np.clip((1.5 - M_sec) / (1.5 - 1.3), 0, 1)
    )
    dOmega_mb_pri = (
        -3.8e30 * (const.rsol**2 / const.secyer)
        * M_pri
        * R_pri**4
        * Omega_pri**3
        / I_pri
        * np.clip((1.5 - M_pri) / (1.5 - 1.3), 0, 1)
    )

    return dOmega_mb_sec, dOmega_mb_pri


def calculate_magnetic_braking_M15(p_pri, p_sec):
    """Use the Matt et al. (2015) prescription.

    Calculate the impact of magnetic braking on a binary using the
    Matt et al. 2015, ApJ, 799, 23 prescription.

    Parameters
    ----------
    p_pri, p_sec : two tuples
        The properties of the primary and secondary stars, respectively. Each
        tuple contains, in order, M_pri (the star's mass in Msun), R_pri
        (the star's radius in Rsun), Omega_pri (the star's angular velocity in
        1/yr), I_pri (the star's moment of inertia in Msun Rsun^2),
        tau_conv_pri (the convective timescale of the primary),
        and Mdot_pri (the mass loss rate of the star in Msun/yr)

    Returns
    -------
    dOmega_mb_sec, dOmega_mb_pri : float
        The spin angular velocity derivative for the secondary and primary
        stars, respectively [units of yr^-2]

    """
    M_pri, R_pri, Omega_pri, I_pri, tau_conv_pri, Mdot_pri = p_pri
    M_sec, R_sec, Omega_sec, I_sec, tau_conv_sec, Mdot_sec = p_sec

    # Torque prescription from Matt et al. 2015, ApJ, 799, L23
    # Constants:
    # [erg] or [g cm^2 s^-2] -> [Msol Rsol^2 yr^-2]
    K = -1.4e30 * const.secyer**2 / (const.msol * const.rsol**2)
    # m = 0.22
    # p = 2.6
    # Above constants were calibrated as in
    # Gossage et al. 2021, ApJ, 912, 65

    # TODO: I am not sure which constants are used from each reference

    # Below, constants are otherwise as assumed as in
    # Matt et al. 2015, ApJ, 799, L23
    omega_sol = 2.6e-6 * const.secyer   # [s^-1] -> [yr^-1]
    # solar rossby = 2
    # solar convective turnover time = 12.9 days
    # Rossby number saturation threshold = 0.14
    chi = 2.0 / 0.14
    tau_conv_sol = 12.9 / 365.25        # 12.9 [days] -> [yr]

    Prot_pri = 2 * np.pi / Omega_pri    # [yr]
    Rossby_number_pri = Prot_pri / tau_conv_pri
    Prot_sec = 2 * np.pi / Omega_sec    # [yr]
    Rossby_number_sec = Prot_sec / tau_conv_sec

    # critical rotation rate in rad/yr
    Omega_crit_pri = np.sqrt(
        const.standard_cgrav * M_pri * const.msol
        / ((R_pri * const.rsol) ** 3)) * const.secyer
    Omega_crit_sec = np.sqrt(
        const.standard_cgrav * M_sec * const.msol
        / ((R_sec * const.rsol) ** 3)) * const.secyer

    # omega/omega_c
    wdivwc_pri = Omega_pri / Omega_crit_pri
    wdivwc_sec = Omega_sec / Omega_crit_sec

    gamma_pri = (1 + (wdivwc_pri / 0.072)**2)**0.5
    T0_pri = K * R_pri**3.1 * M_pri**0.5 * gamma_pri**(-2 * 0.22)
    gamma_sec = (1 + (wdivwc_sec / 0.072)**2)**0.5
    T0_sec = K * R_sec**3.1 * M_sec**0.5 * gamma_sec**(-2 * 0.22)

    if (Rossby_number_sec < 0.14):
        dOmega_mb_sec = (
            T0_sec * (chi**2.6) * (Omega_sec / omega_sol) / I_sec
            * np.clip((1.5 - M_sec) / (1.5 - 1.3), 0, 1)
        )
    else:
        dOmega_mb_sec = (
            T0_sec * ((tau_conv_sec/tau_conv_sol)**2.6)
                   * ((Omega_sec/omega_sol)**(2.6 + 1)) / I_sec
            * np.clip((1.5 - M_sec) / (1.5 - 1.3), 0, 1)
        )

    if (Rossby_number_pri < 0.14):
        dOmega_mb_pri = (
            T0_pri * (chi**2.6) * (Omega_pri / 2.6e-6) / I_pri
            * np.clip((1.5 - M_pri) / (1.5 - 1.3), 0, 1)
        )
    else:
        dOmega_mb_pri = (
            T0_pri * ((tau_conv_pri/tau_conv_sol)**2.6)
                   * ((Omega_pri/omega_sol)**(2.6 + 1)) / I_pri
            * np.clip((1.5 - M_pri) / (1.5 - 1.3), 0, 1)
        )

    return dOmega_mb_sec, dOmega_mb_pri


def calculate_magnetic_braking_G18(p_pri, p_sec):
    """Use the Garraffo et al. (2018) prescription.

    Calculate the impact of magnetic braking on a binary using the
    Garraffo et al. 2018, ApJ, 862, 90 prescription. We adopt the following
    constants:
    a = 0.03
    b = 0.5

    Parameters
    ----------
    p_pri, p_sec : two tuples
        The properties of the primary and secondary stars, respectively. Each
        tuple contains, in order, M_pri (the star's mass in Msun), R_pri
        (the star's radius in Rsun), Omega_pri (the star's angular velocity in
        1/yr), I_pri (the star's moment of inertia in Msun Rsun^2),
        tau_conv_pri (the convective timescale of the primary),
        and Mdot_pri (the mass loss rate of the star in Msun/yr)

    Returns
    -------
    dOmega_mb_sec, dOmega_mb_pri : float
        The spin angular velocity derivative for the secondary and primary
        stars, respectively [units of yr^-2]

    """
    M_pri, R_pri, Omega_pri, I_pri, tau_conv_pri, Mdot_pri = p_pri
    M_sec, R_sec, Omega_sec, I_sec, tau_conv_sec, Mdot_sec = p_sec

    a_const = 0.03
    b_const = 0.5

    # [g cm^2] -> [Msol Rsol^2]
    c = -3e41 / (const.msol * const.rsol**2)
    # Above are as calibrated in Gossage et al. 2021, ApJ, 912, 65

    Prot_pri = 2 * np.pi / Omega_pri            # [yr]
    Rossby_number_pri = Prot_pri / tau_conv_pri
    Prot_sec = 2 * np.pi / Omega_sec            # [yr]
    Rossby_number_sec = Prot_sec / tau_conv_sec

    n_pri = (a_const / Rossby_number_pri) + b_const * Rossby_number_pri + 1.0
    n_sec = (a_const / Rossby_number_sec) + b_const * Rossby_number_sec + 1.0

    Qn_pri = 4.05 * np.exp(-1.4 * n_pri)
    Qn_sec = 4.05 * np.exp(-1.4 * n_sec)

    dOmega_mb_sec = (
            c * Omega_sec**3 * tau_conv_sec * Qn_sec / I_sec
            * np.clip((1.5 - M_sec) / (1.5 - 1.3), 0, 1)
    )

    dOmega_mb_pri = (
            c * Omega_pri**3 * tau_conv_pri * Qn_pri / I_pri
            * np.clip((1.5 - M_sec) / (1.5 - 1.3), 0, 1)
    )

    return dOmega_mb_sec, dOmega_mb_pri


def calculate_magnetic_braking_CARB(p_pri, p_sec):
    """Use the CARB prescription.

    Calculate the impact of magnetic braking on a binary using the
    Van & Ivanova 2019, ApJ, 886, L31 prescription. Our prescription is based
    on files hosted on Zenodo: https://zenodo.org/record/3647683#.Y_TfedLMKUk,
    after converting from cgs to solar units.

    Parameters
    ----------
    p_pri, p_sec : two tuples
        The properties of the primary and secondary stars, respectively. Each
        tuple contains, in order, M_pri (the star's mass in Msun), R_pri
        (the star's radius in Rsun), Omega_pri (the star's angular velocity in
        1/yr), I_pri (the star's moment of inertia in Msun Rsun^2),
        tau_conv_pri (the convective timescale of the primary),
        and Mdot_pri (the mass loss rate of the star in Msun/yr)

    Returns
    -------
    dOmega_mb_sec, dOmega_mb_pri : float
        The spin angular velocity derivative for the secondary and primary
        stars, respectively [units of yr^-2]

    """
    M_pri, R_pri, Omega_pri, I_pri, tau_conv_pri, Mdot_pri = p_pri
    M_sec, R_sec, Omega_sec, I_sec, tau_conv_sec, Mdot_sec = p_sec

    # Constants as assumed in Van & Ivanova 2019, ApJ, 886, L31
    # with units converted from [cm], [g], [s] to [Rsol], [Msol], [yr]
    omega_sol = 3e-6 * const.secyer         # [s^-1] -> [yr^-1]
    tau_conv_sol = 2.8e6 / const.secyer     # [s] -> yr
    K2 = 0.07**2

    tau_ratio_sec = tau_conv_sec / tau_conv_sol
    tau_ratio_pri = tau_conv_pri / tau_conv_sol
    rot_ratio_sec = Omega_sec / omega_sol
    rot_ratio_pri = Omega_pri / omega_sol

    # below in units of [Rsol yr^-1]^2
    v_esc2_sec = ((2 * const.standard_cgrav * M_sec / R_sec)
                  * (const.msol * const.secyer**2 / const.rsol**3))
    v_esc2_pri = ((2 * const.standard_cgrav * M_pri / R_pri)
                  * (const.msol * const.secyer**2 / const.rsol**3))
    v_mod2_sec = v_esc2_sec + (2 * Omega_sec**2 * R_sec**2) / K2
    v_mod2_pri = v_esc2_pri + (2 * Omega_pri**2 * R_pri**2) / K2

    # Van & Ivanova 2019, MNRAS 483, 5595 replace the magnetic field
    # with Omega * tau_conv phenomenology. Thus, the ratios
    # (rot_ratio_* and tau_ratio_*) inherently have units of Gauss
    # [cm^-0.5 g^0.5 s^-1] that needs to be converted to [Rsol],
    # [Msol], [yr]. VI2019 assume the solar magnetic field strength is
    # on average 1 Gauss.
    if (abs(Mdot_sec) > 0):
        R_alfven_div_R3_sec = (
            R_sec**4 * rot_ratio_sec**4 * tau_ratio_sec**4
            / (Mdot_sec**2 * v_mod2_sec)
            * (const.rsol**2 * const.secyer / const.msol**2))
    else:
        R_alfven_div_R3_sec = 0.0

    if (abs(Mdot_pri) > 0):
        R_alfven_div_R3_pri = (
            R_pri**4 * rot_ratio_pri**4 * tau_ratio_pri**4
            / (Mdot_pri**2 * v_mod2_pri)
            * (const.rsol**2 * const.secyer / const.msol**2))
    else:
        R_alfven_div_R3_pri = 0.0

    # Alfven radius in [Rsol]
    R_alfven_sec = R_sec * R_alfven_div_R3_sec**(1./3.)
    R_alfven_pri = R_pri * R_alfven_div_R3_pri**(1./3.)

    dOmega_mb_sec = (
          (2./3.) * Omega_sec * Mdot_sec * R_alfven_sec**2 / I_sec
          * np.clip((1.5 - M_sec) / (1.5 - 1.3), 0, 1)
    )

    dOmega_mb_pri = (
          (2./3.) * Omega_pri * Mdot_pri * R_alfven_pri**2 / I_pri
          * np.clip((1.5 - M_sec) / (1.5 - 1.3), 0, 1)
    )

    return dOmega_mb_sec, dOmega_mb_pri
