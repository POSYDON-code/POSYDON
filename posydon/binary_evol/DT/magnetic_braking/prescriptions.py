import numpy as np

import posydon.utils.constants as const


def RVJ83_braking(primary, secondary, verbose = False):
    """
        Calculates the change in spin of each star due
    to magnetic braking, according to:

        Rappaport, S., Joss, P. C., & Verbunt, F. 1983, ApJ, 275, 713

    Parameters
    ----------
        primary : SingleStar object
            A single star object, representing the primary (more evolved) star
            in the binary and containing its properties.

        secondary : SingleStar object
            A single star object, representing the secondary (less evolved) star
            in the binary and containing its properties.

        verbose : bool
            True if we want to print stuff.

    Returns
    -------
        dOmega_sec : float
            The change in rotation rate of the secondary (less evolved)
        star for a time step in step_detached's solve_ivp(). [rad/yr^2]

        dOmega_pri : float
            The change in rotation rate of the primary (more evolved)
        star for a time step in step_detached's solve_ivp(). [rad/yr^2]

    """

    m1 = primary.latest["mass"]
    R1 = primary.latest["R"]
    omega1 = primary.latest["omega"]
    MOI_1 = primary.latest["inertia"]

    m2 = secondary.latest["mass"]
    R2 = secondary.latest["R"]
    omega2 = secondary.latest["omega"]
    MOI_2 = secondary.latest["inertia"]

    # Torque from Rappaport, Verbunt, and Joss 1983, ApJ, 275, 713
    # The torque is eq.36 of Rapport+1983, with γ = 4. Torque units
    # converted from cgs units to [Msol], [Rsol], [yr] as all stellar
    # parameters are given in units of [Msol], [Rsol], [yr] and so that
    # dOmega_mb/dt is in units of [yr^-2].
    dOmega_mb_sec = (
        -3.8e-30 * (const.rsol**2 / const.secyer)
        * m2
        * R2**4
        * omega2**3
        / MOI_2
        * np.clip((1.5 - m2) / (1.5 - 1.3), 0, 1)
    )
    dOmega_mb_pri = (
        -3.8e-30 * (const.rsol**2 / const.secyer)
        * m1
        * R1**4
        * omega1**3
        / MOI_1
        * np.clip((1.5 - m1) / (1.5 - 1.3), 0, 1)
    )
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

    return dOmega_mb_sec, dOmega_mb_pri

def M15_braking(primary, secondary, verbose = False):
    """
        Calculates the change in spin of each star due
    to magnetic braking, according to:

        Matt et al. 2015, ApJ, 799, L23

        Free parameters calibrated to match the Solar rotation
    rate by Gossage et al. 2021, ApJ, 912, 65.

    Parameters
    ----------
        primary : SingleStar object
            A single star object, representing the primary (more evolved) star
            in the binary and containing its properties.

        secondary : SingleStar object
            A single star object, representing the secondary (less evolved) star
            in the binary and containing its properties.

        verbose : bool
            True if we want to print stuff.

    Returns
    -------
        dOmega_sec : float
            The change in rotation rate of the secondary (less evolved)
        star for a time step in step_detached's solve_ivp(). [rad/yr^2]

        dOmega_pri : float
            The change in rotation rate of the primary (more evolved)
        star for a time step in step_detached's solve_ivp(). [rad/yr^2]

    """

    m1 = primary.latest["mass"]
    R1 = primary.latest["R"]
    omega1 = primary.latest["omega"]
    MOI_1 = primary.latest["inertia"]
    tau_conv_1 = primary.latest["conv_env_turnover_time_l_b"]

    m2 = secondary.latest["mass"]
    R2 = secondary.latest["R"]
    omega2 = secondary.latest["omega"]
    MOI_2 = secondary.latest["inertia"]
    tau_conv_2 = secondary.latest["conv_env_turnover_time_l_b"]

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

    Prot_pri = 2 * np.pi / omega1    # [yr]
    Rossby_number_pri = Prot_pri / tau_conv_1
    Prot_sec = 2 * np.pi / omega2    # [yr]
    Rossby_number_sec = Prot_sec / tau_conv_2

    # critical rotation rate in rad/yr
    Omega_crit_pri = np.sqrt(
        const.standard_cgrav * m1 * const.msol
        / ((R1 * const.rsol) ** 3)) * const.secyer
    Omega_crit_sec = np.sqrt(
        const.standard_cgrav * m2 * const.msol
        / ((R2 * const.rsol) ** 3)) * const.secyer

    # omega/omega_c
    wdivwc_pri = omega1 / Omega_crit_pri
    wdivwc_sec = omega2 / Omega_crit_sec

    gamma_pri = (1 + (wdivwc_pri / 0.072)**2)**0.5
    T0_pri = K * R1**3.1 * m1**0.5 * gamma_pri**(-2 * 0.22)
    gamma_sec = (1 + (wdivwc_sec / 0.072)**2)**0.5
    T0_sec = K * R2**3.1 * m2**0.5 * gamma_sec**(-2 * 0.22)

    if (Rossby_number_sec < 0.14):
        dOmega_mb_sec = (
            T0_sec * (chi**2.6) * (omega2 / omega_sol) / MOI_2
            * np.clip((1.5 - m2) / (1.5 - 1.3), 0, 1)
        )
    else:
        dOmega_mb_sec = (
            T0_sec * ((tau_conv_2/tau_conv_sol)**2.6)
                * ((omega2/omega_sol)**(2.6 + 1)) / MOI_2
            * np.clip((1.5 - m2) / (1.5 - 1.3), 0, 1)
        )

    if (Rossby_number_pri < 0.14):
        dOmega_mb_pri = (
            T0_pri * (chi**2.6) * (omega1 / 2.6e-6) / MOI_1
            * np.clip((1.5 - m1) / (1.5 - 1.3), 0, 1)
        )
    else:
        dOmega_mb_pri = (
            T0_pri * ((tau_conv_1/tau_conv_sol)**2.6)
                * ((omega1/omega_sol)**(2.6 + 1)) / MOI_1
            * np.clip((1.5 - m1) / (1.5 - 1.3), 0, 1)
        )

    return dOmega_mb_sec, dOmega_mb_pri

def G18_braking(primary, secondary, verbose = False):
    """
        Calculates the change in spin of each star due
    to magnetic braking, according to:

        Garraffo et al. 2018, ApJ, 862, 90

        Free parameters calibrated to match the Solar rotation
    rate by Gossage et al. 2021, ApJ, 912, 65.

    Parameters
    ----------
        primary : SingleStar object
            A single star object, representing the primary (more evolved) star
            in the binary and containing its properties.

        secondary : SingleStar object
            A single star object, representing the secondary (less evolved) star
            in the binary and containing its properties.

        verbose : bool
            True if we want to print stuff.

    Returns
    -------
        dOmega_sec : float
            The change in rotation rate of the secondary (less evolved)
        star for a time step in step_detached's solve_ivp(). [rad/yr^2]

        dOmega_pri : float
            The change in rotation rate of the primary (more evolved)
        star for a time step in step_detached's solve_ivp(). [rad/yr^2]

    """

    m1 = primary.latest["mass"]
    omega1 = primary.latest["omega"]
    MOI_1 = primary.latest["inertia"]
    tau_conv_1 = primary.latest["conv_env_turnover_time_l_b"]

    m2 = secondary.latest["mass"]
    omega2 = secondary.latest["omega"]
    MOI_2 = secondary.latest["inertia"]
    tau_conv_2 = secondary.latest["conv_env_turnover_time_l_b"]

    # Torque prescription from Garraffo et al. 2018, ApJ, 862, 90
    # a = 0.03
    # b = 0.5
    # [g cm^2] -> [Msol Rsol^2]
    c = -3e41 / (const.msol * const.rsol**2)
    # Above are as calibrated in Gossage et al. 2021, ApJ, 912, 65

    Prot_pri = 2 * np.pi / omega1            # [yr]
    Rossby_number_pri = Prot_pri / tau_conv_1
    Prot_sec = 2 * np.pi / omega2            # [yr]
    Rossby_number_sec = Prot_sec / tau_conv_2

    n_pri = (0.03 / Rossby_number_pri) + 0.5 * Rossby_number_pri + 1.0
    n_sec = (0.03 / Rossby_number_sec) + 0.5 * Rossby_number_sec + 1.0

    Qn_pri = 4.05 * np.exp(-1.4 * n_pri)
    Qn_sec = 4.05 * np.exp(-1.4 * n_sec)

    dOmega_mb_sec = (
            c * omega2**3 * tau_conv_2 * Qn_sec / MOI_2
            * np.clip((1.5 - m2) / (1.5 - 1.3), 0, 1)
    )

    dOmega_mb_pri = (
            c * omega1**3 * tau_conv_1 * Qn_pri / MOI_1
            * np.clip((1.5 - m1) / (1.5 - 1.3), 0, 1)
    )

    return dOmega_mb_sec, dOmega_mb_pri

def CARB_braking(primary, secondary, verbose = False):
    """
        Calculates the change in spin of each star due
    to magnetic braking, according to:

        Van & Ivanova 2019, ApJ, 886, L31

    Parameters
    ----------
        primary : SingleStar object
            A single star object, representing the primary (more evolved) star
            in the binary and containing its properties.

        secondary : SingleStar object
            A single star object, representing the secondary (less evolved) star
            in the binary and containing its properties.

        verbose : bool
            True if we want to print stuff.

    Returns
    -------
        dOmega_sec : float
            The change in rotation rate of the secondary (less evolved)
        star for a time step in step_detached's solve_ivp(). [rad/yr^2]

        dOmega_pri : float
            The change in rotation rate of the primary (more evolved)
        star for a time step in step_detached's solve_ivp(). [rad/yr^2]

    """

    m1 = primary.latest["mass"]
    R1 = primary.latest["R"]
    omega1 = primary.latest["omega"]
    MOI_1 = primary.latest["inertia"]
    tau_conv_1 = primary.latest["conv_env_turnover_time_l_b"]
    mdot_1 = primary.latest["mdot"]

    m2 = secondary.latest["mass"]
    R2 = secondary.latest["R"]
    omega2 = secondary.latest["omega"]
    MOI_2 = secondary.latest["inertia"]
    tau_conv_2 = secondary.latest["conv_env_turnover_time_l_b"]
    mdot_2 = secondary.latest["mdot"]

    # Torque prescription from Van & Ivanova 2019, ApJ, 886, L31
    # Based on files hosted on Zenodo:
    #         https://zenodo.org/record/3647683#.Y_TfedLMKUk,
    # with units converted from [cm], [g], [s] to [Rsol], [Msol], [yr]

    # Constants as assumed in Van & Ivanova 2019, ApJ, 886, L31
    omega_sol = 3e-6 * const.secyer         # [s^-1] -> [yr^-1]
    tau_conv_sol = 2.8e6 / const.secyer     # [s] -> yr
    K2 = 0.07**2

    tau_ratio_sec = tau_conv_2 / tau_conv_sol
    tau_ratio_pri = tau_conv_1 / tau_conv_sol
    rot_ratio_sec = omega2 / omega_sol
    rot_ratio_pri = omega1 / omega_sol

    # below in units of [Rsol yr^-1]^2
    v_esc2_sec = ((2 * const.standard_cgrav * m2 / R2)
                * (const.msol * const.secyer**2 / const.rsol**3))
    v_esc2_pri = ((2 * const.standard_cgrav * m1 / R1)
                * (const.msol * const.secyer**2 / const.rsol**3))
    v_mod2_sec = v_esc2_sec + (2 * omega2**2 * R2**2) / K2
    v_mod2_pri = v_esc2_pri + (2 * omega1**2 * R1**2) / K2

    # Van & Ivanova 2019, MNRAS 483, 5595 replace the magnetic field
    # with Omega * tau_conv phenomenology. Thus, the ratios
    # (rot_ratio_* and tau_ratio_*) inherently have units of Gauss
    # [cm^-0.5 g^0.5 s^-1] that needs to be converted to [Rsol],
    # [Msol], [yr]. VI2019 assume the solar magnetic field strength is
    # on average 1 Gauss.
    if (abs(mdot_2) > 0):
        R_alfven_div_R3_sec = (
            R2**4 * rot_ratio_sec**4 * tau_ratio_sec**4
            / (mdot_2**2 * v_mod2_sec)
            * (const.rsol**2 * const.secyer / const.msol**2))
    else:
        R_alfven_div_R3_sec = 0.0

    if (abs(mdot_1) > 0):
        R_alfven_div_R3_pri = (
            R1**4 * rot_ratio_pri**4 * tau_ratio_pri**4
            / (mdot_1**2 * v_mod2_pri)
            * (const.rsol**2 * const.secyer / const.msol**2))
    else:
        R_alfven_div_R3_pri = 0.0

    # Alfven radius in [Rsol]
    R_alfven_sec = R2 * R_alfven_div_R3_sec**(1./3.)
    R_alfven_pri = R1 * R_alfven_div_R3_pri**(1./3.)

    dOmega_mb_sec = (
        (2./3.) * omega2 * mdot_2 * R_alfven_sec**2 / MOI_2
        * np.clip((1.5 - m2) / (1.5 - 1.3), 0, 1)
    )

    dOmega_mb_pri = (
        (2./3.) * omega1 * mdot_1 * R_alfven_pri**2 / MOI_1
        * np.clip((1.5 - m2) / (1.5 - 1.3), 0, 1)
    )

    return -abs(dOmega_mb_sec), -abs(dOmega_mb_pri)
