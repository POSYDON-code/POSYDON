import numpy as np
import pandas as pd
import posydon.utils.constants as const

def default_tides(a, e, primary, secondary, verbose=False):

    """
        Calculates the change in orbital separation and eccentricity, 
    plus change in spin of each star due to tides, according to: 
    
        Hut, P. 1981, A&A, 99, 126

    Parameters
    ----------
        a : float
            The current orbital separation. [Rsolar]

        e : float
            The current orbital eccentricity.

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
        da : float
            The change in orbital separation for a time step in 
        step_detached's solve_ivp(). [Rsolar/yr]

        de : float
            The change in orbital eccentricity for a time step in 
        step_detached's solve_ivp().

        dOmega_sec : float
            The change in rotation rate of the secondary (less evolved) 
        star for a time step in step_detached's solve_ivp(). [rad/yr^2]
        
        dOmega_pri : float
            The change in rotation rate of the primary (more evolved) 
        star for a time step in step_detached's solve_ivp(). [rad/yr^2]

    """

    m1 = primary.latest["mass"]
    R1 = primary.latest["R"]
    m_env1 = primary.latest["mass_conv_reg_fortides"]
    dr_env1 = primary.latest["thickness_conv_reg_fortides"]
    rmid_env1 = primary.latest["radius_conv_reg_fortides"]
    omega1 = primary.latest["omega"]
    rtop_env1 = primary.latest["conv_mx1_top_r"]
    rbot_env1 = primary.latest["conv_mx1_bot_r"]
    Xsurf_1 = primary.latest["surface_h1"]
    L1 = primary.latest["L"]
    MOI_1 = primary.latest["inertia"]

    m2 = secondary.latest["mass"]
    R2 = secondary.latest["R"]
    m_env2 = secondary.latest["mass_conv_reg_fortides"]
    dr_env2 = secondary.latest["thickness_conv_reg_fortides"]
    rmid_env2 = secondary.latest["radius_conv_reg_fortides"]
    omega2 = secondary.latest["omega"]
    rtop_env2 = secondary.latest["conv_mx1_top_r"]
    rbot_env2 = secondary.latest["conv_mx1_bot_r"]
    Xsurf_2 = secondary.latest["surface_h1"]
    L2 = secondary.latest["L"]
    MOI_2 = secondary.latest["inertia"]

    q1 = m1/ m2
    q2 = m2 / m1
    # P_orb in years. From 3rd Kepler's law, transforming separation from
    # Rsolar to AU, to avoid using constants
    # TODO: we aleady have this function!
    P_orb = np.sqrt((a / const.aursun) ** 3 / (m1 + m2))
    n = 2.0 * const.pi / P_orb  # mean orbital ang. vel. in rad/year
    f1 = (
            1
            + (31 / 2) * e ** 2
            + (255 / 8) * e ** 4
            + (185 / 16) * e ** 6
            + (25 / 64) * e ** 8
    )
    f2 = 1 + (15 / 2) * e ** 2 + (45 / 8) * e ** 4 + (5 / 16) * e ** 6
    f3 = 1 + (15 / 4) * e ** 2 + (15 / 8) * e ** 4 + (5 / 64) * e ** 6
    f4 = 1 + (3 / 2) * e ** 2 + (1 / 8) * e ** 4
    f5 = 1 + 3 * e ** 2 + (3 / 8) * e ** 4

    # equilibrium timecale
    if ((pd.notna(m_env2) and m_env2 != 0.0)
            and (pd.notna(dr_env2) and dr_env2 != 0.0)
            and (pd.notna(rmid_env2) and rmid_env2 != 0.0)):
        # eq. (31) of Hurley et al. 2002, generalized for convective layers
        # not on surface too
        tau_conv_2 = 0.431 * ((m_env2 * dr_env2 * rmid_env2
                                / (3 * L2)) ** (1.0 / 3.0))
    else:
        if verbose and verbose != 1:
            print("something wrong with M_env/DR_env/Renv_middle",
                m_env2, dr_env2, rmid_env2)
        tau_conv_2 = 1.0e99
    if ((pd.notna(m_env1) and m_env1 != 0.0)
            and (pd.notna(dr_env1) and dr_env1 != 0.0)
            and (pd.notna(rmid_env1) and rmid_env1 != 0.0)):
        # eq. (31) of Hurley et al. 2002, generalized for convective layers
        # not on surface too
        tau_conv_1 = 0.431 * ((m_env1 * dr_env1 * rmid_env1
                                / (3 * L1)) ** (1.0/3.0))
    else:
        if verbose:
            print("something wrong with M_env/DR_env/Renv_middle",
                m_env1, dr_env1, rmid_env1)
        tau_conv_1 = 1.0e99
    P_spin_sec = 2 * np.pi / omega2
    P_spin_pri = 2 * np.pi / omega1
    P_tid_sec = np.abs(1 / (1 / P_orb - 1 / P_spin_sec))
    P_tid_pri = np.abs(1 / (1 / P_orb - 1 / P_spin_pri))
    f_conv_sec = np.min([1, (P_tid_sec / (2 * tau_conv_2)) ** 2])
    f_conv_pri = np.min([1, (P_tid_pri / (2 * tau_conv_1)) ** 2])
    F_tid = 1.  # not 50 as before
    kT_conv_sec = (
            (2. / 21) * (f_conv_sec / tau_conv_2) * (m_env2 / m2)
    )  # eq. (30) of Hurley et al. 2002
    kT_conv_pri = (
            (2. / 21) * (f_conv_pri / tau_conv_1) * (m_env1 / m1)
    )
    if kT_conv_sec is None or not np.isfinite(kT_conv_sec):
        kT_conv_sec = 0.0
    if kT_conv_pri is None or not np.isfinite(kT_conv_pri):
        kT_conv_pri = 0.0

        if verbose:
            print("kT_conv_sec is", kT_conv_sec, ", set to 0.")
            print("kT_conv_pri is", kT_conv_pri, ", set to 0.")
    # this is the 1/timescale of all d/dt calculted below in yr^-1

    if verbose:
        print(
            "Equilibrium tides in deep convective envelope",
            m_env2,
            dr_env2,
            rmid_env2,
            R2,
            m2,
            m_env1,
            dr_env1,
            rmid_env1,
            R1,
            m1
        )
        print("convective tiimescales and efficiencies",
            tau_conv_2, P_orb, P_spin_sec, P_tid_sec,
            f_conv_sec,
            F_tid,
            tau_conv_1, P_orb, P_spin_pri, P_tid_pri,
            f_conv_pri,
            F_tid,
            )

    # dynamical timecale
    F_tid = 1
    # E2 = 1.592e-9*M**(2.84) # eq. (43) of Hurley et al. 2002. Deprecated
    R_conv_sec = rtop_env2 - rbot_env2
    R_conv_pri = rtop_env1 - rbot_env1  # convective core
    # R_conv = conv_mx1_top_r  # convective core
    if (R_conv_sec > R2 or R_conv_sec <= 0.0
            or rbot_env2 / R2 > 0.1):
        # R_conv = 0.5*R
        # if verbose:
        #     print(
        #         "R_conv of the convective core is not behaving well or "
        #         "we are not calculating the convective core, we make it "
        #         "equal to half of Rstar",
        #         R_conv,
        #         R,
        #         conv_mx1_top_r,
        #         conv_mx1_bot_r,
        #     )
        # we switch to Zahn+1975 calculation of E2
        E21 = 1.592e-9 * m2 ** (2.84)
    else:
        if R2 <= 0:
            E21 = 0
        elif Xsurf_2 > 0.4:
            E21 = 10.0 ** (-0.42) * (R_conv_sec / R2) ** (
                7.5
            )  # eq. (9) of Qin et al. 2018, 616, A28
        elif Xsurf_2 <= 0.4:  # "HeStar":
            E21 = 10.0 ** (-0.93) * (R_conv_sec / R2) ** (
                6.7
            )  # eq. (9) of Qin et al. 2018, 616, A28
        else:  # in principle we should not go here
            E21 = 1.592e-9 * m2 ** (
                2.84
            )  # eq. (43) of Hurley et al. 2002 from Zahn+1975 Depreciated
    # kT = 1.9782e4 * np.sqrt(M * R**2 / a**5) * (1 + q)**(5. / 6) * E2
    # eq. (42) of Hurley et al. 2002. Depreciated

    if (R_conv_pri > R1 or R_conv_pri <= 0.0
            or rbot_env1 / R1 > 0.1):
        E22 = 1.592e-9 * m1 ** (2.84)
        if verbose:
            print(
                "R_conv of the convective core is not behaving well or we "
                "are not calculating the convective core, we switch to "
                "Zahn+1975 calculation of E2",
                R_conv_sec,
                R2,
                rtop_env2,
                rbot_env2,
                E21,
                R_conv_pri,
                R1,
                rtop_env1,
                rbot_env1,
                E22
            )
    else:
        if R1 <= 0:
            E22 = 0
        elif Xsurf_1 > 0.4:
            E22 = 10.0 ** (-0.42) * (R_conv_pri / R1) ** (
                7.5
            )  # eq. (9) of Qin et al. 2018, 616, A28
        elif Xsurf_1 <= 0.4:  # "HeStar":
            E22 = 10.0 ** (-0.93) * (R_conv_pri / R1) ** (
                6.7
            )  # eq. (9) of Qin et al. 2018, 616, A28
        else:  # in principle we should not go here
            E22 = 1.592e-9 * m1 ** (
                2.84
            )
    kT_rad_sec = (
        np.sqrt(const.standard_cgrav * (m2 * const.msol)
                * (R2 * const.rsol)**2 / (a * const.rsol)**5)
        * (1 + q1) ** (5.0 / 6)
        * E21
        * const.secyer)
    kT_rad_pri = (
        np.sqrt(const.standard_cgrav * (m1 * const.msol)
                * (R1 * const.rsol)**2 / (a * const.rsol)**5)
        * (1 + q2) ** (5.0 / 6)
        * E22
        * const.secyer)
    # this is the 1/timescale of all d/dt calculted below in yr^-1
    if verbose:
        print(
            "Dynamical tides in radiative envelope",
            rtop_env2,
            rbot_env2,
            R_conv_sec,
            E21,
            rtop_env1,
            rbot_env1,
            R_conv_pri,
            E22,
            F_tid
        )
    kT_sec = max(kT_conv_sec, kT_rad_sec)
    kT_pri = max(kT_conv_pri, kT_rad_pri)
    if verbose:
        print("kT_conv/rad of tides is ", kT_conv_sec, kT_rad_sec,
            kT_conv_pri, kT_rad_pri, "in 1/yr, and we picked the ",
            kT_sec, kT_pri)

    da_tides_sec = (
            -6
            * F_tid
            * kT_sec
            * q1
            * (1 + q1)
            * (R2 / a) ** 8
            * (a / (1 - e ** 2) ** (15 / 2))
            * (f1 - (1 - e ** 2) ** (3 / 2) * f2 * omega2 / n)
    )  # eq. (9) of Hut 1981, 99, 126

    da_tides_pri = (
            -6
            * F_tid
            * kT_pri
            * q2
            * (1 + q2)
            * (R1 / a) ** 8
            * (a / (1 - e ** 2) ** (15 / 2))
            * (f1 - (1 - e ** 2) ** (3 / 2) * f2 * omega1 / n)
    )
    de_tides_sec = (
            -27
            * F_tid
            * kT_sec
            * q1
            * (1 + q1)
            * (R2 / a) ** 8
            * (e / (1 - e ** 2) ** (13 / 2))
            * (f3 - (11 / 18) * (1 - e ** 2) ** (3 / 2) * f4 * omega2/n)
    )  # eq. (10) of Hut 1981, 99, 126

    de_tides_pri = (
            -27
            * F_tid
            * kT_pri
            * q2
            * (1 + q2)
            * (R1 / a) ** 8
            * (e / (1 - e ** 2) ** (13 / 2))
            * (f3 - (11 / 18) * (1 - e ** 2) ** (3 / 2) * f4 * omega1/n)
    )

    dOmega_tides_sec = (
            (3 * F_tid * kT_sec * q1 ** 2 * m2 * R2 ** 2 / MOI_2)
            * (R2 / a) ** 6
            * n
            / (1 - e ** 2) ** 6
            * (f2 - (1 - e ** 2) ** (3 / 2) * f5 * omega2 / n)
    )  # eq. (11) of Hut 1981, 99, 126
    dOmega_tides_pri = (
            (3 * F_tid * kT_pri * q2 ** 2 * m1 * R1 ** 2 / MOI_1)
            * (R1 / a) ** 6
            * n
            / (1 - e ** 2) ** 6
            * (f2 - (1 - e ** 2) ** (3 / 2) * f5 * omega1 / n)
    )
    if verbose:
        print("da,de,dOmega_tides = ",
            da_tides_sec, de_tides_sec, dOmega_tides_sec,
            da_tides_pri, de_tides_pri, dOmega_tides_pri)
        
    da = da_tides_sec + da_tides_pri
    de = de_tides_sec + de_tides_pri
    dOmega_sec = dOmega_tides_sec
    dOmega_pri = dOmega_tides_pri

    return da, de, dOmega_sec, dOmega_pri