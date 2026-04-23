import numpy as np


def default_spin_from_winds(a, e, primary, secondary, verbose = False):

    """
        This function calculates the change in angular rotation rates
    from wind loss.

    Parameters
    ----------
        a : float
            The current orbital separation (unused in calculation,
            and only present in debugging output). [Rsolar]

        e : float
            The current orbital eccentricity (unused in calculation,
            and only present in debugging output).

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

    R1 = primary.latest["R"]
    omega1 = primary.latest["omega"]
    MOI_1 = primary.latest["inertia"]
    mdot_1 = primary.latest["mdot"]
    Idot_1 = primary.latest["Idot"]

    R2 = secondary.latest["R"]
    omega2 = secondary.latest["omega"]
    MOI_2 = secondary.latest["inertia"]
    mdot_2 = secondary.latest["mdot"]
    Idot_2 = secondary.latest["Idot"]

    # Due to the secondary's own evolution, we have:
    # domega_spin/dt = d(Jspin/I)/dt = dJspin/dt * 1/I + Jspin*d(1/I)/dt.
    # These are the two terms calculated below.
    dOmega_spin_wind_sec = (
            2.0 / 3.0 * R2 ** 2 * omega2 * mdot_2 / MOI_2
    )
    # jshell*Mdot/I : specific angular momentum of a thin sperical shell
    # * mdot  / moment of inertia
    # Omega is in rad/yr here, and R, M in solar (Mdot solar/yr).
    dOmega_deformation_sec = np.min(
        [-omega2 * Idot_2 / MOI_2, 100]
    )
    # This is term of Jspin*d(1/I)/dt term of the domega_spin/dt. #
    # We limit its increase due to contraction to 100 [(rad/yr)/yr]
    # increase, otherwise the integrator will fail
    # (usually when we have WD formation).
    dOmega_spin_wind_pri = (
            2.0 / 3.0 * R1 ** 2 * omega1 * mdot_1 / MOI_1
    )
    # jshell*Mdot/I : specific angular momentum of a thin sperical shell
    # * mdot  / moment of inertia
    # Omega is in rad/yr here, and R, M in solar (Mdot solar/yr).
    dOmega_deformation_pri = np.min(
        [-omega1 * Idot_1 / MOI_1, 100]
    )
    #if verbose:
    #    print(
    #        "dOmega_spin_wind , dOmega_deformation = ",
    #        dOmega_spin_wind_sec,
    #        dOmega_deformation_sec,
    #        dOmega_spin_wind_pri,
    #        dOmega_deformation_pri,
    #    )
    dOmega_sec = dOmega_spin_wind_sec + dOmega_deformation_sec
    dOmega_pri = dOmega_spin_wind_pri + dOmega_deformation_pri

    #if verbose:
    #    print("a[Ro],e,Omega[rad/yr] have been =", a, e, omega2, omega1)
        #print("da,de,dOmega (all) in 1yr is = ",
        #    da, de, dOmega_sec, dOmega_pri)

    return dOmega_sec, dOmega_pri

def default_sep_from_winds(a, e, primary, secondary, verbose = False):

    """
        This function calculates the change in orbital separation from
    wind loss, as e.g., in:

        Tauris, T. M., & van den Heuvel, E. 2006,
                Compact stellar X-ray sources, 1, 623

        from simple arguments regarding the balance of orbital angular
    momentum.

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
        da_mt_tot : float
            The change in orbital separation from both stars
        for a time step in step_detached's solve_ivp(). [Rsolar/yr]

    """

    m2 = secondary.latest["mass"] # Msol
    mdot_2 = secondary.latest["mdot"] # Msol/yr

    m1 = primary.latest["mass"] # Msol
    mdot_1 = primary.latest["mdot"] # Msol/yr

    q1 = m2 / m1
    k11 = (1 / (1 + q1)) * (mdot_2 / m2)
    k21 = mdot_2 / m2
    k31 = mdot_2 / (m1 + m2)
    # This is simplified to da_mt = -a * Mdot/(M+Macc), for only (negative)
    # wind Mdot from star M.
    da_mt_sec = a * (2 * k11 - 2 * k21 + k31)

    q2 = m1 / m2
    k12 = (1 / (1 + q2)) * (mdot_1 / m1)
    k22 = mdot_1 / m1
    k32 = mdot_1 / (m1 + m2)
    da_mt_pri = a * (2 * k12 - 2 * k22 + k32)

    #if verbose:
    #    print("da_mt = ", da_mt_sec, da_mt_pri)

    da_mt_tot = da_mt_sec + da_mt_pri

    return da_mt_tot
