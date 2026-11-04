"""Orbit of a binary right after a supernova with a natal kick.

Ported from StepSN.orbital_kick, StepSN.get_combined_tilt and SNCheck on
main. The functions take and return plain numbers; updating the binary and
star objects is done by the SN step.

Geometry
--------
We work in a right-handed coordinate system. The exploding star lies on the
origin. The companion lies on the negative X axis at rest. The relative
velocity of the exploding star with respect to the companion lies in the
X-Y plane, with vY > 0. The orbital angular momentum vector is in the Z
direction, which completes the right-handed coordinate system.

psi
    The angle in the orbital plane between the X axis and the pre-SN
    relative velocity (psi = pi/2 points in the Y direction).
theta
    The polar angle between the kick velocity and the pre-SN relative
    velocity of the exploding star with respect to the companion.
phi
    The corresponding azimuthal angle, such that phi = 0 is on the Z axis.
tilt
    The angle between the pre- and post-SN orbital angular momentum vectors.

References
----------
.. [1] Kalogera, V. 1996, ApJ, 471, 352
.. [2] Wong, T.-W., Valsecchi, F., Fragos, T., & Kalogera, V. 2012, ApJ,
    747, 111
.. [3] Willems, B., Henninger, M., Levin, T., et al. 2005, ApJ, 625, 324
.. [4] Kalogera, V. & Lorimer, D.R. 2000, ApJ, 530, 890
"""

import numpy as np
import scipy as sp

import posydon.utils.constants as const
from posydon.utils.common_functions import rotate


def orbit_after_kick(M_pre, M_post, M_companion, separation, eccentricity,
                     v_kick, cos_theta, phi, mean_anomaly, verbose=False):
    """Compute the orbit right after the supernova [1]_, [2]_.

    Parameters
    ----------
    M_pre : float
        Mass of the exploding star before the SN in Msun.
    M_post : float
        Mass of the compact object left by the SN in Msun.
    M_companion : float
        Mass of the companion in Msun.
    separation : float
        Orbital separation (semi-major axis) before the SN in Rsun.
    eccentricity : float
        Eccentricity before the SN.
    v_kick : float
        Natal kick velocity in km/s.
    cos_theta : float
        Cosine of the polar angle of the kick (see the module docstring).
    phi : float
        Azimuthal angle of the kick in radians.
    mean_anomaly : float
        Mean anomaly of the orbit at the explosion in radians.
    verbose : bool, optional
        If True, print the checks for the disruption of the binary.

    Returns
    -------
    bound : bool
        False if the SN disrupts the binary.
    A_post : float
        Orbital separation after the SN in Rsun.
    e_post : float
        Eccentricity after the SN (nan if the orbit is unbound).
    V_sys : np.ndarray
        Change of the centre-of-mass velocity of the binary in km/s.
    tilt : float
        Angle in radians between the pre- and post-SN orbital angular
        momenta, negative if the kick has a component along -Z.
    true_anomaly : float
        True anomaly of the orbit at the explosion in radians.

    The values other than `bound` are those computed by the equations; for
    a disrupted binary it is up to the caller what to store.
    """
    epre = eccentricity
    Apre = separation
    # Eq 16 in [2]: eccentric anomaly
    E_ma = sp.optimize.brentq(
        lambda x: mean_anomaly - x + epre * np.sin(x), 0, 2 * np.pi
    )
    # Eq 15 in [2]: orbital separation at the time of the explosion
    rpre = Apre * (1.0 - epre * np.cos(E_ma))

    true_anomaly = 2 * np.arctan(
        np.sqrt((1 + epre) / (1 - epre)) * np.tan(E_ma / 2)
    )

    # convert to CGS
    G = const.standard_cgrav
    M_pre = M_pre * const.Msun
    M_companion = M_companion * const.Msun
    M_post = M_post * const.Msun
    Apre = Apre * const.Rsun
    Vkick = v_kick * const.km2cm
    rpre = rpre * const.Rsun
    Mtot_pre = M_pre + M_companion
    Mtot_post = M_post + M_companion

    sin_theta = np.sqrt(1 - (cos_theta ** 2))

    # Eq 1 in [1], extended to Eq 17 in [2]: velocity of the exploding star
    # relative to the companion (from energy conservation); not necessarily
    # along the Y axis if the orbit is eccentric
    Vr = np.sqrt(G * (Mtot_pre) * (2.0 / rpre - 1.0 / Apre))

    # Eq 18 in [2]: psi, the angle between Vr and the X axis (from the
    # specific angular momentum r x Vr = sqrt(G M A (1 - e^2)));
    # sin_psi = 1 for a circular orbit
    sin_psi = np.round(
        np.sqrt(G * (Mtot_pre) * (1 - epre ** 2) * Apre)
        / (rpre * Vr), 5)
    cos_psi = np.sqrt(1 - sin_psi ** 2)
    # allow for -cos_psi (Vr in the -X, +Y quadrant)
    if E_ma > np.pi:
        cos_psi *= -1

    # suppress overflow warnings for extreme kicks that disrupt the binary
    with np.errstate(over='ignore', divide='ignore', invalid='ignore'):
        # Eq 3 in [1], extended to Eq 13 in [2]: separation after the SN
        # (from energy conservation)
        Apost = ((2.0 / rpre)
                 - (((Vkick ** 2) + (Vr ** 2) + (2 * (Vkick * cos_theta) * Vr))
                    / (G * Mtot_post))
                 ) ** -1

        # kick components in the coordinate system
        Vkx = Vkick * (sin_theta * np.sin(phi) * sin_psi + cos_theta * cos_psi)
        Vky = Vkick * (-sin_theta * np.sin(phi) * cos_psi + cos_theta * sin_psi)
        Vkz = Vkick * sin_theta * np.cos(phi)

        # Eq 4 in [1], extended to Eq 14 in [2]: eccentricity after the SN
        # (from the specific angular momentum)
        x = ((Vkz ** 2 + (Vky + Vr * sin_psi) ** 2)
             * rpre ** 2
             / (G * Mtot_post * Apost))

        # catch negative values, i.e. disrupted binaries
        if 1. - x < 0.:
            epost = np.nan
        else:
            epost = np.sqrt(1 - x)

    # centre-of-mass velocity after the SN; it is 0 before the SN in the
    # centre-of-mass frame, so VS_post - VS_pre in our working frame
    VC0x = M_pre * Vr * cos_psi / Mtot_pre
    VC0y = M_pre * Vr * sin_psi / Mtot_pre
    VC0z = 0

    VC1x = M_post * (Vkx + Vr * cos_psi) / Mtot_post
    VC1y = M_post * (Vky + Vr * sin_psi) / Mtot_post
    VC1z = M_post * Vkz / Mtot_post

    VSx = VC1x - VC0x
    VSy = VC1y - VC0y
    VSz = VC1z - VC0z

    # angle between the pre- and post-SN orbital angular momenta:
    # Lpre || Z axis, Lpost || X axis cross the post-SN velocity of the
    # compact object, cos(tilt) = Lpre . Lpost / (|Lpre| |Lpost|);
    # for epre = 0 (sin_psi = 1) this reduces to Eq 4 in [1]
    with np.errstate(over='ignore', invalid='ignore'):
        tilt = np.arccos((Vky + Vr * sin_psi)
                         / np.sqrt(Vkz ** 2 + (Vky + Vr * sin_psi) ** 2))
    # track the direction of the tilt
    if Vkz < 0:
        tilt *= -1

    bound = _orbit_is_valid(M_pre, M_companion, M_post, rpre, Apost, epost,
                            Vr, Vkick, verbose=verbose)

    V_sys = np.array([VSx / const.km2cm, VSy / const.km2cm,
                      VSz / const.km2cm])
    return bound, Apost / const.Rsun, epost, V_sys, tilt, true_anomaly


def _orbit_is_valid(M_pre, M_companion, M_post, rpre, Apost, epost, Vr,
                    Vkick, verbose=False):
    """Check that the binary is not disrupted [3]_, [4]_.

    Parameters
    ----------
    M_pre : float
        Mass of the exploding star before the SN in g.
    M_companion : float
        Companion mass in g.
    M_post : float
        Compact object mass left by the SN in g.
    rpre : float
        Orbital separation at the time of the explosion in cm. If the
        eccentricity before the SN is 0, this corresponds to Apre.
    Apost : float
        Orbital separation after the SN in cm.
    epost : float
        Eccentricity after the SN.
    Vr : float
        Velocity of the exploding star relative to the companion before the
        SN in cm/s.
    Vkick : float
        Kick velocity in cm/s.
    verbose : bool, optional
        If True, print the checks.

    Returns
    -------
    bool
        True if the binary is not disrupted.
    """
    G = const.standard_cgrav
    Mtot_pre = M_pre + M_companion
    Mtot_post = M_post + M_companion

    # machine precision (we can probably lower this number)
    err = const.SNcheck_ERR

    # SNflag1: Eq. 21 in [3] (with typo fixed), from Eq. 10 in Flannery &
    # van den Heuvel (1975), A&A, 39, 61: the post-SN orbit must pass
    # through the pre-SN positions; updated for eccentric orbits, see
    # Eq. 15 in Wong et al. (2012)
    SNflag1 = (1 - epost - rpre / Apost <= err) and (
        rpre / Apost - (1 + epost) <= err
    )

    # SNflag2: Eqs. 22-23 in [3] (see e.g. [4]). The derivation assumes a
    # circular pre-SN orbit, hence the correction for eccentric orbits.
    with np.errstate(divide='ignore', over='ignore', invalid='ignore'):
        eccentric_orbit_correction = Vr**2 * rpre / (G * Mtot_pre)
        tmp1 = 2 - Mtot_pre / Mtot_post * (Vkick / Vr - 1) ** 2\
                   * eccentric_orbit_correction
        tmp2 = 2 - Mtot_pre / Mtot_post * (Vkick / Vr + 1) ** 2\
                   * eccentric_orbit_correction

    SNflag2 = ((rpre / Apost - tmp1 < err)
               and (err > tmp2 - rpre / Apost))

    # SNflag3: e_post must not exceed 1 or be nan
    SNflag3 = not (epost >= 1.0 or np.isnan(epost))

    SNflags = [SNflag1, SNflag2, SNflag3]

    if verbose:
        print()
        print("The orbital checks are:", SNflags)
        print()
        print("1. Post-SN orbit must pass through pre-SN positions.")
        print("2. Lower and upper limits on amount of orbital "
              "contraction or expansion that can take place for a "
              "given amount of mass loss and a given magnitude of the "
              "kick velocity.")
        print("3. Checks that e_post is not larger than 1 or nan.")

    return all(SNflags)


def combined_tilt(tilt_1, tilt_2, true_anomaly_1, true_anomaly_2):
    """Get the spin-orbit tilt after two supernovae.

    Assumes that the spin was not realigned with the orbital angular
    momentum after the first SN.

    Parameters
    ----------
    tilt_1 : float
        Angle, in radians, through which the orbital plane was tilted by
        the first SN.
    tilt_2 : float
        Angle, in radians, through which the orbital plane was tilted by
        the second SN.
    true_anomaly_1 : float
        True anomaly, in radians, at the moment of the first SN.
    true_anomaly_2 : float
        True anomaly, in radians, at the moment of the second SN.

    Returns
    -------
    float
        Angle, in radians, between the spin and the orbital angular
        momentum after the second SN.
    """
    z_prime = rotate((1, 0, 0), tilt_1).dot((0, 0, 1))
    x_prime = rotate(z_prime, true_anomaly_2 - true_anomaly_1).dot((1, 0, 0))

    cos_tilt = np.dot((0, 0, 1), rotate(x_prime, tilt_2).dot(z_prime))
    return np.arccos(cos_tilt)
