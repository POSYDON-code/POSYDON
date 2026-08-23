import numpy as np
from posydon.binary_evol.binarystar import BINARYPROPERTIES


def orbital_kick(binary,verbose):
    """Do the orbital kick.

    This function computes the supernova step of the binary object [1]_,
    [2]_. It checks which binary_state reached the core collapse flag,
    either CC1 or CC2, and runs the step accordingly updating the binary
    object.

    Geometry:
    We work in a right-handed coordinate system. The collapsing helium star,
    here M_he_star, lies on the origin. The companion, here M_companion,
    lies on the negative X axis at rest. The relative velocity of the M_he_star
    with respect to M_companion lies in the X-Y plane, with vY>0.
    The orbital angular momentum vector is in Z direction, which completes the
    right-handed coordinate system.

    psi:
        The angle in the orbital plane between the X axis and the pre-core
        collapse relative velocity. (psi = pi/2 points in Y direction)

    theta :
        The polar angle between the kick velocity and the pre-core collapse
        relative velocity of the M_he_star with respect to M_companion.

    phi :
        The corresponding azimuthal angle such that phi=0 is on the Z axis.

    tilt :
        The angle between pre- and post- supernova orbital angular momentum vectors

    Parameters
    ----------
    binary : object
        Binary object containing the binary properties and the two star
        objects.

    References
    ----------
    .. [1] Kalogera, V. 1996, ApJ, 471, 352

    .. [2] Wong, T.-W., Valsecchi, F., Fragos, T., & Kalogera, V. 2012, ApJ, 747, 111

    """

    # Check that the binary_state is calling correctly the SN_step
    if binary.event != "CC1" and binary.event != "CC2":
        raise ValueError("Something went wrong: invalid call of supernova step!")

    if binary.event == "CC1":
        set_kick(
            star_that_explodes = binary.star_1,
            star_that_doesnt_explode = binary.star_2,
            binary = binary,
            sigma_kick_ECSN = sigma_kick_ECSN,
            mean_kick_ECSN = mean_kick_ECSN,
            sigma_kick_CCSN_NS = sigma_kick_CCSN_NS,
            mean_kick_CCSN_NS = mean_kick_CCSN_NS,
            sigma_kick_CCSN_BH = sigma_kick_CCSN_BH,
            mean_kick_CCSN_BH = mean_kick_CCSN_BH,
            RNG = RNG,
            verbose=verbose,
        )


    elif binary.event == "CC2":
        set_kick(
            star_that_explodes = binary.star_2,
            star_that_doesnt_explode = binary.star_1,
            binary = binary,
            sigma_kick_ECSN = sigma_kick_ECSN,
            mean_kick_ECSN = mean_kick_ECSN,
            sigma_kick_CCSN_NS = sigma_kick_CCSN_NS,
            mean_kick_CCSN_NS = mean_kick_CCSN_NS,
            sigma_kick_CCSN_BH = sigma_kick_CCSN_BH,
            mean_kick_CCSN_BH = mean_kick_CCSN_BH,
            RNG = RNG,
            verbose=verbose,
        )


    # update the orbit
    if binary.state == "disrupted" or binary.state == "initially_single_star" or binary.state == "merged":
        #the binary was already disrupted before the SN

        # update the binary object which was disrupted already before the SN
        for key in BINARYPROPERTIES:
            if key not in  ('nearest_neighbour_distance','state'):
                setattr(binary, key, None)
        #binary.state = "disrupted"
        binary.event = None
        binary.separation = np.nan
        binary.eccentricity = np.nan
        binary.V_sys = np.array([0, 0, 0])
        binary.time = binary.time_history[-1]
        binary.orbital_period = np.nan
        binary.mass_transfer_case = 'None'
        binary.first_SN_already_occurred = True

    elif ((binary.star_1.state == 'massless_remnant') or
          (binary.star_2.state == 'massless_remnant')):
        # the binary should be disrupted when a massless_remnant formed
        # update the tilt
        if not binary.first_SN_already_occurred:
            binary.star_1.spin_orbit_tilt_first_SN = np.nan
            binary.star_2.spin_orbit_tilt_first_SN = np.nan
            binary.first_SN_already_occurred = True
        else:
            binary.star_1.spin_orbit_tilt_second_SN = np.nan
            binary.star_2.spin_orbit_tilt_second_SN = np.nan

        binary.state = "disrupted"
        binary.event = None
        binary.separation = np.nan
        binary.eccentricity = np.nan
        binary.V_sys = np.array([0, 0, 0])
        binary.time = binary.time_history[-1]
        binary.orbital_period = np.nan
        binary.mass_transfer_case = 'None'

    else:

        # The binary exist: flag_binary is True if the binary is not disrupted
        flag_binary = True

        # eccentricity before the SN
        epre = binary.eccentricity
        # the orbital semimajor axis is the orbital separation
        Apre = binary.separation
        # Eq 16, Wong, T.-W., Valsecchi, F., Fragos, T., & Kalogera, V. 2012, ApJ, 747, 111
        # for eccentric anomaly
        E_ma = sp.optimize.brentq(
            lambda x: mean_anomaly - x + epre * np.sin(x), 0, 2 * np.pi
        )
        # Eq 15, Wong, T.-W., Valsecchi, F., Fragos, T., & Kalogera, V. 2012, ApJ, 747, 111
        # orbital separation at the time of the exlosion
        rpre = Apre * (1.0 - epre * np.cos(E_ma))

        true_anomaly = 2 * np.arctan(
            np.sqrt((1 + epre) / (1 - epre)) * np.tan(E_ma / 2)
        )
        # load constants in CGS
        G = const.standard_cgrav

        # Convert inputs to CGS
        M_he_star = M_he_star * const.Msun
        M_companion = M_companion * const.Msun
        M_compact_object = M_compact_object * const.Msun
        Apre = Apre * const.Rsun
        Vkick = Vkick * const.km2cm
        rpre = rpre * const.Rsun
        Mtot_pre = M_he_star + M_companion
        Mtot_post = M_compact_object + M_companion

        # get useful quantity
        sin_theta = np.sqrt(1 - (cos_theta ** 2))

        # Eq 1, in Kalogera, V. 1996, ApJ, 471, 352
        # extended to Eq 17 in Wong, T.-W., Valsecchi, F., Fragos, T., & Kalogera, V. 2012, ApJ, 747, 111
        # Vr is velocity of preSN He core relative to M_companion, NOT necessarily
        # in the direction of the Y axis if eccentric
        # Eq from conservation of energy
        Vr = np.sqrt(G * (Mtot_pre) * (2.0 / rpre - 1.0 / Apre))

        # Eq 18, Wong, T.-W., Valsecchi, F., Fragos, T., & Kalogera, V. 2012, ApJ, 747, 111
        # psi is the polar angle of the position vector of the CO with respect
        # to its pre-SN orbital velocity in the companions frame. i.e. angle between Vr and X axis
        # If epre = 0, sin_psi should be 1
        # Eq from setting specific angular momentum r X Vr = sqrt(G*M*A*(1-e**2))
        sin_psi = np.round(
            np.sqrt(G * (Mtot_pre) * (1 - epre ** 2) * Apre)
            / (rpre * Vr), 5)
        cos_psi = np.sqrt(1 - sin_psi ** 2)
        # Allow for -cos_psi (Vr in the -X, +Y quadrant)
        if E_ma > np.pi: cos_psi *= -1

        # Eq 3, in Kalogera, V. 1996, ApJ, 471, 352
        # extended to Eq 13, in Wong, T.-W., Valsecchi, F., Fragos, T., & Kalogera, V. 2012, ApJ, 747, 111
        # get the orbital separation post SN
        # Eq from conservation of energy
        # Note: Suppress overflow warnings for extreme kick scenarios that lead to
        # disrupted binaries.
        with np.errstate(over='ignore', divide='ignore', invalid='ignore'):
            Apost = ((2.0 / rpre)
                    - (((Vkick ** 2) + (Vr ** 2) + (2 * (Vkick * cos_theta) * Vr)) / (G * Mtot_post))
                    ) ** -1


            # get kicks componets in the coordinate system
            Vkx = Vkick * (sin_theta * np.sin(phi) * sin_psi + cos_theta * cos_psi)
            Vky = Vkick * (-sin_theta * np.sin(phi) * cos_psi + cos_theta * sin_psi)
            Vkz = Vkick * sin_theta * np.cos(phi)


            # Eq 4, in Kalogera, V. 1996, ApJ, 471, 352
            # extended to Eq 14 in Wong, T.-W., Valsecchi, F., Fragos, T., & Kalogera, V. 2012, ApJ, 747, 111
            # get the eccentricity post SN
            # Eq from setting specific angular momentum r X Vr = sqrt(G*M*A*(1-e**2))

            x = ((Vkz ** 2 + (Vky + Vr * sin_psi)** 2)
                * rpre ** 2
                / (G * Mtot_post * Apost))

            # catch negative values, i.e. disrupted binaries
            if 1.-x < 0.:
                epost = np.nan
            else:
                epost = np.sqrt(1 - x)

        # Compute COM velocity, VS, post SN
        # VS_pre in COM frame is 0. So VS_post in COM frame is
        # VS_post - VS_pre in our working frame

        VC0x = M_he_star * Vr * cos_psi / Mtot_pre
        VC0y = M_he_star * Vr * sin_psi / Mtot_pre
        VC0z = 0

        VC1x = M_compact_object * (Vkx + Vr * cos_psi) / Mtot_post
        VC1y = M_compact_object * (Vky + Vr * sin_psi) / Mtot_post
        VC1z = M_compact_object * Vkz / Mtot_post


        VSx = VC1x - VC0x
        VSy = VC1y - VC0y
        VSz = VC1z - VC0z


        # V_sys = np.sqrt(VSx ** 2 + VSy ** 2 + VSz ** 2)

        # Calculate the angle between the pre and post-SN orbital angular momentum vectors
        # Lpre || Z axis
        # Lpost || X axis cross the post SN velocity of the compact object
        # cos(tilt) = Lpre dot Lpost / ||Lpre||||Lpost||
        # For epre=0 (sin_psi=1), reduces to Eq 4, in Kalogera, V. 1996, ApJ, 471, 352

        # Suppress overflow warnings for extreme values
        with np.errstate(over='ignore', invalid='ignore'):
            tilt = np.arccos((Vky + Vr * sin_psi) / np.sqrt( Vkz ** 2 + (Vky + Vr * sin_psi) ** 2 ))

        # Track direction of tilt
        if Vkz < 0: tilt *= -1



        # check if the binary is disrupted
        flag_binary = _SNCheck(M_he_star, M_companion, M_compact_object, rpre,
                              Apost, epost, Vr, Vkick, cos_theta,
                              verbose=self.verbose)

        # update the binary object which was bound at least before the SN
        #Check if this is the first SN
        for key in BINARYPROPERTIES:
            if key not in ['nearest_neighbour_distance','event']:
                setattr(binary, key, None)

        if flag_binary:
            # update the tilt
            if not binary.first_SN_already_occurred:
                # update the tilt
                binary.star_1.spin_orbit_tilt_first_SN = tilt
                binary.star_2.spin_orbit_tilt_first_SN = tilt
                binary.true_anomaly_first_SN = true_anomaly
                binary.first_SN_already_occurred = True
            else:
                if binary.event == 'CC2':
                    # Assume progenitor has aligned with the preSN orbital angular momentum
                    binary.star_2.spin_orbit_tilt_second_SN = tilt
                    binary.star_1.spin_orbit_tilt_second_SN = self.get_combined_tilt(
                        tilt_1 = binary.star_1.spin_orbit_tilt_first_SN,
                        tilt_2 = tilt,
                        true_anomaly_1 = binary.true_anomaly_first_SN,
                        true_anomaly_2 = true_anomaly
                        )
                    binary.true_anomaly_second_SN = true_anomaly
                elif binary.event == 'CC1':
                    # Assume progenitor has aligned with the preSN orbital angular momentum
                    binary.star_1.spin_orbit_tilt_second_SN = tilt
                    binary.star_2.spin_orbit_tilt_second_SN = self.get_combined_tilt(
                        tilt_1 = binary.star_1.spin_orbit_tilt_first_SN,
                        tilt_2 = tilt,
                        true_anomaly_1 = binary.true_anomaly_first_SN,
                        true_anomaly_2 = true_anomaly
                        )
                    binary.true_anomaly_second_SN = true_anomaly
                else:
                    raise ValueError(f"Binary is in SN step but binary state is not CC1 or CC2: {binary.state}")

            # compute new orbital period before reseting the binary properties
            binary.state = "detached"
            binary.event = None
            binary.separation = Apost / const.Rsun
            binary.eccentricity = epost
            binary.V_sys = np.array([VSx / const.km2cm, VSy / const.km2cm, VSz
                                     / const.km2cm])
            binary.time = binary.time_history[-1]
            # in future we will make the orbital period a callable property
            new_orbital_period = orbital_period_from_separation(
                binary.separation, binary.star_1.mass, binary.star_2.mass)
            binary.orbital_period = new_orbital_period
            binary.mass_transfer_case = 'None'

        else:
            # update the tilt
            if not binary.first_SN_already_occurred:
                binary.star_1.spin_orbit_tilt_first_SN = np.nan
                binary.star_2.spin_orbit_tilt_first_SN = np.nan
                binary.first_SN_already_occurred = True
            else:
                binary.star_1.spin_orbit_tilt_second_SN = np.nan
                binary.star_2.spin_orbit_tilt_second_SN = np.nan


            binary.state = "disrupted"
            binary.event = None
            binary.separation = np.nan
            binary.eccentricity = np.nan
            binary.V_sys = np.array([0, 0, 0])
            binary.time = binary.time_history[-1]
            binary.orbital_period = np.nan
            binary.mass_transfer_case = 'None'

def set_kick(
    star_that_explodes,
    star_that_doesnt_explode,
    binary,
    sigma_kick_ECSN,
    mean_kick_ECSN,
    sigma_kick_CCSN_NS,
    mean_kick_CCSN_NS,
    sigma_kick_CCSN_BH,
    mean_kick_CCSN_BH,
    RNG,
    verbose,
):
    if star_that_explodes.SN_type == "WD":
        # compute the new separaiton prior to reseting the binary prop.
        new_separation = separation_evol_wind_loss(
            star_that_explodes.mass, binary.star_1.mass_history[-1],
            star_that_doesnt_explode.mass, binary.separation)
        new_orbital_period = orbital_period_from_separation(
            new_separation, star_that_explodes.mass, star_that_doesnt_explode.mass
        )
        for key in BINARYPROPERTIES:
            if key not in ['V_sys', 'nearest_neighbour_distance', 'state']:
                setattr(binary, key, None)
            # if key is 'nearest_neighbour_distance':
            #     setattr(binary, key, ['None', 'None', 'None'])
        binary.separation = new_separation
        if binary.state != "disrupted" and binary.state != "initially_single_star" and binary.state != "merged":
            binary.state = "detached"

        binary.event = None
        b
    # generate random point in the orbit where the kick happens
    if not star_that_explodes.natal_kick_mean_anomaly is None:
        mean_anomaly = star_that_explodes.natal_kick_mean_anomaly
        # check that ONLY one value is passed and is of type float
        if not isinstance(mean_anomaly, float):
            raise ValueError("mean_anomaly must be a single float value."
                             f"\n mean_anomaly = {mean_anomaly}")
    else:
        mean_anomaly = RNG.uniform(0, 2 * np.pi)
        star_that_explodes.natal_kick_mean_anomaly = mean_anomaly
    returninary.time = binary.time_history[-1]
        binary.eccentricity = binary.eccentricity_history[-1]
        # TODO: in feature we will make the orbital period a callable
        # property
        binary.orbital_period = new_orbital_period
        binary.mass_transfer_case = 'None'
        star_that_explodes.natal_kick_velocity = 0.0
        return

    # load relevant data
    # star1 has already collapsed into a compact object, look in the
    # history of the star to find the properties before supernova
    if star_that_explodes.state_history[-1] in STAR_STATES_CC:
        M_he_star = star_that_explodes.mass_history[-1]
    else:
        raise ValueError(
            "There is no information in the evolutionary history "
            "about STAR_STATES_CC."
        )
    M_compact_object = star_that_explodes.mass
    M_companion = star_that_doesnt_explode.mass

    # check if a kick is passed, otherwise generate it
    if not star_that_explodes.natal_kick_velocity is None:
        Vkick = star_that_explodes.natal_kick_velocity
    else:
        # Draw a random orbital kick
        # Vkick is the kick velocity with components Vkx, Vky, Vkz in
        # the above coordinate system

        if star_that_explodes.SN_type == "ECSN":
            # Kick for electron-capture SN
            Vkick = generate_kick(
                star=star_that_explodes,
                sigma=sigma_kick_ECSN,
                mean=mean_kick_ECSN,
            )
        elif ((star_that_explodes.SN_type == "CCSN")
              or (star_that_explodes.SN_type == "PPISN")
              or (star_that_explodes.SN_type == "PISN")):
            if star_that_explodes.state == 'NS':
                sigma = sigma_kick_CCSN_NS
                mean = mean_kick_CCSN_NS
            elif star_that_explodes.state == 'BH':
                sigma = sigma_kick_CCSN_BH
                mean = mean_kick_CCSN_BH
            elif star_that_explodes.state == 'massless_remnant':
                # No kick on a massless object
                sigma = None
                mean = None
            else:
                raise ValueError("CCSN/PPISN/PISN only for NS/BH.")
            # Kick for core-collapse SN
            Vkick = generate_kick(
                star=star_that_explodes, sigma=sigma, mean=mean
            )

        # removing this bit because it never gets called, added 
        # Vkick=0.0 at the point where it should be

        # elif binary.star_1.SN_type == "WD":
        #     # Kick for white dwarfs (allways f_fb = 1 => Vkick = 0)
        #     Vkick = 0.0

        else:
            Vkick = star_that_explodes.natal_kick_velocity

            raise ValueError("The SN type is not ECSN neither CCSN.")

        star_that_explodes.natal_kick_velocity = Vkick

    if not star_that_explodes.natal_kick_azimuthal_angle is None:
        phi = star_that_explodes.natal_kick_azimuthal_angle
    else:
        phi = RNG.uniform(0, 2 * np.pi)
        star_that_explodes.natal_kick_azimuthal_angle = phi

    if not star_that_explodesstar_1.natal_kick_polar_angle is None:
        cos_theta = np.cos(star_that_explodes.natal_kick_polar_angle)
    else:
        cos_theta = RNG.uniform(-1, 1)
        star_that_explodes.natal_kick_polar_angle = np.arccos(cos_theta)

    # generate random point in the orbit where the kick happens
    if not star_that_explodes.natal_kick_mean_anomaly is None:
        mean_anomaly = star_that_explodes.natal_kick_mean_anomaly
        # check that ONLY one value is passed and is of type float
        if not isinstance(mean_anomaly, float):
            raise ValueError("mean_anomaly must be a single float value."
                             f"\n mean_anomaly = {mean_anomaly}")
    else:
        mean_anomaly = RNG.uniform(0, 2 * np.pi)
        star_that_explodes.natal_kick_mean_anomaly = mean_anomaly
    return

def generate_kick(
    star,
    sigma,
    mean,
    kick_normalisation,
    kick_prescription):
    """Draw a kick from a Maxwellian distribution.

    We follow Hobbs G., Lorimer D. R., Lyne A. G., Kramer M., 2005, MNRAS, 360, 974
    and choose sigma = 265 km/s

    We rescale the kicks by 1 - f_fb as in Eq. 21 of Fryer, C. L., Belczynski, K., Wiktorowicz,
    G., Dominik, M., Kalogera, V., & Holz, D. E. (2012), ApJ, 749(1), 91.

    Parameters
    ----------
    star : object
        Star object containing the star properties.
    sigma : float
        Velocity dispersion for maxwellian or log-normal distribution.
    mean : float
        Mean for the log-normal distribution.

    Returns
    -------
    Vkick : double
        Natal orbital kick in km/s.

    -------------------------------------------------------------------------------------

    The following natal kick prescriptions are based on the mass of the
    ejected envelope during the supernova explosion:

    "asym_ej" - reference [1]
    "linear" - reference [2]

    References
    ----------

    [1] Neutron Star Kicks by the Gravitational Tug-boat Mechanism in Asymmetric
        Supernova Explosions: Progenitor and Explosion Dependence, Janka H.T., 2017,
        ApJ 837(1):84, arXiv:1611.07562 [astro-ph.HE]

    [2] New constraints on the Bray conservation-of-momentum natal kick model from
        multiple distinct observations, Richards S. M., Eldridge J. J., Briel M. M.,
        Stevance H. F., Willcox R., 2022, arXiv e-prints, p. arXiv:2208.02407

    -------------------------------------------------------------------------------------

    The "log_normal" precription draws kicks from a log-normal distribution,
    based on Disberg P., Mandel I., 2025, arXiv e-prints, p. arXiv:2505.22102v1


    """
    if sigma is None:
        Vkick = 0.0
    else:
        norm = _get_kick_normalisation(kick_normalisation,star)
        Vkick_ej = _get_kick_velocity(kick_prescription,star, sigma=sigma, mean=mean)
        Vkick = norm * Vkick_ej

    return Vkick


def _get_kick_velocity(kick_prescription, star, sigma=None, mean=None):
    """Get the kick velocity based on the chosen prescription.

    Parameters
    ----------
    star : object
        Star object containing the star properties.
    sigma : float, optional
        Velocity dispersion for the distribution.
    mean : float, optional
        Mean for the log-normal distribution.

    Returns
    -------
    Vkick_ej : float
        Kick velocity drawn from the chosen distribution.
    """

    if kick_prescription == "maxwellian":
        # sigma==None should never be reached, since in that case Vkick=0
        # in generate_kick function
        # this is a fallback
        if sigma is None:
            sigma = 265.0
        Vkick_ej = sp.stats.maxwell.rvs(loc=0., scale=sigma, size=1, random_state=self.RNG)[0]

    elif kick_prescription == "log_normal":
        # sigma==None should never be reached, since in that case Vkick=0
        # in generate_kick function
        # this is a fallback
        if sigma is None:
            sigma = 0.68
        if mean is None:
            mean = np.exp(5.60)
        Vkick_ej = sp.stats.lognorm.rvs(s=sigma, scale=mean, size=1, random_state=self.RNG)[0]

    elif kick_prescription == "asym_ej":
        f_kin = 0.1         # Fraction of SN explosion energy that is kinetic energy of the gas
        beta = 0.1          # Fraction of ejecta mass that is neutrino heated
        epsilon = 1
        alpha_ej = 0.01
        M_NS = star.mass                                    # Neutron star mass
        M_rembar=(((3*M_NS/20 + 1)**2) - 1)/0.3             # Remnant baryonic mass
        M_ej=abs(star.mass_history[-1] - M_rembar)          # Ejecta mass

        Vkick_ej = 211*(f_kin*beta*epsilon)**(1/2)*(alpha_ej/0.1)*(M_ej/0.1)*(M_NS/1.5)**(-1)

    elif kick_prescription == "linear":
        # DAVID: WHY IS THIS star.co_core_mass_history[-1]-star.mass - star.mass AND NOT abs(star.mass_history[-1] - M_rembar)
        M_ej = star.co_core_mass_history[-1]-star.mass        # Ejecta mass
        M_rem = star.mass                                     # Neutron star mass
        alpha = 115                                           # alpha and beta are best-fit parameters
        beta = 15

        Vkick_ej = alpha * (
M_ej/M_rem) + beta

    else:
        raise ValueError('kick_prescription option not supported!')

    return Vkick_ej


def _get_kick_normalisation(kick_normalisation, star):
    """Get the kick normalisation method.

    Parameters
    ----------
    star : object
        Star object containing the star properties.
    Returns
    -------
    norm : float
        Normalisation factor for the natal kick.

    """

    if kick_normalisation == 'one_minus_fallback':
        # Normalization from Eq. 21, Fryer, C. L., Belczynski, K., Wiktorowicz,
        # G., Dominik, M., Kalogera, V., & Holz, D. E. (2012), ApJ, 749(1), 91.
        norm = (1.0 - star.f_fb)
    elif kick_normalisation == 'one_over_mass':
        if star.state == 'BH':
            norm = 1.4/star.mass
        else:
            norm = 1.0
    elif kick_normalisation == 'NS_one_minus_fallback_BH_one':
        if star.state == 'BH':
            norm = 1.
        else:
            # Normalization from Eq. 21, Fryer, C. L., Belczynski, K., Wiktorowicz,
            # G., Dominik, M., Kalogera, V., & Holz, D. E. (2012), ApJ, 749(1), 91.
            norm = (1.0 - star.f_fb)
    elif kick_normalisation == 'one':
        norm = 1.
    elif kick_normalisation == 'zero':
        norm = 0.
    else:
        raise ValueError('kick_normalisation option not supported!')

    return norm



def _SNCheck(
    M_he_star,
    M_companion,
    M_compact_object,
    rpre,
    Apost,
    epost,
    Vr,
    Vkick,
    cos_theta,
    verbose,
):
    """Check that the binary is not disrupted [1]_, [2]_.

    Parameters
    ----------
    M_he_star : double
        Helium star mass before the SN in g.
    M_companion : double
        Companion star mass in g.
    M_compact_object : double
        Compact object mass left  by the SN in g.
    rpre : double
        Oribtal separation at the time of the exlosion in cm. If the
        eccentricity pre SN is 0 this correpond to Apre.
    Apost : double
        Orbital separtion after the SN in cm.
    epost : double
        Eccentricity after the SN.
    Vr : double
        Velocity of pre-SN He core relative to M_companion, directed
        along the positive y axis in cm/s.
    Vkick : double
        Kick velocity in cm/s.
    cos_theta : double
        The cosine of the angle between pre- & post-SN orbital planes.

    Returns
    -------
    flag_binary : bool
        flag_binary is True if the binary is not disrupted.

    References
    ----------
    .. [1] Willems, B., Henninger, M., Levin, T., et al. 2005, ApJ, 625, 324
    .. [2] Kalogera, V. & Lorimer, D.R. 2000, ApJ, 530, 890

    """
    # flag_binary is True if the binary is not disrupted
    flag_binary = True
    Mtot_pre = M_he_star + M_companion
    Mtot_post = M_compact_object + M_companion

    # Define machine precision (we can probaly lower this number)
    err = const.SNcheck_ERR

    # SNflag1: Eq. 21, Willems, B., Henninger, M., Levin, T., et al. 2005, ApJ, 625, 324 (with typo fixed)
    # from Eq. 10, Flannery, B.P. & van den Heuvel, E.P.J. 1975, A&A, 39, 61
    # Continuity demands post-SN orbit to pass through preSN positions.
    # Updated to work for eccentric orbits,
    # see Eq. 15 in Wong, T.-W., Valsecchi, F., Fragos, T., & Kalogera, V. 2012, ApJ, 747, 111
    SNflag1 = (1 - epost - rpre / Apost <= err) and (
        rpre / Apost - (1 + epost) <= err
    )

    # SNflag2: Equations 22-23, Willems, B., Henninger, M., Levin, T., et al. 2005, ApJ, 625, 324
    # (see, e.g., Kalogera, V. & Lorimer, D.R. 2000, ApJ, 530, 890)
    # The derivation in the papers above assume a circular pre SN
    # orbit. Hence, need a correction for eccentric pre SN orbits:
    # Suppress divide by zero warnings for edge cases
    with np.errstate(divide='ignore', over='ignore', invalid='ignore'):
        eccentric_orbit_correction = Vr**2 * rpre / (G * Mtot_pre)
        tmp1 = 2 - Mtot_pre / Mtot_post * (Vkick / Vr - 1) ** 2\
                   * eccentric_orbit_correction
        tmp2 = 2 - Mtot_pre / Mtot_post * (Vkick / Vr + 1) ** 2\
                   * eccentric_orbit_correction

    SNflag2 = ((rpre / Apost - tmp1 < err)
               and (err > tmp2 - rpre / Apost))

    # SNflag3: check that epost does not exeed 1 or is nan
    if epost >= 1.0 or pd.isna(epost):
        SNflag3 = False
    else:
        SNflag3 = True

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

    # check if the supernova is valid and doesn't disrupt the system
    if not all(SNflags):
        flag_binary = False

    return flag_binary