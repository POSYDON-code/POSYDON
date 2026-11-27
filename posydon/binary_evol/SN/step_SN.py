"""Supernova step.

This step models the end of life of stars. Applied to a binary in which a
star reaches core collapse (event CC1 or CC2), it determines the fate of that
star, from the formation of a white dwarf to an electron-capture supernova,
an Fe core-collapse supernova, a pulsational pair-instability supernova or a
pair-instability supernova, together with its remnant and natal kick, and
updates the orbit of the binary.
"""

__authors__ = [
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Jaime Roman Garza <Jaime.Roman@etu.unige.ch>",
    "Emmanouil Zapartas <ezapartas@gmail.com>",
    "Konstantinos Kovlakas <Konstantinos.Kovlakas@unige.ch>",
    "Devina Misra <devina.misra@unige.ch>",
    "Zepei Xing <Zepei.Xing@unige.ch>",
    "Jeffrey Andrews <jeffrey.andrews@northwestern.edu>",
    "Tassos Fragos <Anastasios.Fragkos@unige.ch>",
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>",
    "Max Briel <max.briel@unige.ch>",
    "Seth Gossage <seth.gossage@northwestern.edu>",
    "Dimitris Souropanis <dsouropanis@ia.forth.gr>",
    "David R. Aguilera-Dena <daguilera@astro.unam.mx>"
]

__credits__ = [
    "Michael Zevin <michael.zevin@ligo.org>",
    "Chase Kimball <charles.kimball@ligo.org",
    "Sam Imperato <samuelimperato2022@u.northwestern.edu>",
]

import numpy as np

from posydon.binary_evol.binarystar import BINARYPROPERTIES
from posydon.binary_evol.flow_chart import (
    STAR_STATES_C_DEPLETED,
    STAR_STATES_CC,
    STAR_STATES_CO,
)
from posydon.binary_evol.SN.explodability.CCSN_explosion_criteria import (
    get_CCSN_prescription,
    replace_deprecated_mechanism,
    replace_unavailable_mechanism,
)
from posydon.binary_evol.SN.explodability.ECSN_prescriptions import (
    get_ECSN_prescription,
)
from posydon.binary_evol.SN.explodability.PISN_prescriptions import (
    get_PISN_prescription,
)
from posydon.binary_evol.SN.interp_collapse import InterpCollapse
from posydon.binary_evol.SN.kicks.kick_prescriptions import (
    get_kick_prescription,
    replace_deprecated_kick_normalisation,
)
from posydon.binary_evol.SN.kicks.orbit import combined_tilt, orbit_after_kick
from posydon.binary_evol.SN.remnant_mass.remnant_finaliser import (
    Remnant_Finaliser,
)
from posydon.binary_evol.SN.remnant_mass.remnant_mass_prescriptions import (
    get_remnant_mass_prescription,
)
from posydon.grids.SN_MODELS import DEFAULT_SN_MODEL
from posydon.utils.common_functions import (
    inspiral_timescale_from_separation,
    is_None_or_nan,
    orbital_period_from_separation,
    separation_evol_wind_loss,
    set_binary_to_failed,
)
from posydon.utils.posydonerror import ModelError

# binary states without an orbit before the SN
NO_ORBIT_STATES = ("disrupted", "initially_single_star", "merged")

def reset_binary_properties(binary,keep):
    """Set every BINARYPROPERTIES key that is not in keep to None."""
    for key in BINARYPROPERTIES:
        if key not in keep:
            setattr(binary,key,None)


class StepSN(object):
    """The supernova step in POSYDON.

    For the collapsing star, the step:

    1. classifies the collapse as a WD, an ECSN, a PPISN or PISN candidate,
       or a CCSN (options ECSN and PISN);
    2. applies the explosion criterion of the mechanism or, with
       use_interp_values=True, uses the outcome interpolated from the grids
       when it is consistent with the star;
    3. computes the baryonic remnant mass and the fallback fraction;
    4. turns them into the remnant: a WD, a NS, a BH or no remnant;
    5. draws the natal kick and updates the orbit of the binary.

    Keyword Arguments
    -----------------
    mechanism : str
        Mechanism of the core collapse, which predicts the remnant. Available
        options:

        * 'Fryer+12-rapid' : the rapid supernova engine of [1]_.
        * 'Fryer+12-delayed' : the delayed supernova engine of [1]_.
        * 'direct' : the collapsing mass (the He core, or the whole star with
          conserve_hydrogen_envelope=True) forms the remnant.
        * 'direct_he_core' : deprecated, the same as 'direct' with
          conserve_hydrogen_envelope=False.
        * 'Sukhbold+16-engine' : the outcomes of the models of [2]_.
        * 'Couch+20-engine' : the outcomes of the models of [6]_.
        * 'Patton&Sukhbold20-engine' : the CO core of the star at He
          depletion is mapped onto the models of [5]_, whose M4 and mu4 are
          judged by the criterion of [9]_.
        * 'Maltsev+25-engine' : not available; 'Patton&Sukhbold20-engine'
          with engine 'N20' is used instead, with a warning.

    engine : str
        Engine of the mechanism:

        * 'Sukhbold+16-engine' : 'N20', 'S19.8', 'W15', 'W18' or 'W20'.
        * 'Patton&Sukhbold20-engine' : 'N20', 'S19.8', 'W15', 'W18', 'W20' or
          'Ertl2020', which set the parameters of the criterion of [9]_.
        * 'Couch+20-engine' : the turbulence strength, as a string: '1.0',
          '1.2', '1.23', '1.25', '1.27', '1.3' or '1.4'.

    PISN : str, float or None
        Prescription for pair instability. Available options:

        * 'Marchant+19' : He-core mass ranges and remnant masses of [3]_.
        * 'Hendriks+23' : CO-core mass range and mass lost in the pulses
          of [8]_.
        * a number : the remnant mass of stars with a He core above this
          mass (Msun) is capped at it.
        * None : no pair instability.

    PISN_CO_shift : float
        Shift (Msun) of the CO-core masses of 'Hendriks+23'.
    PPI_extra_mass_loss : float
        Extra mass (Msun) lost in the pulses of 'Hendriks+23'.
    ECSN : str or None
        Prescription for electron-capture supernovae. Available options:

        * 'Tauris+15' : limits on the CO core mass from [4]_.
        * 'Podsiadlowski+04' : limits on the He core mass from [7]_.
        * 'No_ECSN' or None : no electron-capture supernovae.

    conserve_hydrogen_envelope : bool
        If True, the whole star collapses, including its H envelope;
        otherwise only its He core collapses.
    conserve_hydrogen_PPI : bool
        If True, the pulses of 'Hendriks+23' remove mass from the whole star;
        otherwise from its He core.
    max_neutrino_mass_loss : float
        Maximum baryonic mass (Msun) lost to neutrinos during the collapse.
    max_NS_mass : float
        Maximum NS mass (Msun): remnants with a gravitational mass at or
        above it are BHs.
    use_interp_values : bool
        If True, use the outcome of the collapse interpolated from the
        post-processed grids and stored on the star, when it is consistent
        with the star; otherwise compute it. The settings must then match one
        of the SN models of posydon.grids.SN_MODELS.
    use_profiles : bool
        If True, collapse the stellar profile, when the star has one, to get
        the spin and disk of a BH and the ejecta.
    use_core_masses : bool
        If True, use the core masses at C depletion (the classical population
        synthesis treatment), with the BH spin from the angular momentum of
        the star.
    allow_spin_None : bool
        If True, do not determine the spin during the collapse, while setting
        the other values as with use_core_masses (used to avoid jumps in the
        spin for interpolator training because of missing profiles).
    approx_at_he_depletion : bool
        Only for 'Patton&Sukhbold20-engine'. If True, approximate the CO core
        at He depletion from the CO core mass at C depletion; otherwise take
        it from the star or its history.
    kick : bool
        If False, the natal kick velocity is zero. The mass lost in the
        supernova still changes the orbit.
    kick_prescription : str
        Distribution of the natal kick velocity: 'maxwellian' [10]_,
        'log_normal' [11]_, 'asym_ej' [12]_ or 'linear' [13]_.
    kick_normalisation : str
        Factor multiplying the kick velocity:

        * 'one_minus_fallback' : 1 - f_fb
        * 'one_over_mass' : 1.4 Msun / m_BH for BHs, 1 for NSs
        * 'NS_one_minus_fallback_BH_one' : 1 - f_fb for NSs, 1 for BHs
        * 'one' : 1
        * 'zero' : 0

        The deprecated values 'asym_ej' and 'linear' set kick_prescription
        instead, with a factor of 1.
    sigma_kick_CCSN_NS, sigma_kick_CCSN_BH, sigma_kick_ECSN : float or None
        Dispersion (Maxwellian) or shape (log-normal) of the kick distribution
        of NSs and BHs from core-collapse supernovae, and of the remnants of
        electron-capture supernovae (km/s for the Maxwellian). None means no
        kick.
    mean_kick_CCSN_NS, mean_kick_CCSN_BH, mean_kick_ECSN : float or None
        Median (km/s) of the log-normal kick distribution; if None,
        exp(5.60) km/s.
    no_kick_direct_collapse : bool
        If True, BHs formed by direct collapse (f_fb = 1) get no natal kick.
    RNG : numpy.random.Generator or None
        Random number generator; a new one is created if None.
    verbose : bool
        If True, print messages about each stage of the collapse.

    References
    ----------
    .. [1] Fryer, C. L., Belczynski, K., Wiktorowicz, G., Dominik, M.,
        Kalogera, V., & Holz, D. E. (2012). Compact remnant mass function:
        dependence on the explosion mechanism and metallicity. ApJ, 749, 91.

    .. [2] Sukhbold, T., Ertl, T., Woosley, S. E., Brown, J. M., & Janka, H. T.
        (2016). Core-collapse supernovae from 9 to 120 solar masses based on
        neutrino-powered explosions. ApJ, 821, 38.

    .. [3] Marchant, P., Renzo, M., Farmer, R., Pappas, K. M., Taam, R. E.,
        de Mink, S. E., & Kalogera, V. (2019). Pulsational pair-instability
        supernovae in very close binaries. ApJ, 882, 36.

    .. [4] Tauris, T. M., Langer, N., & Podsiadlowski, P. (2015).
        Ultra-stripped supernovae: progenitors and fate. MNRAS, 451, 2123.

    .. [5] Patton, R. A., & Sukhbold, T. (2020). Towards a realistic explosion
        landscape for binary population synthesis. MNRAS, 499, 2803.

    .. [6] Couch, S. M., Warren, M. L., & O'Connor, E. P. (2020). Simulating
        turbulence-aided neutrino-driven core-collapse supernova explosions in
        one dimension. ApJ, 890, 127.

    .. [7] Podsiadlowski, P., Langer, N., Poelarends, A. J. T., Rappaport, S.,
        Heger, A., & Pfahl, E. (2004). The effects of binary evolution on the
        dynamics of core collapse and neutron star kicks. ApJ, 612, 1044.

    .. [8] Hendriks, D. D., et al. (2023). MNRAS, 526, 4130.

    .. [9] Ertl, T., Janka, H.-T., Woosley, S. E., Sukhbold, T., & Ugliano, M.
        (2016). A two-parameter criterion for classifying the explodability of
        massive stars by the neutrino-driven mechanism. ApJ, 818, 124.

    .. [10] Hobbs, G., Lorimer, D. R., Lyne, A. G., & Kramer, M. (2005).
        A statistical study of 233 pulsar proper motions. MNRAS, 360, 974.

    .. [11] Disberg, P., & Mandel, I. (2025). ApJ, 989, L8.

    .. [12] Janka, H.-T. (2017). Neutron star kicks by the gravitational
        tug-boat mechanism in asymmetric supernova explosions. ApJ, 837, 84.

    .. [13] Richards, S. M., Eldridge, J. J., Briel, M. M., Stevance, H. F., &
        Willcox, R. (2022). arXiv:2208.02407.
    """
    DEFAULT_KWARGS = {
        "kick": True,
        "kick_normalisation": 'one_over_mass',
        "kick_prescription": 'maxwellian',
        "sigma_kick_CCSN_NS": 265.0,
        "mean_kick_CCSN_NS": None,
        "sigma_kick_CCSN_BH": 265.0,
        "mean_kick_CCSN_BH": None,
        "sigma_kick_ECSN": 20.0,
        "mean_kick_ECSN": None,
        "no_kick_direct_collapse": False,
        # other
        "RNG": None,
        "verbose": False
    }

    DEFAULT_KWARGS.update(DEFAULT_SN_MODEL)

    def __init__(self,**kwargs):
        """Set up the stages of the step from the keyword arguments."""
        unknown_kwargs = set(kwargs) - set(self.DEFAULT_KWARGS)
        if unknown_kwargs:
            raise ValueError(
                f"Invalid arguments passed to StepSN: {unknown_kwargs}")
        config = self.DEFAULT_KWARGS.copy()
        config.update(kwargs)
        config = replace_deprecated_mechanism(config)
        config = replace_unavailable_mechanism(config)
        config = replace_deprecated_kick_normalisation(config)
        if config["RNG"] is None:
            config["RNG"] = np.random.default_rng()
        if config["max_neutrino_mass_loss"] is None:
            config["max_neutrino_mass_loss"] = 0
        for key, value in config.items():
            setattr(self, key, value)

        # also tells whether a star is above the pair-instability gap
        self._PISN_check = get_PISN_prescription(**config)
        # checked in order, the first one that is not None gives the SN type
        self._SN_type_checks = (
            get_ECSN_prescription(**config),    # WD or ECSN
            self._PISN_check,                   # PPISN or PISN candidate
            lambda star: "CCSN",                # everything else
        )
        self._explosion_criterion = get_CCSN_prescription(**config)
        self._get_remnant_mass = get_remnant_mass_prescription(**config)
        self._update_params = Remnant_Finaliser(**config)
        self._natal_kick = get_kick_prescription(**config)
        # fails here, at setup, if no SN model matches the settings
        self._interp_collapse = (InterpCollapse(**config)
                                 if config['use_interp_values'] else None)
        if self.verbose:
            print(f"StepSN: {self._SN_type_checks[0]!r}; {self._PISN_check!r};"
                  f" {self._explosion_criterion!r};"
                  f" {self._get_remnant_mass!r}; {self._natal_kick!r}")

    def __call__(self,binary):
        """Apply the supernova to the binary.

        Collapses the star of the event (star_1 for CC1, star_2 for CC2),
        updates the orbit and finalises the step.

        Parameters
        ----------
        binary : BinaryStar
            The binary, with event "CC1" or "CC2".

        Raises
        ------
        ModelError
            If the collapse fails; the binary is then set to failed.
        """
        dying_star,other_star = self._select_dying_star(binary)
        model_err = self.collapse_star(dying_star)
        if model_err is not None:
            set_binary_to_failed(binary)
            raise ModelError(model_err)
        self._reset_other_star_properties(star=other_star)
        binary.update_star_states()
        self.update_orbit(binary,dying_star,other_star)
        self.finalise_step(binary)

    @staticmethod
    def _select_dying_star(binary):
        """Return (collapsing star, companion) for the event of the binary."""
        if binary.event=="CC1":
            return binary.star_1, binary.star_2
        elif binary.event=="CC2":
            return binary.star_2, binary.star_1
        else:
            raise ValueError(
                "Invalid call of step_SN: \n"
                f"Called during binary event {binary.event}. "
                "Step_SN should only be called if binary event is CC1 or CC2"
            )

    def collapse_star(self,star):
        """Collapse the star into its remnant.

        With use_interp_values=True, the outcome interpolated from the grids
        is used when it is consistent with the star. Otherwise the outcome is
        computed: the SN type, the explosion criterion, the remnant mass and
        finally the remnant (state, mass, spin, disk and ejecta), using the
        stellar profile (use_profiles), the core masses (use_core_masses) or
        no spin (allow_spin_None), in that order of preference.

        Parameters
        ----------
        star : SingleStar
            The collapsing star.

        Returns
        -------
        str or None
            An error message if the collapse fails, else None.
        """
        if star.state not in STAR_STATES_CC:
            return f"The star cannot collapse: star state {star.state}."
        if is_None_or_nan(star.co_core_mass):
            # the classification and the remnant masses need the CO core
            raise ValueError("The CO core mass is not correct! CO core = "
                             f"{star.co_core_mass}")
        SN_type = self.collapse_outcome(star)
        if (self._interp_collapse is not None
                and self._interp_collapse(star, SN_type)):
            return None
        explosion = self.explosion_outcome(star, SN_type)
        SN_type, m_rembar, f_fb, state = self._get_remnant_mass(
            star, SN_type, explosion)
        # returns an error message if the collapse fails, else None
        return self._update_params(star, SN_type, m_rembar, f_fb, state)

    def explosion_outcome(self,star,SN_type):
        """Return the explosion outcome of a CCSN, PPISN or PISN.

        A PISN leaves no remnant, but its f_fb comes from the CCSN
        prescription, which needs the explosion outcome.
        Above the pair-instability gap the core collapses directly into a BH
        ("BH"), whatever the mechanism, as in Marchant et al. (2019) and
        Hendriks et al. (2023).
        """
        if SN_type not in ("CCSN", "PPISN", "PISN"):
            return None
        if self._PISN_check.is_above_gap(star):
            if self.verbose:
                print("StepSN: above the pair-instability gap, direct "
                      "collapse into a BH")
            return "BH"
        return self._explosion_criterion(star)

    def collapse_outcome(self,star):
        """Return the SN type of the star, from the classification checks.

        "WD" or "ECSN" from the ECSN check, "PPISN" or "PISN" for
        pair-instability candidates from the PISN check, else "CCSN". The
        final PPISN/PISN label is settled by the remnant-mass stage.
        """
        outcomes = (check(star) for check in self._SN_type_checks)
        return next(outcome for outcome in outcomes if outcome is not None)

    def _reset_other_star_properties(self, star):
        """Reset the properties of the star that is not being collapsed."""
        star.lg_mdot = None
        star.lg_system_mdot = None

    def update_orbit(self,binary,dying_star,other_star):
        """Update the orbit of the binary after the collapse of dying_star.

        A WD loses its envelope gradually and gets no kick, so the orbit slowly
        widens. A SN loses its mass instantaneously, which changes the orbit
        even without a natal kick, so every SN goes through the kick
        equations; with kick=False the natal kick is zero.
        """
        if dying_star.SN_type == "WD":
            self._update_orbit_after_WD(binary,dying_star,other_star)
        else:
            self._update_orbit_with_kick(binary,dying_star,other_star)

    def _update_orbit_with_kick(self,binary,dying_star,other_star):
        """Update the orbit after a supernova: instantaneous mass loss and
        natal kick (see kicks/orbit.py)."""
        if dying_star.state_history[-1] not in STAR_STATES_CC:
            raise ValueError("There is no information in the evolutionary "
                             "history about STAR_STATES_CC.")
        # the kick is drawn, and recorded on the star, even if there is no
        # orbit to update
        v_kick, cos_theta, phi, mean_anomaly = self._natal_kick(dying_star)
        if binary.state in NO_ORBIT_STATES:
            self._clear_orbit(binary)
        elif "massless_remnant" in (binary.star_1.state, binary.star_2.state):
            self._record_tilts(binary,dying_star,other_star)
            self._set_disrupted(binary)
        else:
            (bound, A_post, e_post, V_sys, tilt,
             true_anomaly) = orbit_after_kick(
                M_pre=dying_star.mass_history[-1], M_post=dying_star.mass,
                M_companion=other_star.mass, separation=binary.separation,
                eccentricity=binary.eccentricity, v_kick=v_kick,
                cos_theta=cos_theta, phi=phi, mean_anomaly=mean_anomaly,
                verbose=self.verbose)
            reset_binary_properties(binary,keep=("nearest_neighbour_distance",
                                                 "event"))
            if bound:
                self._record_tilts(binary,dying_star,other_star,tilt,
                                   true_anomaly)
                binary.state = "detached"
                binary.event = None
                binary.separation = A_post
                binary.eccentricity = e_post
                binary.V_sys = V_sys
                binary.time = binary.time_history[-1]
                binary.orbital_period = orbital_period_from_separation(
                    binary.separation, binary.star_1.mass, binary.star_2.mass)
                binary.mass_transfer_case = 'None'
            else:
                self._record_tilts(binary,dying_star,other_star)
                self._set_disrupted(binary)

    def _update_orbit_after_WD(self,binary,dying_star,other_star):
        """A WD gets no kick: the orbit only widens due to the gradually lost
        mass."""
        separation = separation_evol_wind_loss(
            dying_star.mass, dying_star.mass_history[-1], other_star.mass,
            binary.separation)
        orbital_period = orbital_period_from_separation(
            separation, binary.star_1.mass, binary.star_2.mass)
        reset_binary_properties(binary,keep=("V_sys",
                                             "nearest_neighbour_distance",
                                             "state"))
        binary.separation = separation
        if binary.state not in NO_ORBIT_STATES:
            binary.state = "detached"
        binary.event = None
        binary.time = binary.time_history[-1]
        binary.eccentricity = binary.eccentricity_history[-1]
        binary.orbital_period = orbital_period
        binary.mass_transfer_case = 'None'
        dying_star.natal_kick_velocity = 0.0

    @staticmethod
    def _record_tilts(binary,dying_star,other_star,tilt=np.nan,
                      true_anomaly=None):
        """Store the spin-orbit tilts; nan if the binary was disrupted."""
        if not binary.first_SN_already_occurred:
            binary.star_1.spin_orbit_tilt_first_SN = tilt
            binary.star_2.spin_orbit_tilt_first_SN = tilt
            binary.first_SN_already_occurred = True
            if true_anomaly is not None:
                binary.true_anomaly_first_SN = true_anomaly
        elif true_anomaly is None:
            binary.star_1.spin_orbit_tilt_second_SN = np.nan
            binary.star_2.spin_orbit_tilt_second_SN = np.nan
        else:
            # assume the dying star's spin had aligned with the pre-SN orbit
            dying_star.spin_orbit_tilt_second_SN = tilt
            other_star.spin_orbit_tilt_second_SN = combined_tilt(
                tilt_1=other_star.spin_orbit_tilt_first_SN, tilt_2=tilt,
                true_anomaly_1=binary.true_anomaly_first_SN,
                true_anomaly_2=true_anomaly)
            binary.true_anomaly_second_SN = true_anomaly

    @staticmethod
    def _clear_orbit(binary):
        """The system had no orbit before the SN (disrupted, merged or
        single)."""
        reset_binary_properties(binary,keep=("nearest_neighbour_distance",
                                             "state"))
        binary.event = None
        binary.separation = np.nan
        binary.eccentricity = np.nan
        binary.V_sys = np.array([0, 0, 0])
        binary.time = binary.time_history[-1]
        binary.orbital_period = np.nan
        binary.mass_transfer_case = 'None'
        binary.first_SN_already_occurred = True

    @staticmethod
    def _set_disrupted(binary):
        """The SN disrupted the binary."""
        binary.state = "disrupted"
        binary.event = None
        binary.separation = np.nan
        binary.eccentricity = np.nan
        binary.V_sys = np.array([0, 0, 0])
        binary.time = binary.time_history[-1]
        binary.orbital_period = np.nan
        binary.mass_transfer_case = 'None'

    def finalise_step(self,binary):
        """Inspiral time of a double compact object, or the next CC event."""
        state1, state2 = binary.star_1.state, binary.star_2.state
        if (binary.state == "disrupted"
                or "massless_remnant" in (state1, state2)):
            binary.inspiral_time = np.nan
        elif state1 in STAR_STATES_CO and state2 in STAR_STATES_CO:
            binary.inspiral_time = inspiral_timescale_from_separation(
                binary.star_1.mass, binary.star_2.mass, binary.separation,
                binary.eccentricity)
        # the core collapse of the companion may follow immediately
        elif state1 in STAR_STATES_CO and state2 in STAR_STATES_C_DEPLETED:
            binary.event = "CC2"
        elif state1 in STAR_STATES_C_DEPLETED and state2 in STAR_STATES_CO:
            binary.event = "CC1"
        if self.verbose:
            print("End of step SN:\n", binary)
