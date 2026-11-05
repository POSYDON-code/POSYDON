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
    "Dimitris Souropanis <dsouropanis@ia.forth.gr>"
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
)
from posydon.binary_evol.SN.explodability.ECSN_prescriptions import (
    get_ECSN_prescription,
)
from posydon.binary_evol.SN.explodability.PISN_prescriptions import (
    get_PISN_prescription,
)
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
        unknown_kwargs = set(kwargs) - set(self.DEFAULT_KWARGS)
        if unknown_kwargs:
            raise ValueError(f"Invalid arguments passed to step_NS: {unknown_kwargs}")
        config = self.DEFAULT_KWARGS.copy()
        config.update(kwargs)
        config = replace_deprecated_mechanism(config)
        config = replace_deprecated_kick_normalisation(config)
        if config["RNG"] is None:
            config["RNG"] = np.random.default_rng()
        if config["max_neutrino_mass_loss"] is None:
            config["max_neutrino_mass_loss"] = 0
        for key, value in config.items():
            setattr(self, key, value)

        # also tells whether a star is above the pair-instability gap
        self._PISN_check = get_PISN_prescription(**config)
        self._SN_type_checks = (
            get_ECSN_prescription(**config),     # Checks whether a model ends as a WD or ECSN
            self._PISN_check,                    # Checks whether a model ends as a PISN or PPISN
            lambda star: "CCSN",                 # Everything else gets treated as a CCSN directly
        )
        self._explosion_criterion = get_CCSN_prescription(**config)
        self._get_remnant_mass = get_remnant_mass_prescription(**config)
        self._update_params = Remnant_Finaliser(**config)
        self._natal_kick = get_kick_prescription(**config)

    def __call__(self,binary):
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
        if binary.event=="CC1":
            return binary.star_1, binary.star_2
        elif binary.event=="CC2":
            return binary.star_2, binary.star_1
        else:
            raise ValueError(
                "Invalid call of step_SN: \n" + 
                "Called during binary event {event}".format(event=binary.event) + 
                "Step_SN should only be called if binary event is CC1 or CC2"
            )

    def collapse_star(self,star):
        if star.state not in STAR_STATES_CC:
            return f"The star cannot collapse: star state {star.state}."
        SN_type = self.collapse_outcome(star)
        explosion = self.explosion_outcome(star, SN_type)
        SN_type, m_rembar, f_fb, state = self._get_remnant_mass(star, SN_type, explosion)
        # returns an error message if the collapse fails, else None
        return self._update_params(star, SN_type, m_rembar, f_fb, state)

    def explosion_outcome(self,star,SN_type):
        """Return the outcome of the explosion criterion, for CCSN, PPISN and PISN.

        A PISN leaves no remnant, but its f_fb comes from the CCSN
        prescription (as on main), which needs the explosion outcome.
        Above the pair-instability gap the core collapses directly into a BH
        ("BH_direct"), whatever the mechanism, as in Marchant et al. (2019)
        and Hendriks et al. (2023).
        """
        if SN_type not in ("CCSN", "PPISN", "PISN"):
            return None
        if self._PISN_check.is_above_gap(star):
            return "BH_direct"
        return self._explosion_criterion(star)

    def collapse_outcome(self,star):
        return next(outcome for outcome in (check(star) for check in self._SN_type_checks)
                    if outcome is not None)

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
        if dying_star.state_history[-1] not in STAR_STATES_CC:
            raise ValueError("There is no information in the evolutionary "
                             "history about STAR_STATES_CC.")
        # the kick is drawn even if there is no orbit to update, as on main
        v_kick, cos_theta, phi, mean_anomaly = self._natal_kick(dying_star)
        if binary.state in NO_ORBIT_STATES:
            self._clear_orbit(binary)
        elif "massless_remnant" in (binary.star_1.state, binary.star_2.state):
            self._record_tilts(binary,dying_star,other_star)
            self._set_disrupted(binary)
        else:
            bound, A_post, e_post, V_sys, tilt, true_anomaly = orbit_after_kick(
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
        """A WD gets no kick: the orbit only widens due to the gradually lost mass."""
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
        """The system had no orbit before the SN (disrupted, merged, single)."""
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
