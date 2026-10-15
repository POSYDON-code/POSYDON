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
from posydon.grids.SN_MODELS import DEFAULT_SN_MODEL
from posydon.utils.common_functions import (
    set_binary_to_failed,
)
from posydon.utils.posydonerror import ModelError

####################################################################################################
from abc import ABC, abstractmethod
from posydon.binary_evol.SN.explodability.ECSN_prescriptions import get_ECSN_prescription
from posydon.binary_evol.SN.explodability.PISN_prescriptions import get_PISN_prescription
from posydon.binary_evol.SN.explodability.CCSN_explosion_criteria import get_CCSN_prescription
from posydon.binary_evol.SN.remnant_mass.remnant_mass_prescriptions import get_remnant_mass_prescription
from posydon.binary_evol.SN.kicks.kick_prescriptions import get_kick_prescription
####################################################################################################

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
        if config["RNG"] is None:
            config["RNG"] = np.random.default_rng()
        for key, value in config.items():
            setattr(self, key, value)

        self._explosion_checksget_ECSN_prescription = (
            get_ECSN_prescription(**config),     # -> "WD" | "ECSN" | None
            get_CCSN_prescription(**config),     # -> "CCSN"
        )
        self._get_remnant_mass = get_remnant_mass_prescription(**kwargs)
        self._get_kick = get_kick_prescription(**kwargs)

    def __call__(self,binary):
        dying_star,other_star = self._select_dying_star(binary)
        model_err = self.collapse_star(dying_star)
        if model_err is not None:
            set_binary_to_failed(binary)
            raise ModelError(model_err)
        self._reset_other_star_properties(star=other_star)
        binary.update_star_states()
        self.orbital_kick(binary)
        self.finalise_step(binary)

    @staticmethod
    def _select_dying_star(binary):
        if binary.event=="CC1":
            return binary.star_1, binary.star_2
        elif binary.event=="CC2":
            return binary.star_2, binary.star_1
        else:
            raise ValueError(
                "Invalid call of step_SN:"
                "Called during binary event {event}".format(event=binary.event)
                "Step_SN should only be called if binary event is CC1 or CC2"
            )

    def collapse_star(self,star):
        outcome = self.collapse_outcome(star)
        m_rem, f_fb = self._remnant_mass(star, outcome)
        self._update_params(star, outcome, m_rem, f_fb)

    def collapse_outcome(self,star):
        return next(outcome for outcome in (check(star) for check in self._explosion_checks)
                    if outcome is not None)

    def _reset_other_star_properties(self, star):
        """Reset the properties of the star that is not being collapsed."""
        star.lg_mdot = None
        star.lg_system_mdot = None

    def finalise_step(self,binary):
        return