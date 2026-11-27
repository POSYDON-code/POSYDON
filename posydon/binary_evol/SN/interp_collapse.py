"""Find the outcome of core collapse using values interpolated from the post-processed grids.

Before a star collapses, its interpolated outcome is stored on it for every
SN model of SN_MODELS (star.SN_MODEL_v2_XX): state, SN_type, f_fb, mass,
spin, disk masses, M4, mu4 and ejecta. With use_interp_values=True, StepSN
uses that outcome when it is consistent with the star, and computes the
collapse otherwise.
"""

import numpy as np
import pandas as pd

from posydon.binary_evol.singlestar import convert_star_to_massless_remnant
from posydon.binary_evol.SN.remnant_mass.remnant_finaliser import (
    KEEP_COMPACT_OBJECT,
    KEEP_WD,
    reset_star_properties,
)
from posydon.grids.SN_MODELS import get_SN_MODEL_NAME
from posydon.utils.common_functions import CO_radius
from posydon.utils.limits_thresholds import THRESHOLD_CENTRAL_ABUNDANCE
from posydon.utils.posydonwarning import Pwarn

# Interpolated SN types accepted for each SN type of the star. PPISN and PISN
# candidates end as a CCSN, PPISN or PISN depending on their remnant mass, so
# all three are accepted for them.
ACCEPTED_SN_TYPES = {
    "WD": ("WD",),
    "CCSN": ("CCSN",),
    "PPISN": ("CCSN", "PPISN", "PISN"),
    "PISN": ("CCSN", "PPISN", "PISN"),
}


def check_SN_CO_match(SN_type, state):
    """True if the SN type is consistent with the state of the remnant."""
    if state == "PISN":
        state = "massless_remnant"
    if state == "WD":
        return SN_type == "WD"
    if state == "NS":
        return SN_type in ("ECSN", "CCSN")
    if state == "BH":
        return SN_type in ("CCSN", "PPISN")
    if state == "massless_remnant":
        return SN_type == "PISN"
    return True


class InterpCollapse:
    """Apply the interpolated outcome of the collapse to a star."""
    def __init__(self,**kwargs):
        # The ECSN and use_* options are not compared.
        self.SN_MODEL_NAME = get_SN_MODEL_NAME(kwargs,
                                               verbose=kwargs['verbose'])
        if self.SN_MODEL_NAME is None:
            raise ValueError(
                "use_interp_values=True, but no SN model in "
                "posydon.grids.SN_MODELS matches the settings of the SN step "
                "(the ECSN and use_* options are not compared). Change the "
                "settings or set use_interp_values=False.")
        self.use_profiles = kwargs['use_profiles']
        self.use_core_masses = kwargs['use_core_masses']
        self.allow_spin_None = kwargs['allow_spin_None']
        self.verbose = kwargs['verbose']
        if self.verbose:
            print(f"InterpCollapse: interpolated values of "
                  f"{self.SN_MODEL_NAME}")

    def __call__(self,star,SN_type):
        """Apply the interpolated outcome to the star, if it can be used.

        SN_type is the star's own SN type, from the classification of the
        step. Returns True if the outcome was applied, and False if the
        collapse has to be computed: with a warning saying why, except for
        ECSNe, which are always computed.
        """
        outcome = getattr(star, self.SN_MODEL_NAME)
        if outcome is None:
            Pwarn(f"{self.SN_MODEL_NAME}: The collapsed star was not "
                  "interpolated! If use_profiles or use_core_masses is set "
                  "to True, continue with the collapse.",
                  "InterpolationWarning")
            return False
        if SN_type == "ECSN":
            # decided by the ECSN option of the step, not the interpolators
            return False
        problem = self._problem(outcome, SN_type)
        if problem is not None:
            Pwarn(problem + self._alternative(star), "ApproximationWarning")
            return False
        self._apply(star, outcome)
        if self.verbose:
            print(f"InterpCollapse: used the interpolated outcome of "
                  f"{self.SN_MODEL_NAME}: {outcome['SN_type']}, "
                  f"{outcome['state']}, mass {outcome['mass']}")
        return True

    def _problem(self,outcome,SN_type):
        """Why the interpolated outcome cannot be used, or None."""
        if outcome['SN_type'] not in ACCEPTED_SN_TYPES[SN_type]:
            return (f"The SN_type does not match the star: {SN_type}"
                    f"!={outcome['SN_type']}.")
        if not check_SN_CO_match(outcome['SN_type'], outcome['state']):
            return (f"{self.SN_MODEL_NAME}: The SN_type does not match the "
                    "predicted CO.")
        if pd.isna(outcome['mass']):
            return "There is no interpolated remnant mass."
        return None

    def _alternative(self,star):
        """How the collapse continues, for the warnings."""
        if self.use_profiles and star.profile is not None:
            return " Instead use profiles."
        if self.use_core_masses:
            return " Instead use core masses."
        if self.allow_spin_None:
            return " Instead use core mass without spin."
        return ""

    @staticmethod
    def _apply(star,outcome):
        """Copy the interpolated outcome to the star (outcome is unchanged)."""
        for key, value in outcome.items():
            setattr(star, key, value)
        if star.state == "WD":
            star.he_core_mass = star.mass
            if star.center_he4 < THRESHOLD_CENTRAL_ABUNDANCE:
                star.co_core_mass = star.mass
            else:
                star.co_core_mass = 0.0
            reset_star_properties(star, keep=KEEP_WD)
        else:
            star.metallicity = star.metallicity_history[-1]
            reset_star_properties(star, keep=KEEP_COMPACT_OBJECT)
        if star.SN_type == "PISN":
            convert_star_to_massless_remnant(star=star)
        else:
            star.log_R = np.log10(CO_radius(star.mass, star.state))
