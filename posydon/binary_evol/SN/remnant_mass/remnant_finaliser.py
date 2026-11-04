"""Turn the remnant mass of a collapsing star into the final remnant.

The remnant-mass stage gives (SN_type, m_rembar, f_fb, state); this stage
updates the star accordingly: a WD, a massless remnant (PISN), or a NS/BH
with its gravitational mass, spin, disk and ejecta. Ported from
StepSN.collapse_star on main.
"""

import numpy as np

import posydon.utils.constants as const
from posydon.binary_evol.singlestar import (
    STARPROPERTIES,
    convert_star_to_massless_remnant,
)
from posydon.binary_evol.SN.profile_collapse import (
    do_core_collapse_BH,
    get_ejecta_element_mass_at_collapse,
)
from posydon.utils.common_functions import CO_radius
from posydon.utils.limits_thresholds import THRESHOLD_CENTRAL_ABUNDANCE
from posydon.utils.posydonwarning import Pwarn

# Properties kept by a new compact object; all other STARPROPERTIES are reset.
KEEP_COMPACT_OBJECT = ("state", "mass", "spin", "log_R", "metallicity",
                       "m_disk_accreted", "m_disk_radiated", "co_core_mass")
# Properties kept by a new WD. log_R is kept, unlike on main (bug 12).
KEEP_WD = ("state", "mass", "spin", "log_R", "he_core_mass", "co_core_mass",
           "m_disk_accreted", "m_disk_radiated", "center_h1", "center_he4",
           "center_c12", "center_n14", "center_o16")


def reset_star_properties(star,keep):
    """Set every STARPROPERTIES key that is not in keep to None."""
    for key in STARPROPERTIES:
        if key not in keep:
            setattr(star,key,None)


class Remnant_Finaliser:
    """Turn the remnant mass into a WD, NS, BH or massless remnant."""
    def __init__(self,**kwargs):
        self.max_NS_mass = kwargs['max_NS_mass']
        self.max_neutrino_mass_loss = kwargs['max_neutrino_mass_loss']
        self.use_profiles = kwargs['use_profiles']
        self.use_core_masses = kwargs['use_core_masses']
        self.allow_spin_None = kwargs['allow_spin_None']
        self.verbose = kwargs['verbose']
        # Ways to set the spin and disk of a NS/BH, in order of preference.
        # Each returns False if it does not apply to the star.
        self._spin_methods = (
            self._spin_from_profile,
            self._spin_from_core_masses,
            self._spin_None,
        )

    def __call__(self,star,SN_type,m_rembar,f_fb,state):
        """Update the star; return an error message if the collapse fails."""
        for key in ['m_disk_accreted','m_disk_radiated']:
            if getattr(star,key) is None:
                setattr(star,key,np.nan)
        star.SN_type = SN_type
        star.f_fb = f_fb
        if state=="WD":
            return self._finalise_WD(star,m_rembar)
        elif state=="massless_remnant":
            convert_star_to_massless_remnant(star=star)
            return None
        return self._finalise_compact_object(star,SN_type,m_rembar,state)

    def _finalise_WD(self,star,m_rembar):
        star.mass = m_rembar
        star.state = "WD"
        star.spin = 0.0
        star.log_R = np.log10(CO_radius(star.mass,star.state))
        star.he_core_mass = star.mass
        if star.center_he4 < THRESHOLD_CENTRAL_ABUNDANCE:
            star.co_core_mass = star.mass
        else:
            star.co_core_mass = 0.0
        reset_star_properties(star,keep=KEEP_WD)
        return None

    def _finalise_compact_object(self,star,SN_type,m_rembar,state):
        m_grav = self._gravitational_mass(m_rembar)
        state = self._NS_or_BH(state,m_grav)
        for set_spin in self._spin_methods:
            if set_spin(star,SN_type,m_rembar,m_grav,state):
                break
        else:
            reset_star_properties(star,keep=())
            return "FAILED core collapse!"
        star.state = state
        star.metallicity = star.metallicity_history[-1]
        star.log_R = np.log10(CO_radius(star.mass,star.state))
        reset_star_properties(star,keep=KEEP_COMPACT_OBJECT)
        return None

    def _gravitational_mass(self,m_rembar):
        """Lattimer & Yahil (1989), with the neutrino mass loss capped."""
        m_grav = (20.0 / 3.0) * (np.sqrt(1.0 + 0.3 * m_rembar) - 1.0)
        return max(m_grav, m_rembar - self.max_neutrino_mass_loss)

    def _NS_or_BH(self,state,m_grav):
        """The mass decides: a remnant at or above max_NS_mass is a BH.

        If the explosion criterion gave a different state, the mass wins,
        with a warning. On main, an engine BH below max_NS_mass fails the
        binary (profile path) or becomes a NS (core-mass path); this is
        pending a check on the grids.
        """
        mass_state = "BH" if m_grav >= self.max_NS_mass else "NS"
        if state is not None and state != mass_state:
            Pwarn(f"The explosion criterion gave a {state}, but "
                  f"m_grav={m_grav:.3f} Msun (max_NS_mass={self.max_NS_mass}) "
                  f"makes it a {mass_state}.", "ReplaceValueWarning")
        return mass_state

    def _spin_from_profile(self,star,SN_type,m_rembar,m_grav,state):
        """Collapse the stellar profile to get the BH spin and disk."""
        if not (self.use_profiles and star.profile is not None):
            return False
        if state=="BH":
            # m_rembar - m_grav <= max_neutrino_mass_loss by construction
            delta_M = m_rembar - m_grav
            final_BH = do_core_collapse_BH(
                star=star, mass_collapsing=m_rembar,
                mass_central_BH=self.max_NS_mass + delta_M,
                neutrino_mass_loss=delta_M,
                max_neutrino_mass_loss=self.max_neutrino_mass_loss,
                verbose=self.verbose
            )
            for key, value in final_BH.items():
                setattr(star,key,value)
            star.mass = final_BH['M_BH_total']
            star.spin = final_BH['a_BH_total']
            star.m_disk_accreted = final_BH['m_disk_accreted']
            star.m_disk_radiated = final_BH['m_disk_radiated']
        else:
            star.mass = m_grav
            star.spin = 0.0
            star.m_disk_accreted = 0.0
            star.m_disk_radiated = 0.0
        star.h1_mass_ej, star.he4_mass_ej = get_ejecta_element_mass_at_collapse(
            star,star.mass,verbose=self.verbose)
        return True

    def _spin_from_core_masses(self,star,SN_type,m_rembar,m_grav,state):
        """Without a profile: the BH spin comes from the star's angular momentum."""
        # ECSNe always use this method (on main this is bug 10)
        if not (self.use_core_masses or SN_type=="ECSN"):
            return False
        star.mass = m_grav
        if state=="BH":
            if SN_type=="ECSN":
                Pwarn(f"An ECSN should not form a black hole: m_grav={m_grav}.",
                      "InappropriateValueWarning")
            # Eq. 14 of Fryer et al. (2012), in CGS units
            spin = (10**star.log_total_angular_momentum * const.clight
                    / (const.standard_cgrav * (m_grav * const.Msun) ** 2))
            star.spin = min(spin, 1.0)
        else:
            star.spin = 0.0
        star.m_disk_accreted = 0.0
        star.m_disk_radiated = 0.0
        star.h1_mass_ej, star.he4_mass_ej = np.nan, np.nan
        return True

    def _spin_None(self,star,SN_type,m_rembar,m_grav,state):
        """Leave the spin undetermined (for interpolator training)."""
        if not self.allow_spin_None:
            return False
        star.mass = m_grav
        star.spin = None
        star.m_disk_accreted = 0.0
        star.m_disk_radiated = 0.0
        star.h1_mass_ej, star.he4_mass_ej = np.nan, np.nan
        return True
