import copy
import numbers
from abc import ABC, abstractmethod

import numpy as np

import posydon.utils.constants as const
from posydon.binary_evol.SN.explodability.CCSN_explosion_criteria import (
    SN_MECHANISMS,
)
from posydon.utils.posydonerror import ModelError
from posydon.utils.posydonwarning import Pwarn


def get_remnant_mass_prescription(**kwargs):
    option = kwargs['mechanism']
    if option in REMNANT_MASS_PRESCRIPTIONS:
        return REMNANT_MASS_PRESCRIPTIONS[option](**kwargs)
    if option in SN_MECHANISMS:
        raise NotImplementedError(
            f"The remnant mass for mechanism {option!r} has not been ported "
            "to the new step_SN structure yet."
        )
    raise ValueError(
        f"Invalid option {option!r} given to remnant mass prescription. "
        f"Valid options are: {', '.join(SN_MECHANISMS)}."
    )


STATE_FROM_EXPLOSION = {
    "NS": "NS",
    "BH_fallback": "BH",
    "BH_direct": "BH",
    None: None,
}

# NS masses from ECSNe are computed separately from those in the CCSN 
# mechanism. The whole remnant is treated as a PNS, and there is no
# fallback (following Giacobbo & Mapelli 2020).
def ECSN_mass_Tauris(star):
    return star.co_core_mass, 0.0, "NS"

def ECSN_mass_Podsiadlowski(star):
    return 1.38, 0.0, "NS"

def ECSN_mass_No_ECSN(star):
    raise ValueError("A star was classified as ECSN with the ECSN option "
                     "'No_ECSN'.")

ECSN_MASS_OPTIONS = {
    "Tauris+15":ECSN_mass_Tauris,
    "Podsiadlowski+04":ECSN_mass_Podsiadlowski,
    "No_ECSN":ECSN_mass_No_ECSN,
    None:ECSN_mass_No_ECSN,
}

class PPISN_Mass_Base(ABC):
    """Mass left after pulsational pair instability (PPI)."""
    def __init__(self,remnant_mass,**kwargs):
        # the remnant-mass prescription, for PPI prescriptions that need it
        self.remnant_mass = remnant_mass

    @abstractmethod
    def __call__(self,star,explosion):
        """Return the mass left after the pulses, nan if the star is disrupted."""
        raise NotImplementedError

def PISN_mass_Marchant(star):
    """
    Mass left after pulsational pair instability, from Marchant+19.

    8th-order polynomial fit of the BH masses in table 1 of Marchant et al.
    (2019), see the COSMIC paper (Breivik et al. 2020). Valid for the PPISN
    range of PISN_check_Marchant19, 31.99 <= m_He_core <= 61.10 Msun. The
    ranges use He-core masses from table 1 of Marchant et al. (2019) instead
    of the CO-core range 54.48 < M_CO/Msun < 113.29 quoted by Breivik et al.
    (2020), which might cause gaps when switching between core masses.
    """
    m_He_core = star.he_core_mass
    return (
        - 6.29429263e5
        + 1.15957797e5 * m_He_core
        - 9.28332577e3 * m_He_core ** 2.0
        + 4.21856189e2 * m_He_core ** 3.0
        - 1.19019565e1 * m_He_core ** 4.0
        + 2.13499267e-1 * m_He_core ** 5.0
        - 2.37814255e-3 * m_He_core ** 6.0
        + 1.50408118e-5 * m_He_core ** 7.0
        - 4.13587235e-8 * m_He_core ** 8.0
    )

class PPISN_Mass_Marchant19(PPISN_Mass_Base):
    def __call__(self,star,explosion):
        return PISN_mass_Marchant(star)


class PPISN_Mass_Hendriks23(PPISN_Mass_Base):
    """
    Mass left after pulsational pair instability from 
    Hendriks et al. (2023), MNRAS, 526, 4130.

    The pulses remove delta_M_PPI (Eq. 6) from the star, then the
    remnant-mass prescription is applied to what is left. A remnant below
    10 Msun counts as a PISN.
    """
    def __init__(self,remnant_mass,**kwargs):
        super().__init__(remnant_mass,**kwargs)
        shift = kwargs['PISN_CO_shift']
        extra = kwargs['PPI_extra_mass_loss']
        self.delta_M_CO_shift = shift if shift is not None else 0.0
        self.delta_M_PPI_extra_ML = extra if extra is not None else 0.0
        self.conserve_hydrogen_PPI = kwargs['conserve_hydrogen_PPI']

    def __call__(self,star,explosion):
        delta_M_PPI = self.mass_lost_in_pulses(star)
        if delta_M_PPI <= 0.0:
            # no pulses: nothing limits the CCSN remnant mass
            return self.remnant_mass._m_collapsing(star)
        if self.conserve_hydrogen_PPI:
            m_after_pulses = star.mass - delta_M_PPI
        else:
            m_after_pulses = star.he_core_mass - delta_M_PPI
        if m_after_pulses < 0.0:
            return np.nan
        stripped_star = copy.deepcopy(star)
        stripped_star.mass = m_after_pulses
        stripped_star.he_core_mass = min(star.he_core_mass, m_after_pulses)
        stripped_star.co_core_mass = min(star.co_core_mass, m_after_pulses)
        # NOTE: uses the explosion outcome of the original star; on main an
        # engine re-decides for the stripped star (to revisit with engines)
        m_rembar, _, _ = self.remnant_mass._get_rem_mass_core_collapse(
            stripped_star, explosion)
        if m_rembar < 10.0:
            return np.nan
        return m_rembar

    def mass_lost_in_pulses(self,star):
        """delta_M_PPI from Eq. 6 of Hendriks et al. (2023)."""
        # delta_M_PPI -> -inf for Z -> 0, so Z is floored at 1e-4, the
        # lowest (absolute) metallicity in Hendriks et al. (2023)
        Z = max(star.metallicity * const.Zsun, 1e-4)
        x = star.co_core_mass - self.delta_M_CO_shift - 34.8
        return ((0.0006 * np.log10(Z) + 0.0054) * x**3
                - 0.0013 * x**2 + self.delta_M_PPI_extra_ML)

class PPISN_Mass_constant(PPISN_Mass_Base):
    def __init__(self,remnant_mass,**kwargs):
        super().__init__(remnant_mass,**kwargs)
        self.PISN_mass_cap = float(kwargs['PISN'])
    def __call__(self,star,explosion):
        return self.PISN_mass_cap

PPISN_MASS_PRESCRIPTIONS = {
    "Marchant+19": PPISN_Mass_Marchant19,
    "Hendriks+23": PPISN_Mass_Hendriks23,
}

def get_PISN_mass_prescription(remnant_mass,**kwargs):
    option=kwargs['PISN']
    if isinstance(option,numbers.Real) and not isinstance(option,bool):
        return PPISN_Mass_constant(remnant_mass,**kwargs)
    if option in PPISN_MASS_PRESCRIPTIONS:
        return PPISN_MASS_PRESCRIPTIONS[option](remnant_mass,**kwargs)
    return None

class Remnant_Mass_Base(ABC):
    def __init__(self,**kwargs):
        self.mechanism = kwargs['mechanism']
        self.ECSN_option = kwargs['ECSN']
        self.PISN_option = kwargs['PISN']
        self.conserve_hydrogen_envelope = kwargs['conserve_hydrogen_envelope']
        self.verbose = kwargs['verbose']
        # last, because PPI prescriptions may use the settings above
        self._mass_after_pulses = get_PISN_mass_prescription(self,**kwargs)

    def __call__(self,star,SN_type,explosion=None):
        """Return (SN_type, m_rembar, f_fb, state) of the remnant.

        SN_type can differ from the input: a PPISN candidate whose pulses
        do not reduce the remnant mass is a CCSN, as on main.
        state is "WD", "NS", "BH", "massless_remnant", or None if it has
        to be decided from the remnant mass.
        """
        if SN_type=="WD":
            return (SN_type, *WD_mass(star))
        elif SN_type=="ECSN":
            return (SN_type, *self._get_rem_mass_ECSN(star))
        elif SN_type=="CCSN":
            return (SN_type, *self._get_rem_mass_core_collapse(star,explosion))
        elif SN_type=="PPISN":
            return self._get_rem_mass_PPISN(star,explosion)
        elif SN_type=="PISN":
            return (SN_type, *self._get_rem_mass_PISN(star,explosion))
        else:
            raise ValueError(f"Invalid SN type {SN_type!r} passed during remnant mass calculation")

    def _get_rem_mass_PPISN(self,star,explosion):
        """
        Pulsational pair instability pulses remove mass before collapse.

        The CCSN remnant mass is capped at the mass left after the pulses.
        If the cap does not reduce it, the star is labelled CCSN, and if the
        pulses disrupt the star it is labelled PISN.
        """
        m_rembar, f_fb, state = self._get_rem_mass_core_collapse(star,explosion)
        m_after_pulses = self._mass_after_pulses(star,explosion)
        if np.isnan(m_after_pulses):
            # the pulses disrupt the star, label changes to a PISN
            return "PISN", np.nan, f_fb, "massless_remnant"
        if m_rembar > m_after_pulses:
            return "PPISN", m_after_pulses, f_fb, state
        return "CCSN", m_rembar, f_fb, state

    def _get_rem_mass_PISN(self,star,explosion):
        """Pair-instability SNe disrupt the star, leaving no remnant behind"""
        # f_fb is kept from the CCSN prescription
        _, f_fb, _ = self._get_rem_mass_core_collapse(star,explosion)
        return np.nan, f_fb, "massless_remnant"

    def _m_collapsing(self,star):
        """Choose the mass that collapses. Either the whole star, or only its He core"""
        if self.conserve_hydrogen_envelope:
            return star.mass
        return star.he_core_mass

    def _get_rem_mass_direct_collapse(self,star):
        """Failed explosion: the whole collapsing mass forms a BH."""
        return self._m_collapsing(star), 1.0, "BH"

    def _get_rem_mass_core_collapse(self,star,explosion):
        """
        Return (m_rembar, f_fb, state) after the core collapse of the star.

        A failed explosion ("BH_direct", e.g. above the pair-instability gap)
        is a direct collapse for every mechanism; otherwise the mechanism's
        CCSN prescription applies.
        """
        if explosion == "BH_direct":
            return self._get_rem_mass_direct_collapse(star)
        return self._get_rem_mass_CCSN(star,explosion)

    def _get_rem_mass_ECSN(self,star):
        return ECSN_MASS_OPTIONS[self.ECSN_option](star)

    @abstractmethod
    def _get_rem_mass_CCSN(self,star,explosion):
        """Return (m_rembar, f_fb, state) for a star classified as CCSN."""
        raise NotImplementedError

    def __repr__(self):
        return f"{type(self).__name__}(mechanism={self.mechanism!r})"


class Remnant_Mass_constant(Remnant_Mass_Base):
    def __init__(self,**kwargs):
        super().__init__(**kwargs)


class Remnant_Mass_Fryer12(Remnant_Mass_Base):
    """Remnant mass from Fryer et al. (2012), ApJ, 749, 91.

    m_rembar = m_proto + m_fb, where m_fb is the part of the mass outside the
    proto-compact object (m_collapsing - m_proto) that falls back.
    Subclasses give m_proto and m_fb as functions of the CO core mass.
    """
    def _get_rem_mass_CCSN(self,star,explosion):
        m_CO = star.co_core_mass
        m_proto = self._proto_mass(m_CO)
        m_outside_proto = self._m_collapsing(star) - m_proto
        m_fb = self._fallback_mass(m_CO,m_outside_proto)
        f_fb = m_fb / m_outside_proto
        return m_proto + m_fb, f_fb, STATE_FROM_EXPLOSION[explosion]

    @abstractmethod
    def _proto_mass(self,m_CO):
        raise NotImplementedError

    @abstractmethod
    def _fallback_mass(self,m_CO,m_outside_proto):
        raise NotImplementedError


class Remnant_Mass_Fryer_delayed(Remnant_Mass_Fryer12):
    """Delayed mechanism, Eqs. 17-20 of Fryer et al. (2012)."""
    def _proto_mass(self,m_CO):
        if m_CO < 3.5:
            return 1.2
        elif m_CO < 6.0:
            return 1.3
        elif m_CO < 11.0:
            return 1.4
        else:
            return 1.6

    def _fallback_mass(self,m_CO,m_outside_proto):
        if m_CO < 2.5:
            return 0.2
        elif m_CO < 3.5:
            return 0.5 * m_CO - 1.05
        elif m_CO < 11.0:
            a = 0.133 - 0.093 / m_outside_proto
            b = -11.0 * a + 1.0
            return (a * m_CO + b) * m_outside_proto
        else:
            return m_outside_proto


class Remnant_Mass_Fryer_rapid(Remnant_Mass_Fryer12):
    """Rapid mechanism, Eqs. 15-17 of Fryer et al. (2012).

    The proto-compact object mass of 1.1 Msun is from Giacobbo & Mapelli
    (2020), ApJ, 891, 141.
    """
    def _proto_mass(self,m_CO):
        return 1.1

    def _fallback_mass(self,m_CO,m_outside_proto):
        if m_CO < 2.5:
            return 0.2
        elif m_CO < 6.0:
            return 0.286 * m_CO - 0.514
        elif m_CO < 7.0:
            return m_outside_proto
        elif m_CO < 11.0:
            a = 0.25 - 1.275 / m_outside_proto
            b = -11.0 * a + 1.0
            return (a * m_CO + b) * m_outside_proto
        else:
            return m_outside_proto



class Remnant_Mass_Sukhbold16(Remnant_Mass_Base):
    def __init__(self,**kwargs):
        super().__init__(**kwargs)

class Remnant_Mass_Patton20(Remnant_Mass_Base):
    def __init__(self,**kwargs):
        super().__init__(**kwargs)

class Remnant_Mass_direct(Remnant_Mass_Base):
    def _get_rem_mass_CCSN(self, star, explosion):
        return self._m_collapsing(star), 1.0, STATE_FROM_EXPLOSION[explosion]
    def _get_rem_mass_ECSN(self, star):
        return self._m_collapsing(star), 1.0, "NS"

# class Remnant_Mass_Mueller16(Remnant_Mass_Base):
#     def __init__(**kwargs):
#         super().__init__(**kwargs)

# class Remnant_Mass_Boccioli24(Remnant_Mass_Base):
#     def __init__(**kwargs):
#         super().__init__(**kwargs)


REMNANT_MASS_PRESCRIPTIONS = {
    "Fryer+12-rapid": Remnant_Mass_Fryer_rapid,
    "Fryer+12-delayed": Remnant_Mass_Fryer_delayed,
    "direct": Remnant_Mass_direct,
}


def WD_mass(star):
    m_co_core = star.co_core_mass
    m_He_core = star.he_core_mass
    m_star = star.mass
    if m_co_core > 0.:
        # co_core_mass, note there will be no kick
        m_rembar = m_co_core
    elif m_He_core > 0.:
        m_rembar = m_He_core
    else:
        # this is catching H-rich_non_burning stars
        if m_star < 0.5:
            m_rembar = m_star
            if ((m_co_core < 0.)or(m_He_core < 0.)):
                Pwarn('Invalid co/He core masses! '
                                'Setting m_WD=m_star!', "ApproximationWarning")
            else:
                Pwarn('co/He core masses are zero! '
                                'Setting m_WD=m_star!', "ApproximationWarning")
        else:
            raise ModelError('Invalid co/He core masses! Cannot complete SN.')
    f_fb = 1.0
    return m_rembar, f_fb, "WD"

# def NS_mass_Boccioli24(xi_1p5,xi_1p75,xi_2p0,xi_2p25):
#     if xi_1p5<0.6:
#         return 0.191*xi_1p5+1.242
#     elif xi_1p5>=0.6 and xi_1p75<0.6:
#         return 0.578*xi_1p75 + 1.237
#     elif xi_1p5>=0.6 and xi_1p75>=0.6 and xi_2p0<0.6:
#         return 0.914*xi_2p0+1.262
#     else:
#         return 0.665*xi_2p25+1.514
