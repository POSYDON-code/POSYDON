"""Pair-instability checks.

Decide whether a star is a candidate for a pulsational pair-instability
supernova ("PPISN") or a pair-instability supernova ("PISN"), and whether it
is above the pair-instability gap. The final PPISN/PISN label is settled by
the remnant-mass stage.
"""

import numbers
from abc import ABC, abstractmethod


def get_PISN_prescription(**kwargs):
    """Return the PISN check given by kwargs['PISN']."""
    option = kwargs['PISN']
    if option=="Marchant+19":
        return PISN_check_Marchant19(**kwargs)
    elif option=="Hendriks+23":
        return PISN_check_Hendriks23(**kwargs)
    elif isinstance(option,numbers.Real) and not isinstance(option,bool):
        return PISN_check_constant(**kwargs)
    elif option is None:
        return PISN_check_None(**kwargs)
    else:
        raise ValueError(f"Invalid option {option} given to PISN "
                         "prescription")

class PISN_check_base(ABC):
    """Base class of the PISN checks.

    Calling it returns "PPISN", "PISN" or None; subclasses implement
    _classify and, if the prescription has one, the upper end of the gap in
    is_above_gap.
    """
    def __init__(self,**kwargs):
        self.PISN_option = kwargs['PISN']
        self.verbose = kwargs['verbose']

    def __call__(self,star):
        """Return "PPISN" or "PISN" for a candidate, None otherwise."""
        outcome = self._classify(star)
        if self.verbose:
            print(f"{type(self).__name__}: He core {star.he_core_mass:.3f} "
                  f"Msun, CO core {star.co_core_mass:.3f} Msun: "
                  f"{outcome or 'not a PPISN or PISN candidate'}")
        return outcome

    @abstractmethod
    def _classify(self,star):
        """Return "PPISN", "PISN" or None for the star."""
        raise NotImplementedError

    def is_above_gap(self,star):
        """True if the star is above the pair-instability gap.

        Such stars collapse directly into a BH (explosion outcome "BH"),
        whatever the CCSN mechanism.
        """
        return False

    def __repr__(self):
        return (f"Prescription_PISN was initialised with the "
                f"{self.PISN_option} prescription")

class PISN_check_Marchant19(PISN_check_base):
    """He-core mass ranges of Marchant et al. (2019), ApJ, 882, 36.

    PPISN candidates for 31.99-61.10 Msun, PISN for 61.10-124.12 Msun, and
    direct collapse above.
    """
    def _classify(self,star):
        m_He_core = star.he_core_mass
        if m_He_core >= 31.99 and m_He_core <= 61.10:
            return "PPISN"
        if m_He_core > 61.10 and m_He_core < 124.12:
            return "PISN"
        return None

    def is_above_gap(self,star):
        # direct collapse to a BH above the gap, see table 1 of Marchant et
        # al. (2019)
        return star.he_core_mass >= 124.12

class PISN_check_Hendriks23(PISN_check_base):
    """CO-core mass range of Hendriks et al. (2023), MNRAS, 526, 4130.

    PPISN candidates for 38-114 Msun (shifted by PISN_CO_shift), and direct
    collapse above. A candidate becomes a PISN in the remnant-mass stage if
    the pulses leave less than 10 Msun.
    """
    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        self.PISN_CO_shift = kwargs['PISN_CO_shift']
        self.delta_M_CO_shift = (self.PISN_CO_shift
                                 if self.PISN_CO_shift is not None else 0.0)
        self.m_CO_core_PISN_min = 38 + self.delta_M_CO_shift
        self.m_CO_core_PISN_max = 114 + self.delta_M_CO_shift
    def _classify(self,star):
        m_CO_core = star.co_core_mass
        if self.m_CO_core_PISN_min <= m_CO_core <= self.m_CO_core_PISN_max:
            return "PPISN"
        return None

    def is_above_gap(self,star):
        # "For CO core masses greater than 114 Msun we assume direct collapse
        # to a BH" (Hendriks et al. 2023)
        return star.co_core_mass > self.m_CO_core_PISN_max

class PISN_check_constant(PISN_check_base):
    """He cores above the cap (option PISN, in Msun) are PPISN candidates;
    their remnant mass is capped at that value."""
    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        self.PISN_option = float(kwargs['PISN'])
    def _classify(self,star):
        if star.he_core_mass>self.PISN_option:
            return "PPISN"
        return None

class PISN_check_None(PISN_check_base):
    """No pair instability."""
    def _classify(self,star):
        return None
