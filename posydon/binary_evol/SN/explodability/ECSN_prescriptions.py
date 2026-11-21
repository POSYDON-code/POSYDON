"""Electron-capture supernova checks.

From a core mass of the star, decide whether it leaves a WD, explodes as an
electron-capture supernova ("ECSN"), or neither (None).
"""

from abc import ABC, abstractmethod


def get_ECSN_prescription(**kwargs):
    """Return the ECSN check given by kwargs['ECSN']."""
    option = kwargs['ECSN']
    if option=="Tauris+15":
        return ECSN_check_Tauris15(**kwargs)
    elif option=="Podsiadlowski+04":
        return ECSN_check_Podsiadlowski04(**kwargs)
    elif option=="No_ECSN" or (option is None):
        return ECSN_check_None(**kwargs)
    else:
        raise ValueError(
            f"Invalid option {option} given to ECSN prescription. Valid "
            "options are 'Tauris+15', 'Podsiadlowski+04', 'No_ECSN' or None."
        )

class ECSN_check_base(ABC):
    """Base class of the ECSN checks.

    A star with a core below min_core_mass_ECSN leaves a WD, one with a core
    up to max_core_mass_ECSN is an ECSN, and any other gives None.
    Subclasses set the limits, and core_mass gives the core mass compared
    with them.
    """
    def __init__(self,**kwargs):
        self.ECSN_option = kwargs['ECSN']
        self.verbose = kwargs['verbose']

    def __call__(self,star):
        """Return "WD", "ECSN" or None for the star."""
        m_core = self.core_mass(star)
        if m_core < self.min_core_mass_ECSN:
            self._report(f'{m_core} Msun will leave a WD remnant.')
            return "WD"
        if m_core <= self.max_core_mass_ECSN:
            self._report(f'{m_core} Msun is an ECSN (min_core_mass_ECSN = '
                         f'{self.min_core_mass_ECSN}, max_core_mass_ECSN='
                         f'{self.max_core_mass_ECSN})')
            return "ECSN"
        self._report(f'{m_core} Msun is not an ECSN or a WD')
        return None

    def _report(self,outcome):
        """Print the outcome of the check, if verbose."""
        if self.verbose:
            print('The ECSN_check routine has determined that a collapsing '
                  f'star with a {self.ECSN_core_type} core of')
            print(outcome)

    @abstractmethod
    def core_mass(self,star):
        """Return the core mass (Msun) compared with the limits."""
        raise NotImplementedError

    def __repr__(self):
        return (f"Prescription_ECSN was initialised with the "
                f"{self.ECSN_option} prescription")


class ECSN_check_Tauris15(ECSN_check_base):
    """ECSN for CO cores of 1.37 Msun (Takahashi et al. 2013) to 1.43 Msun
    (Tauris et al. 2015, MNRAS, 451, 2123)."""
    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        self.ECSN_core_type = "CO"
        self.min_core_mass_ECSN = 1.37     # Msun from Takahashi et al. (2013)
        self.max_core_mass_ECSN = 1.43     # Msun from Tauris et al. (2015)
    def core_mass(self,star):
        return star.co_core_mass


class ECSN_check_Podsiadlowski04(ECSN_check_base):
    """ECSN for He cores of 1.4 to 2.5 Msun (Podsiadlowski et al. 2004, ApJ,
    612, 1044)."""
    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        self.ECSN_core_type = "He"
        # Msun, from Podsiadlowski et al. (2004)
        self.min_core_mass_ECSN = 1.4
        self.max_core_mass_ECSN = 2.5
    def core_mass(self,star):
        return star.he_core_mass

class ECSN_check_None(ECSN_check_base):
    """No ECSN: CO cores below 1.37 Msun leave a WD, all others give None."""
    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        self.ECSN_core_type = "CO"
        self.min_core_mass_ECSN = 1.37     # Msun from Takahashi et al. (2013)
        self.max_core_mass_ECSN = float("-inf")     # guarantees no ECSN
    def core_mass(self,star):
        return star.co_core_mass
