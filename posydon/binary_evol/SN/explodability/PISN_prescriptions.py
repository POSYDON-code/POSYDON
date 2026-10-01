import numpy as np
import copy
from abc import ABC, abstractmethod
import posydon.utils.constants as const

def get_PISN_prescription(**kwargs):
    option = kwargs['PISN_option']
    if option=="Marchant+19":
        return PISN_check_marchant19(**kwargs)
    elif option=="Hendriks+23":
        return PISN_check_hendricks23(**kwargs)
    else:
        raise ValueError("Invalid option {option} given to PISN prescription".format(option=option))

class PISN_check_base(ABC):
    def __init__(self,**kwargs):
        self.PISN_option = kwargs['PISN']
        self.verbose = kwargs['verbose']
        self.PISN_CO_shift = kwargs['PISN_CO_shift']
        self.PPI_extra_mass_loss = kwargs['PPI_extra_mass_loss']
        self.conserve_hydrogen_PPI = kwargs['conserve_hydrogen_PPI']

    @abstractmethod
    def __call__(self,star):
        raise NotImplementedError

    @abstractmethod
    def __repr__(self):
        raise NotImplementedError
    
    def _pisn_verbose(self,m_He_core,m_PISN):
        if m_PISN is None:
            print("")
            print("The star did NOT lose any mass because of "
                    "PPIN or PISN.")
        elif not np.isna(m_PISN):
            print("")
            print(
                "The star with initial mass {:2.2f}".format(m_He_core),
                "M_sun went through the PISN routine and lost",
                "{:2.2f} M_sun.".format(m_He_core - m_PISN),
                "The new m_rembar mass that will collapse to form a ",
                "CO object is {:2.2f} M_sun.".format(m_PISN))
        else:
            print("The star was disrupted by the PISN prescription!")

        def _pisn_constant(self,star):
            m_He_core = star.he_core_mass
            if m_He_core > self.PISN:
                m_PISN = self.PISN
            elif 0.0 < m_He_core <= self.PISN:
                m_PISN = None
            if self.verbose:
                self._pisn_verbose(m_PISN)
            return m_PISN


class PISN_check_marchant19(PISN_check_base):
    def __call__(self,star):
        m_He_core = star.he_core_mass
        if m_He_core >= 31.99 and m_He_core < 124.12:
            return "PISN"
        elif m_He_core>=124.12:
            return "direct_collapse"
    def __repr__(self):
        return "Prescription_PISN was initialised with the {option} prescription".format(option=self.PISN_option)

class PISN_check_hendricks23(PISN_check_base):
    def __call__(self,star):
        # Hendriks et al. 2023 PISN prescription
        # 10.1093/mnras/stad2857
        # Shifting PPI and PISN gap
        # works by removing delta_M_PPI from the star
        # and then applying any remnant mass prescription
        m_He_core = star.he_core_mass
        m_CO_core = star.co_core_mass
        m_star = star.mass

        delta_M_CO_shift = self.PISN_CO_shift if self.PISN_CO_shift is not None else 0.0
        delta_M_PPI_extra_ML = self.PPI_extra_mass_loss if self.PPI_extra_mass_loss is not None else 0.0

        m_CO_core_PISN_min = 38 + delta_M_CO_shift
        m_CO_core_PISN_max = 114 + delta_M_CO_shift

        if ((m_CO_core >= m_CO_core_PISN_min)
            and m_CO_core <= m_CO_core_PISN_max):

            # delta_PPI -> -inf if Z -> 0
            # limit mass loss to Z = 1e-4 for Z below it.
            # 1e-4 is the lowest metallicity in the Hendriks et al. 2023
            if star.metallicity < 1e-4:
                Z = 1e-4
            else:
                Z = star.metallicity
            # Hendriks et al. 2023 Equation 6
            # 10.1093/mnras/stad2857
            delta_M_PPI = (
                (0.0006 * np.log10(Z * const.Zsun) + 0.0054)
                * (m_CO_core - delta_M_CO_shift - 34.8)**3
                - 0.0013 * (m_CO_core - delta_M_CO_shift - 34.8)**2
                + delta_M_PPI_extra_ML
            )
            if self.verbose:
                print(f"delta_M_PPI: {delta_M_PPI} Msun")
        else:
            delta_M_PPI = 0.0

        if delta_M_PPI <= 0.0:
            # no PPI -> use CCSN prescription
            if self.conserve_hydrogen_envelope:
                m_PISN = m_star
            else:
                m_PISN = m_He_core
        else:
            # PPI occurs
            if self.conserve_hydrogen_PPI:
                m_PISN = m_star - delta_M_PPI
            else:
                m_PISN = m_He_core - delta_M_PPI

            if m_PISN < 0.0:
                m_PISN = np.nan
            else:
                PISN_star = copy.deepcopy(star)
                PISN_star.mass = m_PISN
                if PISN_star.he_core_mass > m_PISN:
                    PISN_star.he_core_mass = m_PISN
                if PISN_star.co_core_mass > m_PISN:
                    PISN_star.co_core_mass = m_PISN
                m_rembar, _, _ = self.compute_m_rembar(PISN_star, m_PISN)

                if m_rembar < 10:
                    m_PISN = np.nan
                else:
                    m_PISN = m_rembar
        if self.verbose:
            self._pisn_verbose(m_He_core,m_PISN)
        return m_PISN
    
    def __repr__(self):
        return "Prescription_PISN was initialised with the {option} prescription".format(option=self.PISN_option)
