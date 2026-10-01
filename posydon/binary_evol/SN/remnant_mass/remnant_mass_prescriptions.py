import numpy as np
from posydon.utils.posydonwarning import Pwarn
from posydon.utils.posydonerror import ModelError
from abc import ABC, abstractmethod

def get_remnant_mass_prescription(**kwargs):
    option = kwargs["NS_mass_prescription"]
    if option=="constant" or isinstance(option,(int,float)):
        return Remnant_Mass_constant(**kwargs)
    elif option=="Fryer+12-delayed":
        return Remnant_Mass_Fryer_delayed(**kwargs)
    elif option=="Fryer+12-rapid":
        return Remnant_Mass_Fryer_rapid(**kwargs)
    elif option=="Sukhbold+16":
        return Remnant_Mass_Sukhbold16(**kwargs)
    elif option=="Patton&Sukhbold20":
        return Remnant_Mass_Patton20(**kwargs)
    elif option=="Mueller+16":
        return Remnant_Mass_Mueller16(**kwargs)
    elif option=="Boccioli&Fragione24":
        return Remnant_Mass_Boccioli24(**kwargs)
    else:
        raise ValueError("Invalid option {option} given to remnant mass prescription".format(option=option))

def ECSN_mass_Tauris(star):
    return star

def ECSN_mass_Podsiadlowski(star):
    return star

def ECSN_mass_No_ECSN(star):
    return star

ECSN_MASS_OPTIONS = {
    "Tauris+15":ECSN_mass_Tauris,
    "Podsiadlowski+04":ECSN_mass_Podsiadlowski,
    "No_ECSN":ECSN_mass_No_ECSN,
}

def PISN_mass_Marchant(star):
    m_He_core = star.he_core_mass
    if m_He_core >= 31.99 and m_He_core <= 61.10:
        # this is the 8th-order polynomial fit of table 1
        # value, see COSMIC paper (Breivik et al. 2020)
        polyfit = (
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
        return polyfit
    elif m_He_core > 61.10 and m_He_core < 124.12:
        # in Breivik et al. (2020) they qoute the CO core mass
        # range as 54.48<M_CO-core/Msun<113.29 here, but this
        # might cause gaps, when switching between core masses,
        # hence take the He-core masses from table 1 of Marchant
        # et al. (2019)
        return np.nan
    
def PISN_mass_Hendriks(star):
    
    return star

PISN_MASS_OPTIONS = {
    "Marchant+19":PISN_mass_Marchant,
    "Hendriks+23":PISN_mass_Hendriks,
}

class Remnant_Mass_Base(ABC):
    def __init__(self,**kwargs):
        self.ECSN_option = kwargs['ECSN']
        self.PISN_option = kwargs['PISN']
        self.conserve_hydrogen_envelope = kwargs['conserve_hydrogen_envelope']
        self.verbose = kwargs['verbose']
        self._get_rem_mass_WD = WD_mass
        self._get_rem_mass_ECSN = ECSN_MASS_OPTIONS[self.ECSN_option]
        self._get_rem_mass_PISN = PISN_MASS_OPTIONS[self.PISN_option]
    def __call__(self,star,end_class):
        if end_class=="WD":
            return self._get_rem_mass_WD(star)
        elif end_class=="ECSN":
            return self._get_rem_mass_ECSN(star)
        elif end_class=="CCSN":
            return self._get_rem_mass_CCSN(star)
        elif end_class=="PISN":
            return self._get_rem_mass_PISN(star)
        else:
            raise ValueError("Invalid end class {end_class} passed during remnant mass calculation".format(end_class=end_class))

    def _get_rem_mass_direct_collapse(self,star):
            if self.conserve_hydrogen_envelope:
                return "BH", 1.0, star.mass
            else:
                return "BH", 1.0, star.he_core_mass

    @abstractmethod
    def _get_rem_mass_CCSN(self,star):
        return NotImplementedError

    @abstractmethod
    def __repr__(self):
        raise NotImplementedError


class Remnant_Mass_constant(Remnant_Mass_Base):
    def __init__(**kwargs):
        super().__init__(**kwargs)

class Remnant_Mass_Fryer_delayed(Remnant_Mass_Base):
    def __init__(**kwargs):
        super().__init__(**kwargs)

class Remnant_Mass_Fryer_rapid(Remnant_Mass_Base):
    def __init__(**kwargs):
        super().__init__(**kwargs)

class Remnant_Mass_Sukhbold16(Remnant_Mass_Base):
    def __init__(**kwargs):
        super().__init__(**kwargs)

class Remnant_Mass_Patton20(Remnant_Mass_Base):
    def __init__(**kwargs):
        super().__init__(**kwargs)

class Remnant_Mass_Mueller16(Remnant_Mass_Base):
    def __init__(**kwargs):
        super().__init__(**kwargs)

class Remnant_Mass_Boccioli24(Remnant_Mass_Base):
    def __init__(**kwargs):
        super().__init__(**kwargs)




# DAVID (2026.09.15): Add functionality for checking if WD, ECSN, PISN or CCSN
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
    return m_rembar, f_fb

def NS_mass_Boccioli24(xi_1p5,xi_1p75,xi_2p0,xi_2p25):
    if xi_1p5<0.6:
        return 0.191*xi_1p5+1.242
    elif xi_1p5>=0.6 and xi_1p75<0.6:
        return 0.578*xi_1p75 + 1.237
    elif xi_1p5>=0.6 and xi_1p75>=0.6 and xi_2p0<0.6:
        return 0.914*xi_2p0+1.262
    else:
        return 0.665*xi_2p25+1.514
