from posydon.utils.posydonwarning import Pwarn
from posydon.utils.posydonerror import ModelError

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

def constant_NS_mass():
    return

def Mueller16_NSmass():
    return

def Boccioli24_NSmass():
    return

CO_MASS_CRITERIA = {
        "constant":constant_NS_mass,
        "Mueller16_NSmass":Mueller16_NSmass,
        "Boccioli24_NSmass":Boccioli24_NSmass,
}

def NS_mass_Boccioli24(xi_1p5,xi_1p75,xi_2p0,xi_2p25):
    if xi_1p5<0.6:
        return 0.191*xi_1p5+1.242
    elif xi_1p5>=0.6 and xi_1p75<0.6:
        return 0.578*xi_1p75 + 1.237
    elif xi_1p5>=0.6 and xi_1p75>=0.6 and xi_2p0<0.6:
        return 0.914*xi_2p0+1.262
    else:
        return 0.665*xi_2p25+1.514
