import numpy as np

        

class PISN_check:
    def __init__(self,PISN_option):
        if PISN_option=="Marchant+19":
            self.criterion = _pisn_marchant19
        elif PISN_option=="Hendriks+23":
            self.criterion = _pisn_hendricks23

def _pisn_marchant19(star,verbose=False):
    m_He_core = star.he_core_mass
    m_star = star.mass
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
        m_PISN = polyfit
    elif m_He_core > 61.10 and m_He_core < 124.12:
        # in Breivik et al. (2020) they qoute the CO core mass
        # range as 54.48<M_CO-core/Msun<113.29 here, but this
        # might cause gaps, when switching between core masses,
        # hence take the He-core masses from table 1 of Marchant
        # et al. (2019)
        m_PISN = np.nan

    else:
        # above the PISN gap we assume direct collapse of the
        # entire star, or the He core, into a BH.
        if self.conserve_hydrogen_envelope:
            m_PISN = m_star
        else:
            m_PISN = m_He_core
    if verbose:
        pisn_verbose(m_PISN)
    return m_PISN

def _pisn_hendricks23(self,star):
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
        pisn_verbose(m_PISN)
    return m_PISN

def pisn_constant(star):
    m_He_core = star.he_core_mass
    if m_He_core > self.PISN:
        m_PISN = self.PISN
    elif 0.0 < m_He_core <= self.PISN:
        m_PISN = None
    if self.verbose:
        pisn_verbose(m_PISN)
    return m_PISN

def pisn_user_defined():
    return

def pisn_verbose(m_PISN):
    if m_PISN is None:
        print("")
        print("The star did NOT lose any mass because of "
                "PPIN or PISN.")
    elif not pd.isna(m_PISN):
        print("")
        print(
            "The star with initial mass {:2.2f}".format(m_He_core),
            "M_sun went through the PISN routine and lost",
            "{:2.2f} M_sun.".format(m_He_core - m_PISN),
            "The new m_rembar mass that will collapse to form a ",
            "CO object is {:2.2f} M_sun.".format(m_PISN))
    else:
        print("The star was disrupted by the PISN prescription!")

PISN_PRESCRIPTIONS = {
    "Marchant+19":pisn_marchant19,
    "Hendriks+23":pisn_hendricks23,
    "User_defined_PISN":pisn_user_defined,
}
