min_M_CO_mass_for_ECSN = {
    "Tauris+15":1.37,        # Msun from Takahashi et al. (2013)
    "Podsiadlowski+04":1.4,  # Msun from Podsiadlowski+2004
    "No_ECSN":1.37,          # Msun from Takahashi et al. (2013)
}

max_M_CO_mass_for_ECSN = {
    "Tauris+15":1.43,        # Msun from Tauris et al. (2015)
    "Podsiadlowski+04":2.5,  # Msun from Podsiadlowski et al. (2015)
    "No_ECSN":1.37,          # Msun from Takahashi et al. (2013)
}


class ECSN_check:
    # DAVID (2026.09.15): All of the currently included ECSN prescriptions work the same way,
    # so I put them all in one function. I had to add an extra if clause to include the 
    # No_ECSN case. Could be avoided by switching either the <= or the >= in the second elif.
    # Please advise!
    # Also: Add verbosity?
    def __init__(self,ECSN_option):
        self.ECSN_option = ECSN_option
        self.min_M_CO_ECSN = min_M_CO_mass_for_ECSN[ECSN_option]
        self.max_M_CO_ECSN = max_M_CO_mass_for_ECSN[ECSN_option]
    def __call__(self,star):
        m_co_core = star.co_core_mass
        if m_co_core < self.min_M_CO_ECSN:
            return "WD"
        elif (m_co_core >= self.min_M_CO_ECSN) and (m_co_core <= self.max_M_CO_ECSN):
            if self.ECSN_option=="No_ECSN":
                return "CCSN"
            else:
                return "ECSN"
        elif m_co_core > self.max_M_CO_ECSN:
            return "CCSN"
        else:
            raise ValueError(
                "The SN step was applied for an on object outside the "
                "domain of electron-capture SN and Fe core-collapse SN."
            )


# DAVID (2026.09.15): Still don't know how to implement this exactly...
def ecsn_user_defined():
    return #m_rembar, f_fb, state, SN_type
