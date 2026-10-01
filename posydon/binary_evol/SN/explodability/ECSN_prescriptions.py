from abc import ABC, abstractmethod

def get_ECSN_prescription(**kwargs):
    option = kwargs['ECSN']
    if option=="Tauris+15":
        return ECSN_check_Tauris(**kwargs)
    elif option=="Podsiadlowski+04":
        return ECSN_check_Podsiadlowski(**kwargs)
    elif option=="No_ECSN":
        return ECSN_check_NoECSN(**kwargs)
    else:
        raise ValueError("Invalid option {option} given to ECSN prescription".format(option=option))

class ECSN_check_base(ABC):
    def __init__(self,**kwargs):
        self.ECSN_option = kwargs['ECSN']
        self.verbose = kwargs['verbose']
    @abstractmethod
    def __call__(self,star):
        raise NotImplementedError
    
    @abstractmethod
    def __repr__(self):
        raise NotImplementedError

    def _ECSN_check(self,star):
        m_co_core = star.co_core_mass
        if m_co_core < self.min_M_CO_ECSN:
            if self.verbose:
                print('The ECSN_check routine has determined that a collapsing star with a CO core of')
                print('{core_mass} Msun will leave a WD remnant.'.format(core_mass=m_co_core))
            return "WD"
        elif (m_co_core >= self.min_M_CO_ECSN) and (m_co_core <= self.max_M_CO_ECSN):
            if self.ECSN_option=="No_ECSN":
                if self.verbose:
                    print('The ECSN_check routine has determined that a collapsing star with a CO core of')
                    print('{core_mass} Msun is not an ECSN (option ECSN is "No_ECSN") or a WD'.format(core_mass=m_co_core))
                return
            else:
                if self.verbose:
                    print('The ECSN_check routine has determined that a collapsing star with a CO core of')
                    print('{core_mass} Msun is an ECSN (min_M_CO_ECSN = {min_CO_ECSN}, max_M_CO_ECSN={max_CO_ECSN})'.format(
                        core_mass=m_co_core,
                        min_CO_ECSN = self.min_M_CO_ECSN,
                        max_CO_ECSN = self.max_M_CO_ECSN,
                        ))
                return "ECSN"
        elif m_co_core > self.max_M_CO_ECSN:
            if self.verbose:
                    print('The ECSN_check routine has determined that a collapsing star with a CO core of')
                    print('{core_mass} Msun is not an ECSN (option ECSN is "No_ECSN") or a WD'.format(core_mass=m_co_core))
            return
        else:
            raise ValueError(
                "Prescription_ECSN has encountered an invalid value of the collapsing star's ",
                "C/O core mass. m_co_core = {m_co_core}, ECSN option = {option},".format(m_co_core=m_co_core,option=self.ECSN_option),
                "min_M_CO_mass_for_ECSN = {min_M_CO}, max_M_CO_mass_for_ECSN = {max_M_CO}".format(min_M_CO=self.min_M_CO_ECSN,max_M_CO=self.max_M_CO_ECSN)
            )


class ECSN_check_Tauris(ECSN_check_base):
    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        self.min_M_CO_ECSN = 1.37        # Msun from Takahashi et al. (2013)
        self.max_M_CO_ECSN = 1.43        # Msun from Tauris et al. (2015)
    def __repr__(self):
        return "Prescription_ECSN was initialised with the {option} prescription".format(option=self.ECSN_option)
    def __call__(self,star):
        return self._ECSN_check(star)

class ECSN_check_Podsiadlowski(ECSN_check_base):
    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        self.min_M_CO_ECSN = 1.4         # Msun from Podsiadlowski+2004
        self.max_M_CO_ECSN = 2.5         # Msun from Podsiadlowski et al. (2015)
    def __repr__(self):
        return "Prescription_ECSN was initialised with the {option} prescription".format(option=self.ECSN_option)
    def __call__(self,star):
        return self._ECSN_check(star)

class ECSN_check_NoECSN(ECSN_check_base):
    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        self.min_M_CO_ECSN = 1.37        # Msun from Takahashi et al. (2013)
        self.max_M_CO_ECSN = 1.37        # Msun from Tauris et al. (2015)
    def __repr__(self):
        return "Prescription_ECSN was initialised with the {option} prescription".format(option=self.ECSN_option)
    def __call__(self,star):
        m_co_core = star.co_core_mass
        if m_co_core <= self.min_M_CO_ECSN:
            if self.verbose:
                print('The ECSN_check routine has determined that a collapsing star with a CO core of')
                print('{core_mass} Msun will leave a WD remnant.'.format(core_mass=m_co_core))
            return "WD"
        elif m_co_core > self.max_M_CO_ECSN:
            if self.verbose:
                    print('The ECSN_check routine has determined that a collapsing star with a CO core of')
                    print('{core_mass} Msun is not an ECSN (option ECSN is "No_ECSN") or a WD'.format(core_mass=m_co_core))
            return
        else:
            raise ValueError(
                "Prescription_ECSN has encountered an invalid value of the collapsing star's ",
                "C/O core mass. m_co_core = {m_co_core}, ECSN option = {option},".format(m_co_core=m_co_core,option=self.ECSN_option),
                "min_M_CO_mass_for_ECSN = {min_M_CO}, max_M_CO_mass_for_ECSN = {max_M_CO}".format(min_M_CO=self.min_M_CO_ECSN,max_M_CO=self.max_M_CO_ECSN)
            )
