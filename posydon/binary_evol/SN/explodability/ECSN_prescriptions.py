from abc import ABC, abstractmethod

def get_ECSN_prescription(**kwargs):
    option = kwargs['ECSN']
    if option=="Tauris+15":
        return ECSN_check_Tauris(**kwargs)
    elif option=="Podsiadlowski+04":
        return ECSN_check_Podsiadlowski(**kwargs)
    elif option=="No_ECSN" or (option is None):
        return ECSN_check_NoECSN(**kwargs)
    else:
        raise ValueError(
            f"Invalid option {option} given to ECSN prescription. "
            "Valid options are 'Tauris+15', 'Podsiadlowski+04', 'No_ECSN' or None."
        )

class ECSN_check_base(ABC):
    def __init__(self,**kwargs):
        self.ECSN_option = kwargs['ECSN']
        self.verbose = kwargs['verbose']
    
    def __call__(self,star):
        m_core = self.core_mass(star)
        if m_core < self.min_core_mass_ECSN:
            if self.verbose:
                print(f'The ECSN_check routine has determined that a collapsing star with a {self.ECSN_core_type} core of')
                print('{core_mass} Msun will leave a WD remnant.'.format(core_mass=m_core))
            return "WD"
        if m_core <= self.max_core_mass_ECSN:
            if self.verbose:
                print(f'The ECSN_check routine has determined that a collapsing star with a {self.ECSN_core_type} core of')
                print('{core_mass} Msun is an ECSN (min_core_mass_ECSN = {min_CO_ECSN}, max_core_mass_ECSN={max_CO_ECSN})'.format(
                    core_mass=m_core,
                    min_CO_ECSN = self.min_core_mass_ECSN,
                    max_CO_ECSN = self.max_core_mass_ECSN,
                    ))
            return "ECSN"
        if self.verbose:
            print(f'The ECSN_check routine has determined that a collapsing star with a {self.ECSN_core_type} core of')
            print('{core_mass} Msun is not an ECSN or a WD'.format(core_mass=m_core))
        return None
        
    @abstractmethod
    def core_mass(self,star):
        raise NotImplementedError
    
    def __repr__(self):
        return "Prescription_ECSN was initialised with the {option} prescription".format(option=self.ECSN_option)


class ECSN_check_Tauris(ECSN_check_base):
    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        self.ECSN_core_type = "CO"
        self.min_core_mass_ECSN = 1.37        # Msun from Takahashi et al. (2013)
        self.max_core_mass_ECSN = 1.43        # Msun from Tauris et al. (2015)
    def core_mass(self,star):
        return star.co_core_mass


class ECSN_check_Podsiadlowski(ECSN_check_base):
    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        self.ECSN_core_type = "He"
        self.min_core_mass_ECSN = 1.4            # Msun from Podsiadlowski et al. (2004)
        self.max_core_mass_ECSN = 2.5            # Msun from Podsiadlowski et al. (2004)
    def core_mass(self,star):
        return star.he_core_mass

class ECSN_check_NoECSN(ECSN_check_base):
    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        self.ECSN_core_type = "CO"
        self.min_core_mass_ECSN = 1.37           # Msun from Takahashi et al. (2013)
        self.max_core_mass_ECSN = float("-inf")  # Guarantees that there are no ECSN
    def core_mass(self,star):
        return star.co_core_mass

