from abc import ABC, abstractmethod
import numbers

def get_PISN_prescription(**kwargs):
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
        raise ValueError("Invalid option {option} given to PISN prescription".format(option=option))

class PISN_check_base(ABC):
    def __init__(self,**kwargs):
        self.PISN_option = kwargs['PISN']
        self.verbose = kwargs['verbose']

    @abstractmethod
    def __call__(self,star):
        raise NotImplementedError

    def __repr__(self):
        return "Prescription_PISN was initialised with the {option} prescription".format(option=self.PISN_option)

class PISN_check_Marchant19(PISN_check_base):
    def __call__(self,star):
        m_He_core = star.he_core_mass
        if m_He_core >= 31.99 and m_He_core <= 61.10:
            return "PPISN"
        if m_He_core > 61.10 and m_He_core < 124.12:
            return "PISN"
        # above the gap (m_He_core >= 124.12) the star is a CCSN, as on main
        return None

class PISN_check_Hendriks23(PISN_check_base):
    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        self.PISN_CO_shift = kwargs.get('PISN_CO_shift')
        self.delta_M_CO_shift = self.PISN_CO_shift if self.PISN_CO_shift is not None else 0.0
        self.m_CO_core_PISN_min = 38 + self.delta_M_CO_shift
        self.m_CO_core_PISN_max = 114 + self.delta_M_CO_shift
    def __call__(self,star):
        m_CO_core = star.co_core_mass
        if m_CO_core >= self.m_CO_core_PISN_min and m_CO_core <= self.m_CO_core_PISN_max:
            return "PPISN"
        # above the gap the star is a CCSN, as on main
        return None

class PISN_check_constant(PISN_check_base):
    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        self.PISN_option = float(kwargs['PISN'])
    def __call__(self,star):
        if star.he_core_mass>self.PISN_option:
            return "PPISN"
        return None

class PISN_check_None(PISN_check_base):
    def __call__(self,star):
        return None