import numpy as np
from abc import ABC, abstractmethod

from posydon.utils.common_functions import (  # CO_radius,; calculate_Patton20_values_at_He_depl,; inspiral_timescale_from_separation,; orbital_period_from_separation,; rotate,; separation_evol_wind_loss,; set_binary_to_failed,
    is_number,
)
from posydon.utils.posydonerror import ModelError
from posydon.utils.posydonwarning import Pwarn

SN_MECHANISM = [
    "Fryer+12-rapid",
    "Fryer+12-delayed",
    "direct",
    "direct_he_core",
    "Sukhbold+16-engine",
    "Patton&Sukhbold20-engine",
    "Couch+20-engine",
    "Maltsev+25-engine",
]

def get_CCSN_prescription(**kwargs):
    option = kwargs['mechanism']
    if option=="Fryer+12-rapid":
        return CCSN_check_Fryer12_rapid(**kwargs)
    elif option=="Fryer+12-delayed":
        return CCSN_check_Fryer12_delayed(**kwargs)
    elif option=="direct":
        return CCSN_check_direct(**kwargs)
    elif option=="direct_he_core":
        return CCSN_check_direct_he_core(**kwargs)
    elif option=="Sukhbold+16-engine":
        return CCSN_check_Sukhbold16(**kwargs)
    elif option=="Patton&Sukhbold20-engine":
        return CCSN_check_Patton20(**kwargs)
    elif option=="Couch+20-engine":
        return CCSN_check_Couch_20(**kwargs)
    elif option=="Maltsev+25-engine":
        return CCSN_check_Maltsev_25(**kwargs)
    else:
        raise ValueError(
            f"Invalid option {option} given to mechanism (i.e. the CCSN prescription). \n" + 
            "Valid options are {valid_options}.".format(
                valid_options="".join([x+", " for x in SN_MECHANISM[:-1]])+"or "+SN_MECHANISM[-1]
            )
        )


class CCSN_check_base(ABC):
    def __init__(self,**kwargs):
        self.mechanism = kwargs['mechanism']
    def __call__(self,star,**kwargs):
        raise ModelError("In step_SN, the explosion criterion is undefined")

class CCSN_check_Fryer12_rapid(CCSN_check_base):
    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        self.expl_crit_name = "Fryer+12-rapid"
    def __call__(self,star,**kwargs):
        return

class CCSN_check_Fryer12_delayed(CCSN_check_base):
    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        self.expl_crit_name = "Fryer+12-rapid"
    def __call__(self,star,**kwargs):
        return

class CCSN_check_direct(CCSN_check_base):
    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        self.expl_crit_name = "Fryer+12-delayed"
    def __call__(self,star,**kwargs):
        return

class CCSN_check_direct_he_core(CCSN_check_base):
    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        self.expl_crit_name = "Fryer+12-delayed"
    def __call__(self,star,**kwargs):
        return

class CCSN_check_Sukhbold16(CCSN_check_base):
    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        self.expl_crit_name = "Fryer+12-delayed"
    def __call__(self,star,**kwargs):
        return

class CCSN_check_Patton20(CCSN_check_base):
    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        self.expl_crit_name = "Fryer+12-delayed"
    def __call__(self,star,**kwargs):
        return

class CCSN_check_Couch_20(CCSN_check_base):
    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        self.expl_crit_name = "Fryer+12-delayed"
    def __call__(self,star,**kwargs):
        return

class CCSN_check_Maltsev_25(CCSN_check_base):
    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        self.expl_crit_name = "Fryer+12-delayed"
    def __call__(self,star,**kwargs):
        return
