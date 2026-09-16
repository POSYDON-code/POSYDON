import numpy as np

from posydon.utils.posydonerror import ModelError

from posydon.utils.common_functions import (
    # CO_radius,
    # calculate_Patton20_values_at_He_depl,
    # inspiral_timescale_from_separation,
    is_number,
    # orbital_period_from_separation,
    # rotate,
    # separation_evol_wind_loss,
    # set_binary_to_failed,
)
from posydon.utils.posydonwarning import Pwarn

class BaseCCSNCrit:
    def __call__(self,binary,**kwargs):
        raise ModelError("In step_SN, the explosion criterion is undefined")

class Fryer12_rapid(BaseCCSNCrit):
    def __init__(self,**kwargs):
        # super().__init__(**kwargs)
        self.expl_crit_name = "Fryer+12-rapid"
    def __call__(self,binary,**kwargs):
        



def Fryer12_delayed():
    return

def direct():
    return

def direct_he_core():
    return

def Sukhbold16_engine():
    return

def Patton20_engine():
    return

def Couch20_engine():
    return

def Maltsev25_detailed():
    return

def Maltsev25_rapid():
    return

def ccsn_user_defined():
    return

EXPLODABILITY_CRITERIA = {
    "Fryer+12-rapid":Fryer12_rapid,
    "Fryer+12-delayed":Fryer12_delayed,
    "direct":direct,
    "direct_he_core":direct_he_core,
    "Sukhbold+16-engine":Sukhbold16_engine,
    "Patton&Sukhbold20-engine":Patton20_engine,
    "Couch+20-engine":Couch20_engine,
    "Maltsev+25-engine":Maltsev25_detailed,
    "Maltsev+25-MCO-rapid":Maltsev25_rapid,
    "User_defined_CCSN":ccsn_user_defined,
}