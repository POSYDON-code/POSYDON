"""
This file contains all core-collapse models
we use in our post processing and, as a 
consequence, what the user will be able to 
select in the initial/final interpolation.
Each model is a dictionary containing the 
properties of the model used by step_SN. 
"""

__authors__ = [
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Max Briel <max.briel@gmail.com>",
]

from posydon.utils.limits_thresholds import (STATE_NS_STARMASS_UPPER_LIMIT,
    NEUTRINO_MASS_LOSS_UPPER_LIMIT)

MODELS = {
    "MODEL01": {
        "mechanism": "direct",
        "engine": "",
        "PISN": "Marchant+19",
        "ECSN": "Podsiadlowski+04",
        "conserve_hydrogen_envelope" : False,
        "max_neutrino_mass_loss": NEUTRINO_MASS_LOSS_UPPER_LIMIT,
        "max_NS_mass": STATE_NS_STARMASS_UPPER_LIMIT,
        "use_interp_values": False,
        "use_profiles": True,
        "use_core_masses": False,
        "approx_at_he_depletion": False,
    },
    "MODEL02": {
        "mechanism": "Fryer+12-rapid",
        "engine": "",
        "PISN": "Marchant+19",
        "ECSN": "Podsiadlowski+04",
        "conserve_hydrogen_envelope" : False,
        "max_neutrino_mass_loss": NEUTRINO_MASS_LOSS_UPPER_LIMIT,
        "max_NS_mass": STATE_NS_STARMASS_UPPER_LIMIT,
        "use_interp_values": False,
        "use_profiles": True,
        "use_core_masses": False,
        "approx_at_he_depletion": False,
    },
    "MODEL03": {
        "mechanism": "Fryer+12-delayed",
        "engine": "",
        "PISN": "Marchant+19",
        "ECSN": "Podsiadlowski+04",
        "conserve_hydrogen_envelope" : False,
        "max_neutrino_mass_loss": NEUTRINO_MASS_LOSS_UPPER_LIMIT,
        "max_NS_mass": STATE_NS_STARMASS_UPPER_LIMIT,
        "use_interp_values": False,
        "use_profiles": True,
        "use_core_masses": False,
        "approx_at_he_depletion": False,
    },
    "MODEL04": {
        "mechanism": "Sukhbold+16-engine",
        "engine": "N20",
        "PISN": "Marchant+19",
        "ECSN": "Podsiadlowski+04",
        "conserve_hydrogen_envelope" : False,
        "max_neutrino_mass_loss": NEUTRINO_MASS_LOSS_UPPER_LIMIT,
        "max_NS_mass": STATE_NS_STARMASS_UPPER_LIMIT,
        "use_interp_values": False,
        "use_profiles": True,
        "use_core_masses": False,
        "approx_at_he_depletion": False,
    },
    "MODEL05": {
        "mechanism": "Patton&Sukhbold20-engine",
        "engine": "N20",
        "PISN": "Marchant+19",
        "ECSN": "Podsiadlowski+04",
        "conserve_hydrogen_envelope" : False,
        "max_neutrino_mass_loss": NEUTRINO_MASS_LOSS_UPPER_LIMIT,
        "max_NS_mass": STATE_NS_STARMASS_UPPER_LIMIT,
        "use_interp_values": False,
        "use_profiles": True,
        "use_core_masses": False,
        "approx_at_he_depletion": False,
    },
    "MODEL06": {
        "mechanism": "direct",
        "engine": "",
        "PISN": "Marchant+19",
        "ECSN": "Podsiadlowski+04",
        "conserve_hydrogen_envelope" : True,
        "max_neutrino_mass_loss": NEUTRINO_MASS_LOSS_UPPER_LIMIT,
        "max_NS_mass": STATE_NS_STARMASS_UPPER_LIMIT,
        "use_interp_values": False,
        "use_profiles": True,
        "use_core_masses": False,
        "approx_at_he_depletion": False,
    },
    "MODEL07": {
        "mechanism": "Fryer+12-rapid",
        "engine": "",
        "PISN": "Marchant+19",
        "ECSN": "Podsiadlowski+04",
        "conserve_hydrogen_envelope" : True,
        "max_neutrino_mass_loss": NEUTRINO_MASS_LOSS_UPPER_LIMIT,
        "max_NS_mass": STATE_NS_STARMASS_UPPER_LIMIT,
        "use_interp_values": False,
        "use_profiles": True,
        "use_core_masses": False,
        "approx_at_he_depletion": False,
    },
    "MODEL08": {
        "mechanism": "Fryer+12-delayed",
        "engine": "",
        "PISN": "Marchant+19",
        "ECSN": "Podsiadlowski+04",
        "conserve_hydrogen_envelope" : True,
        "max_neutrino_mass_loss": NEUTRINO_MASS_LOSS_UPPER_LIMIT,
        "max_NS_mass": STATE_NS_STARMASS_UPPER_LIMIT,
        "use_interp_values": False,
        "use_profiles": True,
        "use_core_masses": False,
        "approx_at_he_depletion": False,
    },
    "MODEL09": {
        "mechanism": "Sukhbold+16-engine",
        "engine": "N20",
        "PISN": "Marchant+19",
        "ECSN": "Podsiadlowski+04",
        "conserve_hydrogen_envelope" : True,
        "max_neutrino_mass_loss": NEUTRINO_MASS_LOSS_UPPER_LIMIT,
        "max_NS_mass": STATE_NS_STARMASS_UPPER_LIMIT,
        "use_interp_values": False,
        "use_profiles": True,
        "use_core_masses": False,
        "approx_at_he_depletion": False,
    },
    "MODEL10": {
        "mechanism": "Patton&Sukhbold20-engine",
        "engine": "N20",
        "PISN": "Marchant+19",
        "ECSN": "Podsiadlowski+04",
        "conserve_hydrogen_envelope" : True,
        "max_neutrino_mass_loss": NEUTRINO_MASS_LOSS_UPPER_LIMIT,
        "max_NS_mass": STATE_NS_STARMASS_UPPER_LIMIT,
        "use_interp_values": False,
        "use_profiles": True,
        "use_core_masses": False,
        "approx_at_he_depletion": False,
    },
    # "MODEL11": {
    #     "mechanism": "Patton&Sukhbold20-engine",
    #     "engine": "N20",
    #     "PISN": "Hendriks+23",
    #     "PISN_CO_shift": 0.0,
    #     "PPI_extra_mass_loss": 0.0,
    #     "ECSN": "Podsiadlowski+04",
    #     "conserve_hydrogen_envelope" : True,
    #     "max_neutrino_mass_loss": NEUTRINO_MASS_LOSS_UPPER_LIMIT,
    #     "max_NS_mass": STATE_NS_STARMASS_UPPER_LIMIT,
    #     "use_interp_values": False,
    #     "use_profiles": True,
    #     "use_core_masses": False,
    #     "approx_at_he_depletion": False,
    # },
    # "MODEL12": {
    #     "mechanism": "Patton&Sukhbold20-engine",
    #     "engine": "N20",
    #     "PISN": "Hendriks+23",
    #     "PISN_CO_shift": 0.0,
    #     "PPI_extra_mass_loss": 0.0,
    #     "ECSN": "Podsiadlowski+04",
    #     "conserve_hydrogen_envelope" : False,
    #     "max_neutrino_mass_loss": NEUTRINO_MASS_LOSS_UPPER_LIMIT,
    #     "max_NS_mass": STATE_NS_STARMASS_UPPER_LIMIT,
    #     "use_interp_values": False,
    #     "use_profiles": True,
    #     "use_core_masses": False,
    #     "approx_at_he_depletion": False,
    # },
}

def get_MODEL_NAME(input_MODEL, verbose=False):
    '''Find MODEL_NAME that matches input_MODEL.
    
    Parameters
    ----------
    input_MODEL : dict
        Dictionary with the properties of the model.
    
    Returns
    -------
    MODEL_NAME_SEL : str
        Name of the model in MODELS.py that matches input_MODEL.
    '''
    MODEL_NAME_SEL = None
    for MODEL_NAME, MODEL in MODELS.items():
        tmp = MODEL_NAME
        for key, val in MODEL.items():
            if "use_" in key or key=="ECSN":
                # escape values, which are allowed to differ
                continue
            if input_MODEL[key] != val:
                if verbose:
                    print(tmp, 'mismatch:', key,
                            input_MODEL[key], val)
                tmp = None
                break
        if tmp is not None:
            if verbose:
                print('matched to model:', tmp)
            MODEL_NAME_SEL = tmp
            
    return MODEL_NAME_SEL