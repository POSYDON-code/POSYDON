"""
This file contains all supernova models we use in our post processing and, as a
consequence, what the user will be able to select in the initial/final
interpolation. There is a default supernova model, which is used to fill not
specified parameters in the set of models. Each model is a dictionary
containing the properties of the model used by step_SN. The functions
`get_SN_MODEL` provides the supernova model with all parameters and
`get_SN_MODEL_NAME` findes the matching model to a given set of supernova
parameters of the step_SN.
"""

__authors__ = [
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>",
    "Max Briel <max.briel@gmail.com>",
]

from posydon.utils.limits_thresholds import (STATE_NS_STARMASS_UPPER_LIMIT,
                                             NEUTRINO_MASS_LOSS_UPPER_LIMIT)

DEFAULT_SN_MODEL = {
    "mechanism": "Fryer+12-delayed",
    "engine": "",
    "PISN": "Hendriks+23",
    "PISN_CO_shift": 0.0,
    "PPI_extra_mass_loss": -20.0,
    "ECSN": "Tauris+15",
    "conserve_hydrogen_envelope" : False,
    "conserve_hydrogen_PPI" : False,
    "max_neutrino_mass_loss": NEUTRINO_MASS_LOSS_UPPER_LIMIT,
    "max_NS_mass": STATE_NS_STARMASS_UPPER_LIMIT,
    "use_interp_values": True,
    "use_profiles": True,
    "use_core_masses": True,
    "allow_spin_None" : False,
    "approx_at_he_depletion": False,
    }

# pre-defined supernova models values not changed to the default model are put
# as comments
SN_MODELS = {
    "SN_MODEL_v2_01": {
#        "mechanism": "Fryer+12-delayed",
#        "engine": "",
#        "PISN": "Hendriks+23",
#        "PISN_CO_shift": 0.0,
#        "PPI_extra_mass_loss": -20.0,
#        "ECSN": "Tauris+15",
#        "conserve_hydrogen_envelope" : False,
#        "conserve_hydrogen_PPI" : False,
#        "max_neutrino_mass_loss": NEUTRINO_MASS_LOSS_UPPER_LIMIT,
#        "max_NS_mass": STATE_NS_STARMASS_UPPER_LIMIT,
        "use_interp_values": False,
#        "use_profiles": True,
        "use_core_masses": False,
#        "allow_spin_None" : False,
#        "approx_at_he_depletion": False,
    },
    "SN_MODEL_v2_02": {
#        "mechanism": "Fryer+12-delayed",
#        "engine": "",
#        "PISN": "Hendriks+23",
#        "PISN_CO_shift": 0.0,
#        "PPI_extra_mass_loss": -20.0,
#        "ECSN": "Tauris+15",
        "conserve_hydrogen_envelope" : True,
#        "conserve_hydrogen_PPI" : False,
#        "max_neutrino_mass_loss": NEUTRINO_MASS_LOSS_UPPER_LIMIT,
#        "max_NS_mass": STATE_NS_STARMASS_UPPER_LIMIT,
        "use_interp_values": False,
#        "use_profiles": True,
        "use_core_masses": False,
#        "allow_spin_None" : False,
#        "approx_at_he_depletion": False,
    },
    "SN_MODEL_v2_03": {
#        "mechanism": "Fryer+12-delayed",
#        "engine": "",
#        "PISN": "Hendriks+23",
#        "PISN_CO_shift": 0.0,
        "PPI_extra_mass_loss": 0.0,
#        "ECSN": "Tauris+15",
#        "conserve_hydrogen_envelope" : False,
#        "conserve_hydrogen_PPI" : False,
#        "max_neutrino_mass_loss": NEUTRINO_MASS_LOSS_UPPER_LIMIT,
#        "max_NS_mass": STATE_NS_STARMASS_UPPER_LIMIT,
        "use_interp_values": False,
#        "use_profiles": True,
        "use_core_masses": False,
#        "allow_spin_None" : False,
#        "approx_at_he_depletion": False,
    },
    "SN_MODEL_v2_04": {
#        "mechanism": "Fryer+12-delayed",
#        "engine": "",
#        "PISN": "Hendriks+23",
#        "PISN_CO_shift": 0.0,
        "PPI_extra_mass_loss": 0.0,
#        "ECSN": "Tauris+15",
        "conserve_hydrogen_envelope" : True,
#        "conserve_hydrogen_PPI" : False,
#        "max_neutrino_mass_loss": NEUTRINO_MASS_LOSS_UPPER_LIMIT,
#        "max_NS_mass": STATE_NS_STARMASS_UPPER_LIMIT,
        "use_interp_values": False,
#        "use_profiles": True,
        "use_core_masses": False,
#        "allow_spin_None" : False,
#        "approx_at_he_depletion": False,
    },
    "SN_MODEL_v2_05": {
        "mechanism": "Fryer+12-rapid",
#        "engine": "",
#        "PISN": "Hendriks+23",
#        "PISN_CO_shift": 0.0,
#        "PPI_extra_mass_loss": -20.0,
#        "ECSN": "Tauris+15",
#        "conserve_hydrogen_envelope" : False,
#        "conserve_hydrogen_PPI" : False,
#        "max_neutrino_mass_loss": NEUTRINO_MASS_LOSS_UPPER_LIMIT,
#        "max_NS_mass": STATE_NS_STARMASS_UPPER_LIMIT,
        "use_interp_values": False,
#        "use_profiles": True,
        "use_core_masses": False,
#        "allow_spin_None" : False,
#        "approx_at_he_depletion": False,
    },
    "SN_MODEL_v2_06": {
        "mechanism": "Fryer+12-rapid",
#        "engine": "",
#        "PISN": "Hendriks+23",
#        "PISN_CO_shift": 0.0,
#        "PPI_extra_mass_loss": -20.0,
#        "ECSN": "Tauris+15",
        "conserve_hydrogen_envelope" : True,
#        "conserve_hydrogen_PPI" : False,
#        "max_neutrino_mass_loss": NEUTRINO_MASS_LOSS_UPPER_LIMIT,
#        "max_NS_mass": STATE_NS_STARMASS_UPPER_LIMIT,
        "use_interp_values": False,
#        "use_profiles": True,
        "use_core_masses": False,
#        "allow_spin_None" : False,
#        "approx_at_he_depletion": False,
    },
    "SN_MODEL_v2_07": {
        "mechanism": "Fryer+12-rapid",
#        "engine": "",
#        "PISN": "Hendriks+23",
#        "PISN_CO_shift": 0.0,
        "PPI_extra_mass_loss": 0.0,
#        "ECSN": "Tauris+15",
#        "conserve_hydrogen_envelope" : False,
#        "conserve_hydrogen_PPI" : False,
#        "max_neutrino_mass_loss": NEUTRINO_MASS_LOSS_UPPER_LIMIT,
#        "max_NS_mass": STATE_NS_STARMASS_UPPER_LIMIT,
        "use_interp_values": False,
#        "use_profiles": True,
        "use_core_masses": False,
#        "allow_spin_None" : False,
#        "approx_at_he_depletion": False,
    },
    "SN_MODEL_v2_08": {
        "mechanism": "Fryer+12-rapid",
#        "engine": "",
#        "PISN": "Hendriks+23",
#        "PISN_CO_shift": 0.0,
        "PPI_extra_mass_loss": 0.0,
#        "ECSN": "Tauris+15",
        "conserve_hydrogen_envelope" : True,
#        "conserve_hydrogen_PPI" : False,
#        "max_neutrino_mass_loss": NEUTRINO_MASS_LOSS_UPPER_LIMIT,
#        "max_NS_mass": STATE_NS_STARMASS_UPPER_LIMIT,
        "use_interp_values": False,
#        "use_profiles": True,
        "use_core_masses": False,
#        "allow_spin_None" : False,
#        "approx_at_he_depletion": False,
    },
    "SN_MODEL_v2_09": {
        "mechanism": "Sukhbold+16-engine",
        "engine": "N20",
#        "PISN": "Hendriks+23",
#        "PISN_CO_shift": 0.0,
#        "PPI_extra_mass_loss": -20.0,
#        "ECSN": "Tauris+15",
#        "conserve_hydrogen_envelope" : False,
#        "conserve_hydrogen_PPI" : False,
#        "max_neutrino_mass_loss": NEUTRINO_MASS_LOSS_UPPER_LIMIT,
#        "max_NS_mass": STATE_NS_STARMASS_UPPER_LIMIT,
        "use_interp_values": False,
#        "use_profiles": True,
        "use_core_masses": False,
#        "allow_spin_None" : False,
#        "approx_at_he_depletion": False,
    },
    "SN_MODEL_v2_10": {
        "mechanism": "Sukhbold+16-engine",
        "engine": "N20",
#        "PISN": "Hendriks+23",
#        "PISN_CO_shift": 0.0,
#        "PPI_extra_mass_loss": -20.0,
#        "ECSN": "Tauris+15",
        "conserve_hydrogen_envelope" : True,
#        "conserve_hydrogen_PPI" : False,
#        "max_neutrino_mass_loss": NEUTRINO_MASS_LOSS_UPPER_LIMIT,
#        "max_NS_mass": STATE_NS_STARMASS_UPPER_LIMIT,
        "use_interp_values": False,
#        "use_profiles": True,
        "use_core_masses": False,
#        "allow_spin_None" : False,
#        "approx_at_he_depletion": False,
    },
    "SN_MODEL_v2_11": {
        "mechanism": "Sukhbold+16-engine",
        "engine": "N20",
#        "PISN": "Hendriks+23",
#        "PISN_CO_shift": 0.0,
        "PPI_extra_mass_loss": 0.0,
#        "ECSN": "Tauris+15",
#        "conserve_hydrogen_envelope" : False,
#        "conserve_hydrogen_PPI" : False,
#        "max_neutrino_mass_loss": NEUTRINO_MASS_LOSS_UPPER_LIMIT,
#        "max_NS_mass": STATE_NS_STARMASS_UPPER_LIMIT,
        "use_interp_values": False,
#        "use_profiles": True,
        "use_core_masses": False,
#        "allow_spin_None" : False,
#        "approx_at_he_depletion": False,
    },
    "SN_MODEL_v2_12": {
        "mechanism": "Sukhbold+16-engine",
        "engine": "N20",
#        "PISN": "Hendriks+23",
#        "PISN_CO_shift": 0.0,
        "PPI_extra_mass_loss": 0.0,
#        "ECSN": "Tauris+15",
        "conserve_hydrogen_envelope" : True,
#        "conserve_hydrogen_PPI" : False,
#        "max_neutrino_mass_loss": NEUTRINO_MASS_LOSS_UPPER_LIMIT,
#        "max_NS_mass": STATE_NS_STARMASS_UPPER_LIMIT,
        "use_interp_values": False,
#        "use_profiles": True,
        "use_core_masses": False,
#        "allow_spin_None" : False,
#        "approx_at_he_depletion": False,
    },
    "SN_MODEL_v2_13": {
        "mechanism": "Patton&Sukhbold20-engine",
        "engine": "N20",
#        "PISN": "Hendriks+23",
#        "PISN_CO_shift": 0.0,
#        "PPI_extra_mass_loss": -20.0,
#        "ECSN": "Tauris+15",
#        "conserve_hydrogen_envelope" : False,
#        "conserve_hydrogen_PPI" : False,
#        "max_neutrino_mass_loss": NEUTRINO_MASS_LOSS_UPPER_LIMIT,
#        "max_NS_mass": STATE_NS_STARMASS_UPPER_LIMIT,
        "use_interp_values": False,
#        "use_profiles": True,
        "use_core_masses": False,
#        "allow_spin_None" : False,
#        "approx_at_he_depletion": False,
    },
    "SN_MODEL_v2_14": {
        "mechanism": "Patton&Sukhbold20-engine",
        "engine": "N20",
#        "PISN": "Hendriks+23",
#        "PISN_CO_shift": 0.0,
#        "PPI_extra_mass_loss": -20.0,
#        "ECSN": "Tauris+15",
        "conserve_hydrogen_envelope" : True,
#        "conserve_hydrogen_PPI" : False,
#        "max_neutrino_mass_loss": NEUTRINO_MASS_LOSS_UPPER_LIMIT,
#        "max_NS_mass": STATE_NS_STARMASS_UPPER_LIMIT,
        "use_interp_values": False,
#        "use_profiles": True,
        "use_core_masses": False,
#        "allow_spin_None" : False,
#        "approx_at_he_depletion": False,
    },
    "SN_MODEL_v2_15": {
        "mechanism": "Patton&Sukhbold20-engine",
        "engine": "N20",
#        "PISN": "Hendriks+23",
#        "PISN_CO_shift": 0.0,
        "PPI_extra_mass_loss": 0.0,
#        "ECSN": "Tauris+15",
#        "conserve_hydrogen_envelope" : False,
#        "conserve_hydrogen_PPI" : False,
#        "max_neutrino_mass_loss": NEUTRINO_MASS_LOSS_UPPER_LIMIT,
#        "max_NS_mass": STATE_NS_STARMASS_UPPER_LIMIT,
        "use_interp_values": False,
#        "use_profiles": True,
        "use_core_masses": False,
#        "allow_spin_None" : False,
#        "approx_at_he_depletion": False,
    },
    "SN_MODEL_v2_16": {
        "mechanism": "Patton&Sukhbold20-engine",
        "engine": "N20",
#        "PISN": "Hendriks+23",
#        "PISN_CO_shift": 0.0,
        "PPI_extra_mass_loss": 0.0,
#        "ECSN": "Tauris+15",
        "conserve_hydrogen_envelope" : True,
#        "conserve_hydrogen_PPI" : False,
#        "max_neutrino_mass_loss": NEUTRINO_MASS_LOSS_UPPER_LIMIT,
#        "max_NS_mass": STATE_NS_STARMASS_UPPER_LIMIT,
        "use_interp_values": False,
#        "use_profiles": True,
        "use_core_masses": False,
#        "allow_spin_None" : False,
#        "approx_at_he_depletion": False,
    },
    "SN_MODEL_v2_17": {
        "mechanism": "Sukhbold+16-engine",
        "engine": "W20",
#        "PISN": "Hendriks+23",
#        "PISN_CO_shift": 0.0,
#        "PPI_extra_mass_loss": -20.0,
#        "ECSN": "Tauris+15",
#        "conserve_hydrogen_envelope" : False,
#        "conserve_hydrogen_PPI" : False,
#        "max_neutrino_mass_loss": NEUTRINO_MASS_LOSS_UPPER_LIMIT,
#        "max_NS_mass": STATE_NS_STARMASS_UPPER_LIMIT,
        "use_interp_values": False,
#        "use_profiles": True,
        "use_core_masses": False,
#        "allow_spin_None" : False,
#        "approx_at_he_depletion": False,
    },
    "SN_MODEL_v2_18": {
        "mechanism": "Sukhbold+16-engine",
        "engine": "W20",
#        "PISN": "Hendriks+23",
#        "PISN_CO_shift": 0.0,
#        "PPI_extra_mass_loss": -20.0,
#        "ECSN": "Tauris+15",
        "conserve_hydrogen_envelope" : True,
#        "conserve_hydrogen_PPI" : False,
#        "max_neutrino_mass_loss": NEUTRINO_MASS_LOSS_UPPER_LIMIT,
#        "max_NS_mass": STATE_NS_STARMASS_UPPER_LIMIT,
        "use_interp_values": False,
#        "use_profiles": True,
        "use_core_masses": False,
#        "allow_spin_None" : False,
#        "approx_at_he_depletion": False,
    },
    "SN_MODEL_v2_19": {
        "mechanism": "Sukhbold+16-engine",
        "engine": "W20",
#        "PISN": "Hendriks+23",
#        "PISN_CO_shift": 0.0,
        "PPI_extra_mass_loss": 0.0,
#        "ECSN": "Tauris+15",
#        "conserve_hydrogen_envelope" : False,
#        "conserve_hydrogen_PPI" : False,
#        "max_neutrino_mass_loss": NEUTRINO_MASS_LOSS_UPPER_LIMIT,
#        "max_NS_mass": STATE_NS_STARMASS_UPPER_LIMIT,
        "use_interp_values": False,
#        "use_profiles": True,
        "use_core_masses": False,
#        "allow_spin_None" : False,
#        "approx_at_he_depletion": False,
    },
    "SN_MODEL_v2_20": {
        "mechanism": "Sukhbold+16-engine",
        "engine": "W20",
#        "PISN": "Hendriks+23",
#        "PISN_CO_shift": 0.0,
        "PPI_extra_mass_loss": 0.0,
#        "ECSN": "Tauris+15",
        "conserve_hydrogen_envelope" : True,
#        "conserve_hydrogen_PPI" : False,
#        "max_neutrino_mass_loss": NEUTRINO_MASS_LOSS_UPPER_LIMIT,
#        "max_NS_mass": STATE_NS_STARMASS_UPPER_LIMIT,
        "use_interp_values": False,
#        "use_profiles": True,
        "use_core_masses": False,
#        "allow_spin_None" : False,
#        "approx_at_he_depletion": False,
    },
    "SN_MODEL_v2_21": {
        "mechanism": "Patton&Sukhbold20-engine",
        "engine": "W20",
#        "PISN": "Hendriks+23",
#        "PISN_CO_shift": 0.0,
#        "PPI_extra_mass_loss": -20.0,
#        "ECSN": "Tauris+15",
#        "conserve_hydrogen_envelope" : False,
#        "conserve_hydrogen_PPI" : False,
#        "max_neutrino_mass_loss": NEUTRINO_MASS_LOSS_UPPER_LIMIT,
#        "max_NS_mass": STATE_NS_STARMASS_UPPER_LIMIT,
        "use_interp_values": False,
#        "use_profiles": True,
        "use_core_masses": False,
#        "allow_spin_None" : False,
#        "approx_at_he_depletion": False,
    },
    "SN_MODEL_v2_22": {
        "mechanism": "Patton&Sukhbold20-engine",
        "engine": "W20",
#        "PISN": "Hendriks+23",
#        "PISN_CO_shift": 0.0,
#        "PPI_extra_mass_loss": -20.0,
#        "ECSN": "Tauris+15",
        "conserve_hydrogen_envelope" : True,
#        "conserve_hydrogen_PPI" : False,
#        "max_neutrino_mass_loss": NEUTRINO_MASS_LOSS_UPPER_LIMIT,
#        "max_NS_mass": STATE_NS_STARMASS_UPPER_LIMIT,
        "use_interp_values": False,
#        "use_profiles": True,
        "use_core_masses": False,
#        "allow_spin_None" : False,
#        "approx_at_he_depletion": False,
    },
    "SN_MODEL_v2_23": {
        "mechanism": "Patton&Sukhbold20-engine",
        "engine": "W20",
#        "PISN": "Hendriks+23",
#        "PISN_CO_shift": 0.0,
        "PPI_extra_mass_loss": 0.0,
#        "ECSN": "Tauris+15",
#        "conserve_hydrogen_envelope" : False,
#        "conserve_hydrogen_PPI" : False,
#        "max_neutrino_mass_loss": NEUTRINO_MASS_LOSS_UPPER_LIMIT,
#        "max_NS_mass": STATE_NS_STARMASS_UPPER_LIMIT,
        "use_interp_values": False,
#        "use_profiles": True,
        "use_core_masses": False,
#        "allow_spin_None" : False,
#        "approx_at_he_depletion": False,
    },
    "SN_MODEL_v2_24": {
        "mechanism": "Patton&Sukhbold20-engine",
        "engine": "W20",
#        "PISN": "Hendriks+23",
#        "PISN_CO_shift": 0.0,
        "PPI_extra_mass_loss": 0.0,
#        "ECSN": "Tauris+15",
        "conserve_hydrogen_envelope" : True,
#        "conserve_hydrogen_PPI" : False,
#        "max_neutrino_mass_loss": NEUTRINO_MASS_LOSS_UPPER_LIMIT,
#        "max_NS_mass": STATE_NS_STARMASS_UPPER_LIMIT,
        "use_interp_values": False,
#        "use_profiles": True,
        "use_core_masses": False,
#        "allow_spin_None" : False,
#        "approx_at_he_depletion": False,
    },
}

def get_SN_MODEL(name):
    """Get predefined supernova model with all properties.

    Parameters
    ----------
    name : str
        Name of a pre-defined supernova model.

    Returns
    -------
    SN_MODEL : dict
        Dictionary with the properties of the supernova model. If the given
        name does not exist in the pre-defined models, the default supernova
        model is returned.
    """
    if name in SN_MODELS:
        SN_MODEL = DEFAULT_SN_MODEL.copy()
        SN_MODEL.update(SN_MODELS[name])
        return SN_MODEL
    else:
        return DEFAULT_SN_MODEL.copy()

def get_SN_MODEL_NAME(input_SN_MODEL, verbose=False):
    """Find SN_MODEL_NAME that matches input_SN_MODEL.

    Parameters
    ----------
    input_SN_MODEL : dict
        Dictionary with the properties of the supernova model.
    verbose : bool (default: False)
        Enables additional output.

    Returns
    -------
    SN_MODEL_NAME_SEL : str
        Name of the model in SN_MODELS that matches input_SN_MODEL.

    """
    SN_MODEL_NAME_SEL = None
    for SN_MODEL_NAME, SN_MODEL in SN_MODELS.items():
        tmp = SN_MODEL_NAME
        for key in DEFAULT_SN_MODEL.keys():
            val = SN_MODEL.get(key, DEFAULT_SN_MODEL[key])
            if "use_" in key or key=="ECSN":
                # escape values, which are allowed to differ
                continue
            if key not in input_SN_MODEL:
                if verbose:
                    print(tmp, 'missing key:', key, input_SN_MODEL)
                tmp = None
                break
            elif input_SN_MODEL[key] != val:
                if verbose:
                    print(tmp, 'mismatch:', key, input_SN_MODEL[key], val)
                tmp = None
                break
        if tmp is not None:
            if verbose:
                print('matched to supernova model:', tmp)
            SN_MODEL_NAME_SEL = tmp

    return SN_MODEL_NAME_SEL
