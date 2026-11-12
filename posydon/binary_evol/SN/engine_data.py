"""Tables of the supernova engines.

They are shared by the explosion criteria and the remnant-mass prescriptions,
and each table is read only once per process.
"""

import os
import json
from functools import lru_cache

import numpy as np
import pandas as pd

from posydon.config import PATH_TO_POSYDON_DATA
from posydon.utils.data_download import data_download
from posydon.utils.interpolators import interp1d

PATH_TO_SUKHBOLD16 = os.path.join(PATH_TO_POSYDON_DATA, "Sukhbold+16")
SUKHBOLD16_ENGINES = ("N20", "S19.8", "W15", "W18", "W20")

PATH_TO_COUCH20 = os.path.join(PATH_TO_POSYDON_DATA, "Couch+2020")
COUCH20_TURBULENCE = ("1.0", "1.2", "1.23", "1.25", "1.27", "1.3", "1.4")

# code for remnant types in the tables of Sukhbold et al. (2016)
NS_STATE = 13
BH_STATE = 14


def _auxiliary_file(filename):
    """Return filename, downloading the auxiliary data first if it is missing."""
    if not os.path.exists(filename):
        data_download(set_name='auxiliary', confirm=True,
                      missing_file=filename)
    return filename


@lru_cache(maxsize=None)
def Sukhbold16_table(engine):
    """Return the NS and BH models of Sukhbold et al. (2016) for an engine.

    BHs formed by fallback are left out, so every BH in the table comes from
    a failed explosion. The rows are sorted by He core mass ('He_c_mass');
    'Rem_mass' is the remnant mass and 'stellar_state' the remnant type
    (NS_STATE or BH_STATE). The table is shared: do not modify it.
    """
    if engine not in SUKHBOLD16_ENGINES:
        raise ValueError(
            f"Invalid engine {engine!r} for the mechanism "
            f"'Sukhbold+16-engine'. Valid engines are: "
            f"{', '.join(SUKHBOLD16_ENGINES)}."
        )
    filename = os.path.join(PATH_TO_SUKHBOLD16,
                            f"results_{engine}_table.csv")
    data = pd.read_csv(_auxiliary_file(filename))
    no_fallback = ((data["stellar_state"] == NS_STATE)
                   | ((data["stellar_state"] == BH_STATE)
                      & (data["fallback_mass"] == 0)))
    return data[no_fallback].sort_values("He_c_mass", ignore_index=True)


@lru_cache(maxsize=None)
def Couch20_table(turbulence):
    """Return the models of Couch et al. (2020) for a turbulence strength.

    A model without explosion energy ('Eexp' = 0) forms a BH, the others a
    NS. The Couch+20 data has no He core masses; they come from the table of
    Sukhbold et al. (2016) by ZAMS mass ('Mzams'), as in v2.4. Columns:
    'Mzams', 'Eexp', 'He_c_mass' and 'stellar_state' (NS_STATE or BH_STATE).

    The rows keep the order of the data file: some models share a He core
    mass (rounded in Sukhbold's table) but not the outcome, and a
    one-neighbour classifier then takes the first of them, as in v2.4.
    The table is shared: do not modify it.
    """
    if turbulence not in COUCH20_TURBULENCE:
        raise ValueError(
            f"Invalid turbulence strength {turbulence!r} (option 'engine') "
            f"for the mechanism 'Couch+20-engine'. Valid values are: "
            f"{', '.join(COUCH20_TURBULENCE)}."
        )
    filename = os.path.join(PATH_TO_COUCH20, "explDatsSTIR2.json")
    with open(_auxiliary_file(filename)) as file:
        models = json.load(file)[turbulence]
    data = pd.DataFrame({"Mzams": [float(m) for m in models],
                         "Eexp": [model["Eexp"] for model in models.values()]})
    filename = os.path.join(PATH_TO_COUCH20, "Sukhbold_Mzams_He_c_core.csv")
    Sukhbold = pd.read_csv(_auxiliary_file(filename),
                           usecols=["Mzams", "He_c_mass"])
    data["He_c_mass"] = interp1d(Sukhbold["Mzams"],
                                 Sukhbold["He_c_mass"])(data["Mzams"])
    data["stellar_state"] = np.where(data["Eexp"] == 0.0, BH_STATE, NS_STATE)
    return data
