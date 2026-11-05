"""Tables of the supernova engines.

They are shared by the explosion criteria and the remnant-mass prescriptions,
and each table is read only once per process.
"""

import os
from functools import lru_cache

import pandas as pd

from posydon.config import PATH_TO_POSYDON_DATA
from posydon.utils.data_download import data_download

PATH_TO_SUKHBOLD16 = os.path.join(PATH_TO_POSYDON_DATA, "Sukhbold+16")
SUKHBOLD16_ENGINES = ("N20", "S19.8", "W15", "W18", "W20")
# remnant types in the tables of Sukhbold et al. (2016)
SUKHBOLD16_NS = 13
SUKHBOLD16_BH = 14


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
    (SUKHBOLD16_NS or SUKHBOLD16_BH). The table is shared: do not modify it.
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
    no_fallback = ((data["stellar_state"] == SUKHBOLD16_NS)
                   | ((data["stellar_state"] == SUKHBOLD16_BH)
                      & (data["fallback_mass"] == 0)))
    return data[no_fallback].sort_values("He_c_mass", ignore_index=True)
