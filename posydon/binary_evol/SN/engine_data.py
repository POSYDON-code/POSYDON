"""Tables of the supernova engines, and the tools that place stars on them.

They are shared by the explosion criteria and the remnant-mass prescriptions.
Each table is read only once per process; objects fitted to the tables are
created by the stages that use them.
"""

import json
import os
from functools import lru_cache

import numpy as np
import pandas as pd
from sklearn.neighbors import KNeighborsRegressor

from posydon.config import PATH_TO_POSYDON_DATA
from posydon.utils.common_functions import calculate_Patton20_values_at_He_depl
from posydon.utils.data_download import data_download
from posydon.utils.interpolators import interp1d
from posydon.utils.posydonerror import ModelError

PATH_TO_SUKHBOLD16 = os.path.join(PATH_TO_POSYDON_DATA, "Sukhbold+16")
PATH_TO_PATTON20 = os.path.join(PATH_TO_POSYDON_DATA, "Patton+Sukhbold20")
PATH_TO_COUCH20 = os.path.join(PATH_TO_POSYDON_DATA, "Couch+2020")

SUKHBOLD16_ENGINES = ("N20", "S19.8", "W15", "W18", "W20")
COUCH20_TURBULENCE = ("1.0", "1.2", "1.23", "1.25", "1.27", "1.3", "1.4")
PATTON20_QUANTITIES = ("M4", "mu4", "Xi", "sc")

# code for remnant types in the tables of Sukhbold et al. (2016)
NS_STATE = 13
BH_STATE = 14

N_NEIGHBORS = 1

def _auxiliary_file(filename):
    """Return filename, downloading the auxiliary data first if it is missing."""
    if not os.path.exists(filename):
        data_download(set_name='auxiliary', confirm=True,
                      missing_file=filename)
    return filename


@lru_cache(maxsize=None)
def Sukhbold16_table(engine):
    """
    Return the NS and BH models of Sukhbold et al. (2016) for an engine.

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
def Patton20_table(quantity):
    """
    Return the pre-computed properties of pre-SN models of Patton & Sukhbold
    (2020), MNRAS, 499, 2803.

    Patton & Sukhbold evolved a grid of CO cores (without envelopes) to core
    collapse using the KEPLER stellar evolution code. Their grid covers CO
    core masses in the range of 2.5 - 10.0 solar masses, and with carbon mass
    fractions at core He depletion between 0.05 and 0.5. For each model, the
    values of M4, mu4, xi_2.5 and central entropy sc are given, so that the
    explodability criterion of Ertl et al. (2016) can be applied.
    """
    if quantity not in PATTON20_QUANTITIES:
        raise ValueError(
            f"Invalid quantity {quantity!r} of the Patton & Sukhbold (2020) "
            f"models. Valid quantities are: {', '.join(PATTON20_QUANTITIES)}."
        )
    filename = os.path.join(PATH_TO_PATTON20, f"Kepler_{quantity}_table.dat")
    data = pd.read_csv(_auxiliary_file(filename), comment="#", sep=r"\s+",
                       index_col=0)
    table = data.stack().reset_index()
    table.columns = ["X_C", "M_CO", quantity]
    table["M_CO"] = table["M_CO"].astype(float)
    return table


def _missing(value):
    """True for a value that was not set (None or NaN)."""
    return value is None or (isinstance(value, float) and np.isnan(value))


def CO_core_at_He_depletion(star, approximate=False):
    """
    Return the carbon mass fraction and CO core mass of a star at core He
    depletion, the properties that place it among the CO-core models of
    Patton & Sukhbold (2020).

    Without approximation, they come from the star's
    avg_c_in_c_core_at_He_depletion and co_core_mass_at_He_depletion, or
    from its history if those are missing. With approximate=True, the CO core
    mass at the end of the star's evolution stands in for the one at He
    depletion, and the carbon mass fraction comes from fits to it, for H-rich
    and stripped stars.

    Returns
    -------
    X_C, M_CO : float
    """
    if approximate:
        M_CO = star.co_core_mass
        if "H-rich" in star.state:
            return 0.20 / M_CO + 0.15, M_CO
        if "stripped_He" in star.state:
            return -0.084 * np.log(M_CO) + 0.4, M_CO
        raise ValueError(f"Cannot approximate the CO core at He depletion "
                         f"of a star in state {star.state!r}: it must be "
                         f"H-rich or stripped_He.")
    if (_missing(star.avg_c_in_c_core_at_He_depletion)
            or _missing(star.co_core_mass_at_He_depletion)):
        calculate_Patton20_values_at_He_depl(star)
    X_C = star.avg_c_in_c_core_at_He_depletion
    M_CO = star.co_core_mass_at_He_depletion
    if _missing(X_C) or _missing(M_CO):
        raise ModelError(f"The star has no CO core at He depletion, neither "
                         f"as attributes nor in its history: M_CO = {M_CO}, "
                         f"X_C = {X_C}.")
    return X_C, M_CO


class Patton20_CO_core_models:
    """
    Pre-supernova quantities of the CO-core models of Patton & Sukhbold
    (2020), for a CO core at He depletion.

    When created, fits one k-nearest-neighbours regressor per quantity on the
    carbon mass fraction and CO core mass (X_C, M_CO) of the models. Calling
    the object returns the quantities of the model(s) matching a CO core: with
    n_neighbors = 1, those of the closest model (outside the grid, of the
    closest model on its edge). X_C and M_CO are not rescaled, which only
    matters for n_neighbors > 1.
    """
    def __init__(self, quantities):
        self.n_neighbors = N_NEIGHBORS
        self._maps = {}
        for quantity in quantities:
            table = Patton20_table(quantity)
            regressor = KNeighborsRegressor(self.n_neighbors, weights="distance")
            self._maps[quantity] = regressor.fit(
                table[["X_C", "M_CO"]].to_numpy(), table[quantity].to_numpy())

    def __call__(self, X_C, M_CO):
        """Return {quantity: value} for a CO core at He depletion."""
        point = [[X_C, M_CO]]
        return {quantity: float(regressor.predict(point)[0])
                for quantity, regressor in self._maps.items()}

    def __repr__(self):
        return (f"{type(self).__name__}(quantities={tuple(self._maps)}, "
                f"n_neighbors={self.n_neighbors})")


@lru_cache(maxsize=None)
def Couch20_table(turbulence):
    """
    Return the models of Couch et al. (2020) for a turbulence strength.

    A model without explosion energy ('Eexp' = 0) forms a BH, the others a
    NS. The Couch+20 data has no He core masses; they come from the table of
    Sukhbold et al. (2016) by ZAMS mass ('Mzams'). Columns:
    'Mzams', 'Eexp', 'He_c_mass' and 'stellar_state' (NS_STATE or BH_STATE).

    The rows keep the order of the data file: some models share a He core
    mass (rounded in Sukhbold's table) but not the outcome, and a
    one-neighbour classifier then takes the first of them, so the order
    matters.
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
