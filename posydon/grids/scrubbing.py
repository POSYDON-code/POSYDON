"""Methods for scrubbing MESA histories from backsteps etc."""


__authors__ = [
    "Konstantinos Kovlakas <Konstantinos.Kovlakas@unige.ch>",
]


import numpy as np
import warnings


def scrub(tables, models, ages):
    """Scrub binary and star histories of a MESA run in one go.

    Parameters
    ----------
    tables : list
        All history tables to be returned scrubbed.
    models : list
        Model number columns corresponding to each table in `tables`.
    ages : list
        Age columns corresponding to each table in `tables`.

    Returns
    -------
    list
        The list of scrubbed tables, ensuring monotonicity in age, and
        overlap in model numbers.

    """
    # scrub each history separately
    scr_tables = []
    scr_models = []
    for table, model, age in zip(tables, models, ages):
        if table is None:
            scr_tables.append(None)
            scr_models.append(None)
        else:
            n = len(model)
            userow = np.zeros(n, dtype=bool)
            userow[-1] = True
            last_m, last_t = model[-1], age[-1]
            for i in range(n-2, -1, -1):
                t_i, m_i = age[i], model[i]
                if t_i < last_t and m_i < last_m:
                    userow[i] = True
                    last_m = m_i
                    last_t = t_i
            scr_tables.append(table[userow])
            scr_models.append(model[userow])

    # find common models in scrubbed histories
    common_models = None
    for scr_model in scr_models:
        if scr_model is not None:
            if common_models is None:
                common_models = set(scr_model)
            else:
                common_models &= set(scr_model)

    # return only the rows for which models appear in all histories
    final_tables = []
    for scr_table, scr_model in zip(scr_tables, scr_models):
        if scr_table is None:
            final_tables.append(None)
            continue
        userow = np.array([m in common_models for m in scr_model], dtype=bool)
        final_tables.append(scr_table[userow])

    return final_tables


def keep_after_RLO(bh, h1, h2):
    """Scrub histories from the initial steps without RLO mass transfer.

    Parameters
    ----------
    bh : array
        The `binary_history` array.
    h1 : array
        The `history1` array.
    h2 : array
        The `history2` array.

    Returns
    -------
    tuple or arrays, or None
        The binary, history1 and history2 arrays after removing the leading
        steps without RLO overflow mass transfer from star1 or star2. If an
        RLO phase is not detected, a warning is raised, and returns `None`.

    """
    if bh is None:
        return bh, h1, h2

    bh_colnames = bh.dtype.names

    if "lg_mtransfer_rate" in bh_colnames:
        rate = bh["lg_mtransfer_rate"] >= -12
    else:
        raise ValueError("No `lg_mtransfer_rate` in binary history.")

    rlo1 = (bh["rl_relative_overflow_1"] >= -0.05
            if "rl_relative_overflow_1" in bh_colnames else None)

    rlo2 = (bh["rl_relative_overflow_2"] >= -0.05
            if "rl_relative_overflow_2" in bh_colnames else None)

    if rlo1 is None:
        rlo_1_or_2 = rlo2
    elif rlo2 is None:
        rlo_1_or_2 = rlo1
    else:
        rlo_1_or_2 = rlo1 | rlo2

    if rlo_1_or_2 is None:
        raise ValueError("No `rl_relative_overflow` in any star history.")

    conditions_met = rlo_1_or_2 & rate
    where_conditions_met = np.where(conditions_met)[0]

    if len(where_conditions_met) == 0:
        warnings.warn("No RLO overflow for this binary.")
        return None

    first_index = where_conditions_met[0]

    new_bh = None if bh is None else bh[first_index:]
    new_h1 = None if h1 is None else h1[first_index:]
    new_h2 = None if h2 is None else h2[first_index:]

    age_to_remove = new_bh["age"][0]
    new_ages = new_bh["age"] - age_to_remove

    # check if numerical precision was lost, and fix the issue
    if len(new_ages) > 1 and min(np.diff(new_ages)) == 0.0:
        min_dt = min(np.diff(new_bh["age"]))
        new_age_to_remove = (age_to_remove // min_dt) * min_dt
        relative_error = abs(new_age_to_remove - age_to_remove) / age_to_remove
        if relative_error > 0.01:
            raise Exception("Numerical precision fix too aggressive.")
        age_to_remove = new_age_to_remove
        new_ages = new_bh["age"] - age_to_remove
        if np.any(np.diff(new_ages) <= 0.0):
            raise Exception("Numerical precision fix failed.")

    new_bh["age"] = new_ages
    if new_h1 is not None and "star_age" in new_h1.dtype.names:
        new_h1["star_age"] -= age_to_remove
    if new_h2 is not None and "star_age" in new_h2.dtype.names:
        new_h2["star_age"] -= age_to_remove

    return new_bh, new_h1, new_h2
