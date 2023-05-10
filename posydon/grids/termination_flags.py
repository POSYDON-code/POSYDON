"""Module for computing termination flags of MESA runs.

There are four flags:
    Flag 1: Indicates why the run terminated.
    Flag 2: Denotes the mass transfer phase(s) the system has gone through.
    Flag 3: Marks the evolutionary stage of the primary star.
    Flag 4: Marks the evolutionary stage of the seconary star.

The funtion `get_flags_from_MESA_run` is responsible for getting all flags for
a run, using the more specialized functions defined in the module.
"""


__authors__ = [
    "Konstantinos Kovlakas <Konstantinos.Kovlakas@unige.ch>",
    "Ying Qin <<yingqin2013@hotmail.com>",
    "Devina Misra <devina.misra@unige.ch>",
    "Jeffrey Andrews <jeffrey.andrews@northwestern.edu>",
    "Emmanouil Zapartas <ezapartas@gmail.com>",
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>",
]


import os
import gzip
import warnings
import numpy as np

from posydon.utils.common_functions import (
    infer_star_state, cumulative_mass_transfer_flag, infer_mass_transfer_case,
    RL_RELATIVE_OVERFLOW_THRESHOLD, LG_MTRANSFER_RATE_THRESHOLD)
from posydon.visualization.combine_TF import (
    TF1_POOL_STABLE, TF1_POOL_UNSTABLE,
    TF1_POOL_INITIAL_RLO, TF1_POOL_ERROR, TF2_POOL_NO_RLO
)


# variables needed for inferring star states
STAR_HISTORY_VARIABLES = ["surface_h1", "center_h1", "center_he4",
                          "center_c12", "log_LH", "log_LHe", "log_Lnuc"]


def get_flag_from_MESA_output(MESA_log_path):
    """Return a flag about evolutionary outcome based on MESA output file.

    Parameters
    ----------
    MESA_log_path : string
        The output file of MESA.

    Returns
    -------
    str
        The termination flag inferred from the last line containing the strings
        `min_timestep_limit`, or `termination code:`, or `Terminate:`.

    """
    if MESA_log_path is not None and os.path.isfile(MESA_log_path):
        if MESA_log_path.endswith(".gz"):
            with gzip.open(MESA_log_path, "rt") as log_file:
                log_lines = log_file.readlines()
        else:
            with open(MESA_log_path, "r") as log_file:
                log_lines = log_file.readlines()

        for line in reversed(log_lines):
            has_term_code = line.startswith("termination code: ")
            has_terminate = line.startswith("Terminate: ")
            reached_tpagb = "Reached TPAGB" in line
            min_timestep = "min_timestep_limit" in line

            if not (has_term_code or has_terminate):
                if min_timestep and reached_tpagb:
                    return "Reached TPAGB"
            else:
                truncate_at = len(
                    "termination code: " if has_term_code else "Terminate: ")
                return line[truncate_at:].strip()

    return "reach cluster timelimit"


def get_mass_transfer_flag(binary_history, history1, history2,
                           start_at_RLO=False, mesa_flag=None):
    """Read flag from MESA history.

    In case of contact during MS, RLOF1, RLOF2 or no_RLOF.

    Parameters
    ----------
    binary_history : np.array
        The binary history file of MESA.
        Can also be the history of star 1 or 2 if binary history columns.
    history1 : np.array
        The history file of MESA for star1.
    history2 : np.array
        The history file of MESA for star2.

    Returns
    -------
    flag_system_evolution_history : string
        Possible flags are: "contact during MS", "no_RLOF", a cumulative MT
        flag + "_from_star1/2", or in case of missing binary history data,
        "undetermined_flag_mass_transfer_with_no_binary_history".

    """
    if mesa_flag is not None:
        if "overflow from L2" in mesa_flag:
            if "at ZAMS" in mesa_flag:
                return "initial_RLOF"
            else:
                return "L2_RLOF"

    if binary_history is None:
        return "undetermined_flag_mass_transfer_with_no_binary_history"

    rel1 = binary_history["rl_relative_overflow_1"]
    rel2 = binary_history["rl_relative_overflow_2"]
    rate = binary_history["lg_mtransfer_rate"]

    where_rl_rel_1 = rel1 > RL_RELATIVE_OVERFLOW_THRESHOLD
    where_rl_rel_2 = rel2 > RL_RELATIVE_OVERFLOW_THRESHOLD
    if np.any(where_rl_rel_1 & where_rl_rel_2):
        return "contact_during_MS"

    where_transfer = rate > LG_MTRANSFER_RATE_THRESHOLD
    where_rlof_1 = where_rl_rel_1 & where_transfer
    where_rlof_2 = where_rl_rel_2 & where_transfer

    if np.any(where_rlof_1):
        star = 1
        star_history = history1
        star_mass = binary_history["star_1_mass"]
        where_rlof = where_rlof_1
        rel_overflow = rel1
    elif np.any(where_rlof_2):
        star = 2
        star_history = history2
        star_mass = binary_history["star_2_mass"]
        where_rlof = where_rlof_2
        rel_overflow = rel2
    else:
        return "no_RLOF"

    indices_with_rlo = np.arange(len(where_rlof))[where_rlof]

    if not start_at_RLO and indices_with_rlo[0] == 0:
        return "initial_RLOF"

    mass_transfer_cases = []
    for index in indices_with_rlo:
        star_state = check_state_from_history(
            history=star_history, mass=star_mass, model_index=index)
        mt_case = infer_mass_transfer_case(
            rl_relative_overflow=rel_overflow[index],
            lg_mtransfer_rate=rate[index], donor_state=star_state)
        mass_transfer_cases.append(mt_case)

    flag = cumulative_mass_transfer_flag(mass_transfer_cases)
    return flag + "_from_star{}".format(star)


def check_state_from_history(history, mass, model_index=-1):
    """Return the final state of the star from the star history data.

    Parameters
    ----------
    history: np.array
        MESA history of the star
    model_index: int
        Index of the model in history for which the state will be computed.
        By default it is the end of the evolution (last model).

    Returns
    -------
    state : str
        The state of the star.

    """
    if history is None:
        return infer_star_state(star_mass=mass[model_index], star_CO=True)

    for col in STAR_HISTORY_VARIABLES:
        if col not in history.dtype.names:
            warnings.warn(
                "The data column {} is not in the history file. It is needed "
                "for the determination of the star state.".format(col))

    variables_to_pass = {varname: history[varname][model_index]
                         for varname in STAR_HISTORY_VARIABLES}
    return infer_star_state(star_CO=False, **variables_to_pass)


def get_flags_from_MESA_run(MESA_log_path, binary_history=None,
                            history1=None, history2=None, start_at_RLO=False):
    """Return the four termination flags.

    Parameters
    ----------
    MESA_log_path: str
        path to the MESA terminal output
    binary_history, history1, history2: np.array
        MESA output histories.

    Returns
    -------
    flag_MESAout, flag_mass_transfer, final_state_1, final_state_2 : str
        flag_MESAout: actual line of MESA terminal output_path
        flag_mass_transfer: describes the mass transfer (e.g., case A, case B).
        final_state_1, final_state_2: describe the final evolutionary
            state of the two stars. None if history star is not provided.
    """
    flag_out = get_flag_from_MESA_output(MESA_log_path)
    final_state_1 = check_state_from_history(history1,
                                             binary_history["star_1_mass"])
    final_state_2 = check_state_from_history(history2,
                                             binary_history["star_2_mass"])
    flag_mass_transfer = get_mass_transfer_flag(binary_history,
                                                history1, history2,
                                                start_at_RLO=start_at_RLO,
                                                mesa_flag=flag_out)
    return flag_out, flag_mass_transfer, final_state_1, final_state_2


def infer_interpolation_class(tf1, tf2, binary_history, history1,history2):
    """Use the first two termination flags to infer the interpolation class."""
    if tf1 in TF1_POOL_INITIAL_RLO:
        return "initial_MT"
    if tf1 in TF1_POOL_ERROR:
        return "not_converged"
    if tf2 in TF2_POOL_NO_RLO:
        return "no_MT"
    if tf1 in TF1_POOL_STABLE:
        if (binary_history is not None and history1 is not None):
            rl_overflow1 = binary_history['rl_relative_overflow_1']
            rl_overflow2 = binary_history['rl_relative_overflow_2']
            center_h1 = history1['center_h1']
            forward = np.logical_and(rl_overflow1>RL_RELATIVE_OVERFLOW_THRESHOLD,
                                     rl_overflow2<RL_RELATIVE_OVERFLOW_THRESHOLD)
            reverse = np.logical_and(rl_overflow1<RL_RELATIVE_OVERFLOW_THRESHOLD,
                                     rl_overflow2>RL_RELATIVE_OVERFLOW_THRESHOLD)
            reverse_postMS = np.logical_and(reverse, &
                                           center_h1<THRESHOLD_CENTRAL_ABUNDANCE)
            if (np.max(forward) == True and np.max(reverse_postMT) == True):
                return "stable_reverse_MT"
            else:
                return "stable_MT"
        return "stable_MT"
    if tf1 in TF1_POOL_UNSTABLE:
        if (binary_history is not None and history1 is not None):
            rl_overflow1 = binary_history['rl_relative_overflow_1'][-1]
            rl_overflow2 = binary_history['rl_relative_overflow_2'][-1]
            forward = np.logical_and(rl_overflow1>RL_RELATIVE_OVERFLOW_THRESHOLD,
                                     rl_overflow2<RL_RELATIVE_OVERFLOW_THRESHOLD)
            reverse = np.logical_and(rl_overflow1<RL_RELATIVE_OVERFLOW_THRESHOLD,
                                     rl_overflow2>RL_RELATIVE_OVERFLOW_THRESHOLD)
            contact = np.logical_and(rl_overflow1>RL_RELATIVE_OVERFLOW_THRESHOLD,
                                     rl_overflow2>RL_RELATIVE_OVERFLOW_THRESHOLD)
            final_state_1 = check_state_from_history(history1,
                                                   binary_history["star_1_mass"])
            final_state_2 = check_state_from_history(history2,
                                                   binary_history["star_2_mass"])
            if (forward or contact):
                if final_state_2 not in [
                        "H-rich_Core_H_burning",
                        "stripped_He_Core_He_burning", "WD", "NS", "BH"]:
                    return "DoubleCE"
                else:
                    return "unstable_MT"

            elif reverse:
                if final_state_1 not in [
                        "H-rich_Core_H_burning",
                        "stripped_He_Core_He_burning", "WD", "NS", "BH"]:
                    return "DoubleCE"
                else:
                    return "unstable_MT"
            else:
                return "unstable_MT"
        return "unstable_MT"
    return "unknown"


def get_detected_initial_RLO(grid):
    """Generates a list of already detected initial RLO
    
    Parameters
    ----------
    grid : a PSyGrid
        The grid to check.
    
    Retruns
    -------
    list
        A list containing systems already detected to be initial_MT based on
        termination_flag_1. For each system there is a dictionary with the
        impartant data, e.g. initial masses, periods, termination_flags
    """
    #new list
    detected = []
    #go through grid
    N_runs = len(grid.initial_values)
    for i in range(N_runs):
        flag1 = grid.final_values[i]["termination_flag_1"]
        #find systems with termination because of initial overflow
        if flag1 == "Terminate because of overflowing initial model":
            mass1 = grid.initial_values[i]["star_1_mass"]
            mass2 = grid.initial_values[i]["star_2_mass"]
            period = grid.initial_values[i]["period_days"]
            exists_already = False
            #check for already existing entries of same mass combination
            for j, d in enumerate(detected):
                if (abs(d["star_1_mass"]-mass1)<1.0e-5 and
                    abs(d["star_2_mass"]-mass2)<1.0e-5):
                    exists_already = True
                    #update values if new one has a larger period
                    if d["period_days"]<period:
                        detected[j]=({"star_1_mass": mass1,
                                      "star_2_mass": mass2,
                                      "period_days": period,
                                      "termination_flag_3":
                            grid.final_values[i]["termination_flag_3"],
                                      "termination_flag_4":
                            grid.final_values[i]["termination_flag_4"],
                                     })
            #add masses, period, and termination flags 3 and 4 of detected
            # system to the list
            if not exists_already:
                detected.append({"star_1_mass": mass1,
                                 "star_2_mass": mass2,
                                 "period_days": period,
                                 "termination_flag_3":
                                  grid.final_values[i]["termination_flag_3"],
                                 "termination_flag_4":
                                  grid.final_values[i]["termination_flag_4"],
                                 })
    return detected


def get_nearest_known_initial_RLO(mass1, mass2, known_initial_RLO):
    #default values
    d2min = 1.0e+99
    nearest = {"star_1_mass": 0.0,
               "star_2_mass": 0.0,
               "period_days": 0.0,
              }
    for sys in known_initial_RLO:
        #search for a known system with closest mass combination
        #use distance^2=(delta mass1)^2+(delta mass2)^2
        d2 = (mass1-sys["star_1_mass"])**2 + (mass2-sys["star_2_mass"])**2
        if d2<d2min:
            #update nearest system
            d2min = d2
            nearest = sys
    return nearest

