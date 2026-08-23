"""Common functions to be used while running populations."""


__authors__ = [
    "Konstantinos Kovlakas <Konstantinos.Kovlakas@unige.ch>",
    "Devina Misra <devina.misra@unige.ch>",
    "Emmanouil Zapartas <ezapartas@gmail.com>",
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Nam Tran <tranhn03@gmail.com>",
    "Ying Qin <<yingqin2013@hotmail.com>",
    "Jeffrey Andrews <jeffrey.andrews@northwestern.edu>",
    "Tassos Fragos <Anastasios.Fragkos@unige.ch>",
    "Scott Coughlin <scottcoughlin2014@u.northwestern.edu>",
    "Kyle Akira Rocha <kylerocha2024@u.northwestern.edu>",
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>",
    "Camille Liotine <cliotine@u.northwestern.edu>",
]


import copy
import os

import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.optimize import newton

from posydon.binary_evol.flow_chart import STAR_STATES_H_RICH, STAR_STATES_HE_RICH
from posydon.utils import constants as const
from posydon.utils.interpolators import interp1d
from posydon.utils.limits_thresholds import (
    LG_MTRANSFER_RATE_THRESHOLD,
    LOG10_BURNING_THRESHOLD,
    REL_LOG10_BURNING_THRESHOLD,
    RL_RELATIVE_OVERFLOW_THRESHOLD,
    STATE_NS_STARMASS_LOWER_LIMIT,
    STATE_NS_STARMASS_UPPER_LIMIT,
    STATE_WD_STARMASS_UPPER_LIMIT,
    THRESHOLD_CENTRAL_ABUNDANCE,
    THRESHOLD_HE_NAKED_ABUNDANCE,
)
from posydon.utils.posydonwarning import Pwarn

# Constants related to inferring star states
STATE_UNDETERMINED = "undetermined_evolutionary_state"

# ALL POSSIBLE STAR STATES
BURNING_STATES = ["Core_H_burning", "Core_He_burning",
                  "Shell_H_burning", "Core_He_depleted",
                  "Core_C_depleted"]
RICHNESS_STATES = ["H-rich", "stripped_He", "accreted_He"]
COMPACT_OBJECTS = ["WD", "NS", "BH","massless_remnant"]

ALL_STAR_STATES = COMPACT_OBJECTS + [STATE_UNDETERMINED]
ALL_STAR_STATES.extend(["{}_{}".format(rich_in, burning)
                        for rich_in in RICHNESS_STATES
                        for burning in BURNING_STATES])

# Mass-transfer cases in form of integer flags
MT_CASE_NO_RLO = 0
MT_CASE_A = 1
MT_CASE_B = 2
MT_CASE_C = 3
MT_CASE_BA = 4
MT_CASE_BB = 5
MT_CASE_BC = 6
MT_CASE_NONBURNING = 8
MT_CASE_UNDETERMINED = 9

# All cases meaning RLO is happening
ALL_RLO_CASES = set([MT_CASE_A, MT_CASE_B, MT_CASE_C,
                     MT_CASE_BA, MT_CASE_BB, MT_CASE_BC,
                     MT_CASE_NONBURNING])

# Conversion of integer mass-transfer flags to strings
MT_CASE_TO_STR = {
    MT_CASE_NO_RLO: "no_RLO",
    MT_CASE_A: "A",
    MT_CASE_B: "B",
    MT_CASE_C: "C",
    MT_CASE_BA: "BA",
    MT_CASE_BB: "BB",
    MT_CASE_BC: "BC",
    MT_CASE_NONBURNING: "nonburning",
    MT_CASE_UNDETERMINED: "undetermined_MT"
}

# Conversion of strings to integer mass-transfer flags
MT_STR_TO_CASE = {string: integer for integer, string
                  in MT_CASE_TO_STR.items()}

DEFAULT_CE_OPTION_FOR_LAMBDA = \
    "lambda_from_profile_gravitational_plus_internal_minus_recombination"

###################################################################

def rotate(axis, angle):
    """Generate rotation matrix to rotate a vector about an arbitrary axis by
        a given angle

    Parameters
    ----------
    axis : array of length 3
        Axis to rotate about
    angle : float
        Angle, in radians, through which to rotate about axis

    Returns
    -------
    rotation_matrix : 3x3 array
        Array such that rotation_matrix.dot(vector) rotates vector
        about the given axis by the given angle
    """

    if len(axis)!=3:
        raise ValueError("axis should be of dimension 3")
    # normalize the axis vector
    norm = np.linalg.norm(axis)
    if norm==0:
        raise ValueError("axis is a point")
    axis = axis / norm

    # calculate the cosine and sine of the angle
    cos_theta = np.cos(angle)
    sin_theta = np.sin(angle)

    # construct the rotation matrix
    rotation_matrix = np.array([
            [cos_theta + axis[0]**2 * (1 - cos_theta),
            axis[0] * axis[1] * (1 - cos_theta) - axis[2] * sin_theta,
            axis[0] * axis[2] * (1 - cos_theta) + axis[1] * sin_theta],
            [axis[1] * axis[0] * (1 - cos_theta) + axis[2] * sin_theta,
            cos_theta + axis[1]**2 * (1 - cos_theta),
            axis[1] * axis[2] * (1 - cos_theta) - axis[0] * sin_theta],
            [axis[2] * axis[0] * (1 - cos_theta) - axis[1] * sin_theta,
            axis[2] * axis[1] * (1 - cos_theta) + axis[0] * sin_theta,
            cos_theta + axis[2]**2 * (1 - cos_theta)]
        ])

    return rotation_matrix

def linear_interpolation_between_two_cells(array_y, array_x, x_target,
                                           top=None, bot=None, verbose=False):
    """Interpolate quantities between two star profile shells."""
    if (pd.isna(top) and pd.isna(bot)):
        top = np.argmax(array_x >= x_target)
        bot = top - 1
    elif pd.isna(bot):
        bot = top - 1
    elif pd.isna(top):
        top = bot + 1

    if top >= len(array_y):
        Pwarn("top={} is too large, use last element in array_y".format(top),
              "ReplaceValueWarning")
        top = len(array_y)-1
    if top >= len(array_x):
        Pwarn("array_x too short, use y at top={}".format(top),
              "InterpolationWarning")
        return array_y[top]
    if bot < 0:
        Pwarn("bot={} is too small, use first element".format(bot),
              "ReplaceValueWarning")
        bot = 0

    if top == bot:
        y_target = array_y[top]
        Pwarn("linear interpolation occured between the same point: x_target,"
              " top, bot, len(array_x), y_target = {}, {}, {}, {}, {}".format(\
              x_target, top, bot, len(array_x), y_target),
              "InterpolationWarning")
    elif bot > top:
        y_target = array_y[top]
        Pwarn("bot={} is too large: use y at top={}".format(bot, top),
              "InterpolationWarning")
    else:
        x_top = array_x[top]
        x_bot = array_x[bot]

        y_top = array_y[top]
        y_bot = array_y[bot]

        slope = (y_top - y_bot) / (x_top - x_bot)
        const = (y_top*x_bot - y_bot*x_top) / (x_top - x_bot)
        y_target = slope * x_target - const

        if verbose:
            print("linear interpolation")
            print("x_target, top, bot, len(array_x)",
                  x_target, top, bot, len(array_x))
            print("x_top, x_bot, y_top, y_bot, y_target",
                  x_top, x_bot, y_top, y_bot, y_target)

    return y_target

def is_number(s):
    """Check if the input can be converted to a float."""
    try:
        float(s)
        return True
    except ValueError:
        return False

def zero_negative_values(arr, key): # pragma no cover
    """
        Set negative values in the array to zero.

        Parameters
        ----------
        arr : np.ndarray
            The input array to process.
        key : string
            The name of the array column

        Returns
        -------
        np.ndarray
            The processed array with negative values set to zero.
    """
    arr = np.atleast_1d(arr)

    if np.any(arr < 0):
        Pwarn("A " + key + " value is negative. Setting to zero.",
              "ReplaceValueWarning")

    arr[arr < 0] = 0.0
    return arr

##################################
####### I/O  #####################
##################################

def initialize_empty_array(arr):
    """Initialize an empty record array with NaNs and empty strings."""
    res = arr.copy()
    for colname in res.dtype.names:
        if np.issubsctype(res[colname], float):
            res[colname] = np.nan
        if np.issubsctype(res[colname], str):
            res[colname] = np.nan
        #TODO: handle h5py.string_dtype()
    return res

def read_histogram_from_file(path):
    """Read a histogram from a CSV file.

    The expected format is:

    # comment line
    x[0], x[1], ..., x[n], x[n+1]
    y[0], y[1], ..., y[n]

    where # denotes a comment line (also empty lines are ignored), and `n` is
    the number of bins (notice that the first line contains n+1 elements.)

    Usage: bin_edges, bin_counts = read_histogram_from_file("a_histogram.csv").


    Parameters
    ----------
    path : str
        The path of the CSV file containing the histogram information.

    Returns
    -------
    list of arrays
        The bin edges and bin counts of the histogram.

    """
    with open(path, "r") as f:
        arrays = []
        for line in f:
            line = line.strip()
            if len(line) == 0 or line.startswith("#"):
                continue
            arrays.append(np.fromstring(line.strip(), dtype=float, sep=","))
            if len(arrays) > 2:
                raise IndexError("More than two lines found in the histogram"
                                 " document.")
    if len(arrays) < 2:
        raise IndexError("Less than two lines found in the histogram"
                         " document.")
    if len(arrays[0]) - 1 != len(arrays[1]):
        raise IndexError("The number of elements in the second data line is"
                         " not one less than the number in the first data"
                         " line.")

    return arrays