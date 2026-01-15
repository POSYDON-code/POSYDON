"""
Module implementing preprocessing for IF Interpolation
"""

__authors__ = [
    "Philipp Moura Srivastava <philipp.msrivastava@northwestern.edu>",
]

import numpy as np

eps = 1.0e-16

IN_SCALING_OPTIONS = [
    "min-max",
    "standard",
    "log_min-max",
    "log_standard"
]

OUT_SCALING_OPTIONS = [
    # "log_min-max",
    # "log_standard"
    "min-max"
]


def normalize(data, shift, scale, log = False):
    if log:
        if np.isnan(data).any() or np.isinf(data).any():
            raise ValueError("nans or infs detected in crucial grid parameters")

        if not (data < 0).any():
            data = np.log10(data + eps)

    return (data - shift) / (scale + eps)

def unnormalize(data, shift, scale, log = False):
    data = (data * (scale + eps)) + shift

    if log:
        data = 10**data - eps

    return data

def compute_statistics(data, scaling):

    computations = [
        lambda data: [data.min(axis = 0), data.max(axis = 0) - data.min(axis = 0)],
        lambda data: [data.mean(axis = 0), data.std(axis = 0)],
        lambda data: [np.log10(data + eps).min(axis = 0), np.log10(data + eps).max(axis = 0) - np.log10(data + eps).min(axis = 0)],
        lambda data: [np.log10(data + eps).mean(axis = 0), np.log10(data + eps).std(axis = 0)],
    ]
    compute = dict(zip(IN_SCALING_OPTIONS, computations)) # this line assumes that all other options are a subset of IN_SCALING_OPTION
    return compute[scaling](data)


def find_normalization_evaluation_matrix(eval_fnc, kwarg_fnc, kwargs):
    # eval_fnc - test_classifier?, what's changing
    # kwarg_fnc - input to eval_fnc
    # kwargs - inputs we iterate over

    normalization_eval_matrix = []
    normalization_stat_matrix = []
    
    for row in kwargs["input_matrix"]:
        eval_row = []
        stat_row = []

        for col in row:
            acc, stat = eval_fnc(**kwarg_fnc(**{"item": col, "kwargs": kwargs}))

            eval_row.append(
                acc
            )
            stat_row.append(stat)
        
        normalization_eval_matrix.append(eval_row)
        normalization_stat_matrix.append(stat_row)

    normalization_eval_matrix = np.array(normalization_eval_matrix)
    normalization_stat_matrix = np.array(normalization_stat_matrix)

    return normalization_eval_matrix, normalization_stat_matrix
