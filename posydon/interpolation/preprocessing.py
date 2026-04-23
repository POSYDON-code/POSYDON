"""
Module implementing preprocessing for IF Interpolation
"""

__authors__ = [
    "Philipp Moura Srivastava <philipp.msrivastava@northwestern.edu>",
]

import numpy as np

eps = 1.0e-32

IN_SCALING_OPTIONS = [
    "none",
    "min-max",
    "standard",
    "log_min-max",
    "log_standard"
]

OUT_SCALING_OPTIONS = [
    "none",
    "log_min-max",
    "log_standard",
    "min-max",
    "standard"
]

class Transformer:

    def __init__(self, data, scaling):
        """
            If a dimension contains negative values we assume that it is in log space and unlog it.
            This is an assumption that we know doesn't hold since things like rates can be negative, but it
            simplifies the preprocessing code for now.
        """
        data = data.copy()
        self.logged = (data < 0.0).any(axis = 0)
        data[:, self.logged] = 10**data[:, self.logged]

        computations = [
            lambda data: [0, 1],
            lambda data: [data.min(axis = 0), data.max(axis = 0) - data.min(axis = 0)],
            lambda data: [data.mean(axis = 0), data.std(axis = 0)],
            lambda data: [np.log10(data + eps).min(axis = 0), np.log10(data + eps).max(axis = 0) - np.log10(data + eps).min(axis = 0)],
            lambda data: [np.log10(data + eps).mean(axis = 0), np.log10(data + eps).std(axis = 0)],
        ]
        compute = dict(zip(IN_SCALING_OPTIONS, computations)) # this line assumes that all other options are a subset of IN_SCALING_OPTION

        self.log = "log" in scaling
        self.shift, self.scale = compute[scaling](data)

    def normalize(self, data):

        data = data.copy()

        if self.logged.any() and data.shape[1] == self.logged.shape[0]:
            data[:, self.logged] = 10**data[:, self.logged]

        if self.log:
            data = np.log10(data + eps)

        return (data - self.shift) / (self.scale + eps)

    # def unnormalize(self, data):

    #     data = data.copy()

    #     data = (data * (self.scale + eps)) + self.shift

    #     if self.log:
    #         data = 10**data - eps
    #         if self.logged.any() and data.shape[1] == self.logged.shape[0]:
    #             data[:, self.logged] = np.log10(data[:, self.logged])
    #     else:
    #         if self.logged.any() and data.shape[1] == self.logged.shape[0]:
    #             data[:, self.logged] = np.log10(data[:, self.logged])

    #     return data
    def unnormalize(self, data):
        data = data.copy()
        data = (data * (self.scale + eps)) + self.shift
        if self.log:
            non_logged = ~self.logged
            if non_logged.any() and data.shape[1] == self.logged.shape[0]:
                data[:, non_logged] = 10**data[:, non_logged] - eps
            elif data.shape[1] != self.logged.shape[0]:
                data = 10**data - eps  # fallback: no column info available
        if self.logged.any() and data.shape[1] == self.logged.shape[0]:
            data[:, self.logged] = np.log10(np.maximum(data[:, self.logged], eps))
        return data



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
