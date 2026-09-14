"""
Module implementing preprocessing for IF Interpolation

This module provides the normalization utilities used by ``IFInterpolator``
during both training and inference. It centers on two pieces:

1. **``Transformer``**: a reusable normalizer/denormalizer for a matrix of
   input or output quantities. Given a scaling scheme — one of ``"none"``,
   ``"min-max"``, ``"standard"``, ``"log_min-max"``, or ``"log_standard"``
   (see ``IN_SCALING_OPTIONS`` / ``OUT_SCALING_OPTIONS``) — a ``Transformer``
   computes the shift/scale statistics needed to normalize data into that
   space and back again via ``normalize`` and ``unnormalize``. A single
   scaling string can be supplied for all columns, or a list/array of
   per-column scalings. Columns whose names contain ``"log"`` or ``"lg"``
   are assumed to already be stored in log space and are un-logged before
   any further scaling is applied; columns containing negative values are
   assumed to be unsuitable for log scaling and are silently downgraded
   from a ``"log_*"`` scheme to its non-log counterpart, since a negative
   quantity (e.g. a rate) cannot be logged safely.

2. **``find_normalization_evaluation_matrix``**: a small grid-search
   helper used by ``IFInterpolator`` during training. Given an evaluation
   function, an argument-building function, and a 2D ``input_matrix`` of
   parameter combinations to try (e.g. every pairing of neighbor count
   ``k`` with each input scaling option, or every pairing of class label
   with each output scaling option), it evaluates every cell and returns
   two matrices of the same shape: one of scores (e.g. balanced accuracy
   or interpolation error) and one of the fitted objects (e.g. classifiers
   or ``Transformer`` instances) produced along the way. ``IFInterpolator``
   uses the resulting score matrix to pick the best-performing combination
   for each classifier and interpolation target.

Together these let ``IFInterpolator`` automatically search over
normalization schemes rather than requiring them to be hand-tuned per
quantity.
"""

__authors__ = [
    "Philipp Moura Srivastava <philipp.msrivastava@northwestern.edu>",
]

import sys

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

    def __init__(self, data, scaling, keys):
        """
            If a dimension contains negative values we assume that it is in log space and unlog it.
            This is an assumption that we know doesn't hold since things like rates can be negative, but it
            simplifies the preprocessing code for now.
        """
        data = data.copy()

        self.keys = np.array(keys)
        self.logged = np.array(["log" in k or "lg" in k for k in keys])
        data[:, self.logged] = 10**data[:, self.logged]

        self.negative = (data < 0).any(axis = 0)

        computations = [
            lambda data: [0, 1],
            lambda data: [data.min(axis = 0), data.max(axis = 0) - data.min(axis = 0)],
            lambda data: [data.mean(axis = 0), data.std(axis = 0)],
            lambda data: [np.log10(data + eps).min(axis = 0), np.log10(data + eps).max(axis = 0) - np.log10(data + eps).min(axis = 0)],
            lambda data: [np.log10(data + eps).mean(axis = 0), np.log10(data + eps).std(axis = 0)],
        ]
        compute = dict(zip(IN_SCALING_OPTIONS, computations)) # this line assumes that all other options are a subset of IN_SCALING_OPTION
        og = scaling[:]

        # if there's any in the data which are less than 0, it cannot be logged and a linear scaling is applied instead
        if self.negative.any() and type(scaling) != list and "log" in scaling:
            list_scaling = np.array([scaling] * len(self.keys))
            list_scaling[self.negative] = scaling.replace("log_", "")
            scaling = list_scaling
        elif self.negative.any() and type(scaling) == list or type(scaling) == np.ndarray:
            for i, s in enumerate(scaling):
                if self.negative[i] and "log" in s:
                    scaling[i] = scaling[i].replace("log_", "")

        self.scaling = scaling

        if type(scaling) == list or type(scaling) == np.ndarray:
            self.log = ["log" in s for s in scaling]


            self.shift = np.array([compute[s](data[:, i])[0] for i, s in enumerate(scaling)])
            self.scale = np.array([compute[s](data[:, i])[1] for i, s in enumerate(scaling)])
        else:
            self.log = "log" in scaling
            self.shift, self.scale = compute[scaling](data)

    def normalize(self, data):

        data = data.copy()

        if self.logged.any() and data.shape[1] == self.logged.shape[0]:
            data[:, self.logged] = 10**data[:, self.logged]


        if type(self.log) == bool and self.log:
            data = np.log10(data + eps)


        elif type(self.log) == list:
            for i, l in enumerate(self.log):

                if l:
                    data[:, i] = np.log10(data[:, i] + eps)



        data = (data - self.shift) / (self.scale + eps)

        return data


    def unnormalize(self, data):
        data = data.copy()


        data = (data * (self.scale + eps)) + self.shift

        if type(self.log) == bool and self.log:
            non_logged = ~self.logged

            if non_logged.any() and data.shape[1] == self.logged.shape[0]:
                data[:, non_logged] = 10**data[:, non_logged] - eps

            elif data.shape[1] != self.logged.shape[0]:
                data = 10**data - eps  # fallback: no column info available
        elif type(self.log) == list:
            for i, l in enumerate(self.log):
                if l:
                    data[:, i] = 10**data[:, i] - eps


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
