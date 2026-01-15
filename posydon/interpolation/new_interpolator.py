"""
Module implementing initial-final (IF) interpolation.

"""

__authors__ = [
    "Philipp Moura Srivastava <philipp.msrivastava@northwestern.edu>",
]


import numpy as np
import os
import pickle
from datetime import date

from scipy.spatial import Delaunay
# POSYDON
from posydon.grids.psygrid import PSyGrid
from posydon.interpolation.data_scaling import DataScaler, SCALING_OPTIONS
from posydon.interpolation.preprocessing import (
    normalize,
    unnormalize,
    find_normalization_evaluation_matrix, 
    compute_statistics,
    IN_SCALING_OPTIONS,
    OUT_SCALING_OPTIONS)

from posydon.utils.posydonwarning import Pwarn
from posydon.interpolation.constraints import (
    find_constraints_to_apply, sanitize_interpolated_quantities)

# ML Imports
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import balanced_accuracy_score

import sys
import time

eps = 1.0e-16


class IFInterpolator:

    def __init__(self, grids, in_keys, out_keys, max_k):

        if type(grids) != list:
            sys.exit("Please provide a list of PSyGrids containing both a training and validation grid to train the interpolator")
        else:

            self.in_keys = in_keys

            self.out_key_dict = out_keys
            self.continuous_out_keys = sum(list(out_keys.values()), []) # keys to be interpolated which correspond to numerical quantities
            self.discrete_out_keys = list(out_keys.keys()) # keys to be interpolated which correspond to discrete quantities
            self.constraints = find_constraints_to_apply(self.continuous_out_keys)

            # ============= checks =============
            if "interpolation_class" not in self.discrete_out_keys:
                sys.exit("The key \"interpolation_class\" needs to be provided as one of the interpolation keys")

            self.max_k = max_k

            self.training_grid = self.preprocess_grid(grids[0], training_grid = True)
            self.validation_grid = self.preprocess_grid(grids[1])

            self.triangulate(self.training_grid)
            # =============== usage statistics variables ============
            self.outside_convex_hull = dict(zip(self.discrete_out_keys, [0] * len(self.discrete_out_keys)))
            self.inside_convex_hull = dict(zip(self.discrete_out_keys, [0] * len(self.discrete_out_keys)))

    def stats(self, _print = False):
        percentages = []
        
        for key in self.discrete_out_keys:
            percentages.append(
                self.outside_convex_hull[key] / (self.outside_convex_hull[key] + self.inside_convex_hull[key])
            )
        if _print:
            print(f"Total of {sum(percentages) / len(percentages):.2f} outside of hull")

        return dict(zip(self.discrete_out_keys, percentages))

    def train(self):
        self.classifiers = dict(
            zip(
                self.discrete_out_keys,
                [self.find_hyperparameters(key) for key in self.discrete_out_keys]
            )
        )
        self.out_scalers = dict(
            zip(
                self.discrete_out_keys,
                [self.optimize_normalization(key) for key in self.discrete_out_keys]
            )
        )
        

    def interpolate(self, iv, klass, sn_model):

        interpolated = []
        ics = {}
        weights = {}

        interpolation_class_ind = self.discrete_out_keys.index("interpolation_class")
        sn_class_ind = self.discrete_out_keys.index(sn_model)
        klass = [klass[interpolation_class_ind], klass[sn_class_ind]]
        classification_schemes = ["interpolation_class", sn_model]

        for key, c in zip(classification_schemes, klass):

            triangulation = self.training_grid["triangulations"][key][c]

            simplex = triangulation.find_simplex(iv)            

            if simplex == -1:
                interpolated.extend(
                    self.get_nearest_neighbor(iv, key)
                )
                self.outside_convex_hull[key] += 1
                continue
            else:
                self.inside_convex_hull[key] += 1

            vertices = triangulation.simplices[simplex]
            ics[key] = triangulation.points[vertices]

            class_inds = self.training_grid["class_inds"][key][c]

            final_values = np.array(self.training_grid["final_values"][key][class_inds][vertices].tolist())
            print("Before", final_values)
            if np.isnan(final_values).any():
                print("Do final values have NaNs?", np.isnan(final_values).any())
            # if False or self.training:
            #     if not np.isnan(final_values).any():
            #         final_values = normalize(final_values, self._stats[0], self._stats[1], True)
            print("After", final_values)
            barycentric_weights = self.compute_barycentric_coordinates(iv, triangulation.points[vertices])[..., np.newaxis]
            if np.isnan(barycentric_weights).any():
                print("Do barycentric weights have NaNs?", np.isnan(barycentric_weights).any())

            weights[key] = barycentric_weights
            print(final_values, barycentric_weights)
            # if self.training:
            # if not np.isnan(final_values).any():
            #     final_values = unnormalize(
            #         np.sum(final_values * barycentric_weights, axis = 0),
            #         self._stats[0], self._stats[1], True
            #     )
            # else:
            #     final_values = np.sum(final_values * barycentric_weights, axis = 0)

            if np.isnan(final_values).any():
                print("Output has nans")
            # denormalized = self.denormalize_output(
            #     np.sum(final_values * barycentric_weights, axis = 0), 
            #     klass
            # )
            interpolated.extend(final_values)


        meta_data = {
            "weights": weights, 
            "ics": ics,
            "ic": iv,
            "interpolated": interpolated
        }

        return interpolated, meta_data

    def evaluate(self, initial_values, sn_model = "S1_SN_MODEL_v2_01_SN_type"):

        if self.classifiers is None:
            sys.exit("Please find classifier hyperparameters before using interpolator")
        
        interpolation_class_ind = self.discrete_out_keys.index("interpolation_class")

        classes = np.array([
            cl["classifier"].predict(normalize(initial_values, cl["stats"][0], cl["stats"][1], cl["log"])) 
            for cl in self.classifiers.values()]).T

        interpolated_values = []
        n = []
        
        for iv, klass in zip(initial_values, classes):
            if klass[interpolation_class_ind] == "initial_MT":
                continue

            interpolated, meta_data = self.interpolate(iv, klass, sn_model)
            interpolated = self.apply_continuous_constraints(interpolated, sn_model)

            interpolated_values.append(interpolated)
            n.append(meta_data)

        interpolated_values = np.array(interpolated_values)
        classes = np.array(classes)

        return interpolated_values, classes, n

    def find_hyperparameters(self, klass):
        
        input_matrix = []

        for k in range(1, self.max_k):
            row = []
            for opt in IN_SCALING_OPTIONS:
                row.append(
                    [k, opt]
                )
            input_matrix.append(row)

        kwargs = {
            "input_matrix": input_matrix,
            "self": self,
            "klass": klass
        }

        def kwargs_fnc(**kwargs):

            kwargs = {
                "self": kwargs["kwargs"]["self"],
                "k": kwargs["item"][0],
                "scaling": kwargs["item"][1]
            }

            return kwargs

        def eval_fnc(self, k, scaling):

            validation_classifier = KNeighborsClassifier(n_neighbors = k, weights = "distance")

            training_initial_values = self.training_grid["initial_values"]

            stats = compute_statistics(training_initial_values, scaling)
            training_initial_values = normalize(
                training_initial_values, stats[0], stats[1], "log" in scaling)

            validation_classifier.fit(
                training_initial_values, 
                self.training_grid["final_classes"][klass]
            )

            validation_initial_values = self.validation_grid["initial_values"]
            validation_initial_values = normalize(
                validation_initial_values, stats[0], stats[1], "log" in scaling)
            predicted_classes = validation_classifier.predict(validation_initial_values)

            bacc = balanced_accuracy_score(
                self.validation_grid["final_classes"][klass],
                predicted_classes                
            )

            return bacc, stats

        eval_matrix, stat_matrix = find_normalization_evaluation_matrix(eval_fnc, kwargs_fnc, kwargs)

        k_star = list(np.unravel_index(eval_matrix.argmax(), eval_matrix.shape))

        classifier = KNeighborsClassifier(n_neighbors = k_star[0] + 1, weights = "distance")

        training_initial_values = self.training_grid["initial_values"]

        scaling = IN_SCALING_OPTIONS[k_star[1]]

        stats = compute_statistics(training_initial_values, scaling)
        training_initial_values = normalize(
            training_initial_values, stats[0], stats[1], "log" in scaling)

        classifier.fit(
            training_initial_values, 
            self.training_grid["final_classes"][klass]
        )

        return {
            "classifier": classifier,
            "stats": stat_matrix[*k_star],
            "log": "log" in IN_SCALING_OPTIONS[k_star[1]],
            "k_star": k_star, 
            "eval_matrix": eval_matrix
        }

    def optimize_normalization(self, key):

        input_matrix = []

        labels = np.unique(self.training_grid["final_classes"][key])
        labels = np.delete(labels, np.where(labels == "initial_MT")[0])

        for label in labels:
            row = []
            for opt in OUT_SCALING_OPTIONS:
                row.append(
                    [label, opt]
                )
            input_matrix.append(row)
    
        kwargs = {
            "input_matrix": input_matrix,
            "self": self,
            "key": key
        }

        def kwargs_fnc(**kwargs):

            kwargs = {
                "self": kwargs["kwargs"]["self"],
                "key": kwargs["kwargs"]["key"],
                "klass": kwargs["item"][0],
                "scaling": kwargs["item"][1]
            }

            return kwargs

        def eval_fnc(self, key, klass, scaling):
            self.training = True
            self.scaling = scaling

            klass_inds = np.where(self.validation_grid["final_classes"][key] == klass)[0]

            training_final_values = self.training_grid["final_values"][key][klass_inds]
            self._stats = compute_statistics(training_final_values, scaling)
            print("Does input have NaNs?", np.isnan(self.validation_grid["initial_values"][klass_inds]).any())
            interpolated, classes, _ = self.evaluate(self.validation_grid["initial_values"][klass_inds])

            classes = classes[np.where(classes[:, 0] != "initial_MT")[0]]
            predicted_klass_inds = np.where((classes[:, 0] == klass) | (classes[:, 1] == klass))[0]

            # needs to be fixed to include any arbitrary SN model but this will do for now
            ground_truth = np.concatenate(
                [self.validation_grid["final_values"]["interpolation_class"], self.validation_grid["final_values"]["S1_SN_MODEL_v2_01_SN_type"]], axis = 1
            )

            errors = np.abs(
                (interpolated[predicted_klass_inds] - ground_truth[klass_inds][predicted_klass_inds]) /
                (ground_truth[klass_inds][predicted_klass_inds] + eps)
            )

            self.training = False

            return errors.mean(), self._stats

        eval_matrix, stat_matrix = find_normalization_evaluation_matrix(eval_fnc, kwargs_fnc, kwargs)

        opt = list(np.unravel_index(eval_matrix.argmax(), eval_matrix.shape))

        return {
            "stats": stat_matrix[opt],
            "log": "log" in OUT_SCALING_OPTIONS[opt[1]],
            "scaling": OUT_SCALING_OPTIONS[opt[1]], 
            "eval_matrix": eval_matrix
        }

    # =================== helper methods below ===========================

    def preprocess_grid(self, grid, training_grid = False):

        final_values = np.array(grid.final_values[self.continuous_out_keys].tolist())

        valid_inds = np.where(
            (grid.final_values["interpolation_class"] != "not_converged") &
            (grid.final_values["interpolation_class"] != "ignored_no_RLO") &
            (grid.final_values["interpolation_class"] != "ignored_no_binary_history")
        )[0]

        initial_values = np.array(grid.initial_values[self.in_keys][valid_inds].tolist())
        # determining if should interp in q
        if training_grid:
            m1, m2 = 10**initial_values[:, 0], 10**initial_values[:, 1]
            self.interp_in_q = (m2[m1 > 0.95 * m1.max()].min() / m2[m1 < 1.05 * m1.min()].min() > 2)
            self.interp_in_q = False
            
        initial_values = np.log10(initial_values + eps)

        if self.interp_in_q:
            initial_values[:, 1] = (10**initial_values[:, 1] - eps) / (10**initial_values[:, 0] - eps)

        if training_grid:
            self.iv_min = initial_values.min(axis = 0, keepdims = True)
            self.iv_max = initial_values.max(axis = 0, keepdims = True)

        class_inds = {}

        for key in self.discrete_out_keys:
            class_labels = np.unique(grid.final_values[valid_inds][key])
            class_inds[key] = dict(zip(
                class_labels, 
                [np.where(grid.final_values[valid_inds][key] == label)[0] for label in class_labels]
            ))

        return {
            "initial_values": 10**initial_values,
            "final_values": dict(zip(self.out_key_dict.keys(), [np.array(grid.final_values[valid_inds][keys].tolist()) for keys in self.out_key_dict.values()])), # np.array(grid.final_values[self.continuous_out_keys][valid_inds].tolist()),
            "final_classes": dict(zip(self.discrete_out_keys, np.array(grid.final_values[self.discrete_out_keys][valid_inds].tolist()).T)),
            "class_inds": class_inds,
        }
    
    def triangulate(self, grid_dict):

        triangulations = {}

        for label_name in self.discrete_out_keys:
            classes = np.unique(grid_dict["final_classes"][label_name]).tolist()
            if "initial_MT" in classes:
                classes.remove("initial_MT")

            class_triangulations = {}

            for klass in classes:

                class_inds = grid_dict["class_inds"][label_name][klass]
                class_triangulations[klass] = Delaunay(grid_dict["initial_values"][class_inds])

            triangulations[label_name] = class_triangulations

        grid_dict["triangulations"] = triangulations

    def compute_barycentric_coordinates(self, point, coords):

        T = np.array([
            coords[0] - coords[3],
            coords[1] - coords[3],
            coords[2] - coords[3]
        ]) # our matrix
        T = T.T
        T_I = np.linalg.inv(T)

        r_a = point - coords[3]

        weights = (T_I @ r_a).tolist()

        weights.append(1 - weights[0] - weights[1] - weights[2])

        weights = np.array(weights) / sum(weights)

        return weights

    def get_nearest_neighbor(self, iv, key):

        dists = np.sqrt(np.square(self.training_grid["initial_values"] - iv).sum(axis = 1))
        sorted_inds = dists.argsort()

        return np.array(self.training_grid["final_values"][key][sorted_inds[0]].tolist())

    def apply_continuous_constraints(self, interpolated, sn_model):
        keys = self.out_key_dict["interpolation_class"] + self.out_key_dict[sn_model]

        sanitized = sanitize_interpolated_quantities(
            dict(zip(keys, interpolated)),
            self.constraints, verbose=False
        )
        return np.array([sanitized[key] for key in keys])

    def normalize_output(self, input, klass):
        class_ind = self.classes.index(klass)
        scalers = self.out_scalers[class_ind]

        ret_value = np.zeros_like(input)

        for dim, scaler in enumerate(scalers):
            ret_value[:, dim] = scaler.transform(input[:, dim])

        return ret_value


    def denormalize_output(self, input, klass):
        class_ind = self.classes.index(klass)
        scalers = self.out_scalers[class_ind]

        ret_value = np.zeros_like(input)

        for dim, scaler in enumerate(scalers):
            ret_value[dim] = scaler.inv_transform(input[dim])

        return ret_value

    def normalize_triangulations(self, out_keys):

        new_final_values = np.zeros_like(self.training_grid["final_values"][out_keys])

        for i, (klass, fv) in enumerate(zip(self.training_grid["classes"], self.training_grid["final_values"][out_keys])):
            if klass == "initial_MT": # should be taken out in preprocessing
                continue
            new_final_values[i] = self.normalize_output(np.array([fv]), klass)

        self.training_grid["final_values"][out_keys] = new_final_values

        