"""
Module implementing initial-final (IF) interpolation.

"""

import numpy as np
import os
import pickle
from datetime import date

from scipy.spatial import Delaunay
# POSYDON
from posydon.grids.psygrid import PSyGrid
from posydon.interpolation.data_scaling import DataScaler, SCALING_OPTIONS
from posydon.utils.posydonwarning import Pwarn
from posydon.interpolation.constraints import (
    find_constraints_to_apply, sanitize_interpolated_quantities)

# ML Imports
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import balanced_accuracy_score

import sys

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
            print(f"Total of {sum(pecentages) / len(percentages):.2f} outside of hul")

        return dict(zip(self.discrete_out_keys, percentages))

    def train(self):
        self.classifiers = dict(
            zip(
                self.discrete_out_keys,
                [self.find_hyperparameters(key) for key in self.discrete_out_keys]
            )
        )
        # self.out_scalers = dict(
        #     zip(
        #         self.discrete_out_keys,
        #         [self.optimize_normalization(key) for key in self.discrete_out_keys]
        #     )
        # )
        self.training = False
        
    def interpolate(self, iv, klass):

        interpolated = []
        ics = {}
        weights = {}

        for key, c in zip(self.discrete_out_keys, klass):

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

            # if self.training:
            #     final_values = self.normalize_output(final_values, c)

            barycentric_weights = self.compute_barycentric_coordinates(iv, triangulation.points[vertices])[..., np.newaxis]
            weights[key] = barycentric_weights

            # denormalized = self.denormalize_output(
            #     np.sum(final_values * barycentric_weights, axis = 0), 
            #     klass
            # )
            
            interpolated.extend(np.sum(final_values * barycentric_weights, axis = 0))

        meta_data = {
            "weights": weights, 
            "ics": ics,
            "ic": iv,
            "interpolated": interpolated
        }

        return interpolated, meta_data

    def evaluate(self, initial_values):

        if self.classifiers is None:
            sys.exit("Please find classifier hyperparameters before using interpolator")
        
        interpolation_class_ind = self.discrete_out_keys.index("interpolation_class")

        if not self.training:
            initial_values = ( np.log10(initial_values + eps) - self.iv_min ) / (self.iv_max - self.iv_min + eps)

        classes = np.array([clasifier["classifier"].predict(initial_values) for clasifier in self.classifiers.values()]).T

        interpolated_values = []
        n = []

        for iv, klass in zip(initial_values, classes):
            if klass[interpolation_class_ind] == "initial_MT":
                continue

            interpolated, meta_data = self.interpolate(iv, klass)
            interpolated = self.apply_continuous_constraints(interpolated)

            interpolated_values.append(interpolated)
            n.append(meta_data)

        interpolated_values = np.array(interpolated_values)
        classes = np.array(classes)
        print("interp_values", interpolated_values.shape)
        return interpolated_values, classes, n

    def find_hyperparameters(self, klass):

        k_accuracies = []
        
        for k in range(1, self.max_k):

            validation_classifier = KNeighborsClassifier(n_neighbors = k, weights = "distance")
            validation_classifier.fit(
                self.training_grid["initial_values"], 
                self.training_grid["final_classes"][klass]
            )

            predicted_classes = validation_classifier.predict(self.validation_grid["initial_values"])
            
            k_accuracies.append(
                balanced_accuracy_score(
                    predicted_classes, 
                    self.validation_grid["final_classes"][klass]
                )
            )

        self.k_accuracies = np.array(k_accuracies)

        k_star = self.k_accuracies.argmax()

        classifier = KNeighborsClassifier(n_neighbors = k_star, weights = "distance")
        classifier.fit(
            self.training_grid["initial_values"], 
            self.training_grid["final_classes"][klass]
        )

        return {"classifier": classifier, "k_star": k_star, "k_accuracies": k_accuracies}

    def optimize_normalization(self, key):

        output_dim = len(self.out_key_dict[key])
        classes = np.unique(self.training_grid["final_classes"][key])

        all_scalers = []
        scaler_accuracies = np.zeros(shape = (len(SCALING_OPTIONS), len(classes), output_dim))

        self.training = True
        self.scaler_types = []

        for i, option in enumerate(SCALING_OPTIONS):
            scalers = [[DataScaler() for _ in range(output_dim)] for _ in range(len(classes))]

            for j, klass in enumerate(classes):
                for parameter in range(output_dim):
                    param_min = np.nanmin(self.training_grid["final_values"][keys][:, parameter])

                    if param_min < 0 and "log" in option:
                        scalers[j][parameter].fit(self.training_grid["final_values"][keys][:, parameter], "none")
                    else:
                        scalers[j][parameter].fit(self.training_grid["final_values"][keys][:, parameter], option)

            self.out_scalers = scalers

            for j, klass in enumerate(self.classes):

                class_mask = self.validation_grid["class_inds"][klass]
                class_ivs = self.validation_grid["initial_values"][class_mask]
                
                interpolated_values, predicted_classes, _ = self.evaluate(class_ivs)
                predicted_class_mask = np.where(predicted_classes != "initial_MT") # taking out mispredictions of initial MT

                relative_errors = np.abs(
                    (interpolated_values - self.validation_grid["final_values"][keys][class_mask][predicted_class_mask]) 
                    / (self.validation_grid["final_values"][keys][class_mask][predicted_class_mask] + eps)
                ).mean(axis = 0)

                scaler_accuracies[i][j] = relative_errors

            all_scalers.append(scalers)

        all_scalers = np.array(all_scalers, dtype = object)

        star_scalers = []
        star_scaler_inds = scaler_accuracies.argmax(axis = 0)
        self.scaler_accuracies = scaler_accuracies

        for c in range(len(classes)):
            star_scalers.append(
                all_scalers[:, c][star_scaler_inds[c], np.arange(len(keys))]
            )

        self.out_scalers = np.array(star_scalers, dtype = object)
        self.normalize_triangulations()

        self.training = False
        return np.array(star_scalers, dtype = object)

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
            "initial_values": (initial_values - self.iv_min) / (self.iv_max - self.iv_min + eps),
            "final_values": dict(zip(self.out_key_dict.keys(), [grid.final_values[valid_inds][keys] for keys in self.out_key_dict.values()])), # np.array(grid.final_values[self.continuous_out_keys][valid_inds].tolist()),
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

    def apply_continuous_constraints(self, interpolated):
        sanitized = sanitize_interpolated_quantities(
            dict(zip(self.continuous_out_keys, interpolated)),
            self.constraints, verbose=False
        )
        
        return np.array([sanitized[key] for key in self.continuous_out_keys])

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
        