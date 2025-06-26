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

# INITIAL-FINAL INTERPOLATOR
class IFInterpolator:
    """Class handling initial-final interpolation.

    Class handling initial-final interpolation with support to interpolate
    certain keys by differing classes.
    """

    def __init__(self, grids = None, interpolators = None):
        """Initialize the IFInterpolator class.

        Parameters
        ----------
        grid : PSyGrid
            The training grid or tuple of PSyGrid (training, validation)
        interpolators : list of dictionaries
            Contain parameters for the BaseIFIneterpolators to be constructed

        """
        self.interpolators = []
        self.interpolator_parameters = interpolators
        self.grids = grids

        if self.interpolator_parameters is not None:
            out_keys = []

            for params in self.interpolator_parameters:

                if ("out_keys" not in params
                        and len(self.interpolator_parameters) > 1):
                    raise ValueError("Overlapping out keys between different "
                                    "interpolators are not permited!")
                elif "out_keys" not in params:
                    continue

                for key in params["out_keys"]:
                    if(key in out_keys):
                        raise ValueError(
                            f"Overlapping out keys between different "
                            f"interpolators are not permited! ({key} in more "
                            f"than one set of out_keys)")
                    else:
                        out_keys.append(key)


    def train(self):
        """Train the interpolator(s) on the PSyGrid used for construction."""
        for interpolator in self.interpolator_parameters:
            interp = BaseIFInterpolator(grids=self.grids,
                                                         **interpolator)
            interp.train()
            self.interpolators.append(interp)

        self.interp_in_q = self.interpolators[0].interp_in_q

    def evaluate(self, binary, sanitization_verbose=False):
        """Get the output vector approximation.

        Get the output vector approximation as well as all of the
        classifications given a Binary Star

        Parameters
        ----------
        binary : BinaryStar
            A class which is an input vector which are to be classified and for
            which an output vector will be approximated.

        Returns
        -------
        Output space approximation as a tuple containing two dictionary

        """
        ynums = []
        ycats = []
        n = []

        for interpolator in self.interpolators:
            ynum, ycat, ns = interpolator.evaluate(binary)

            ynums.extend(ynum)
            ycats.extend(ycat)
            n.extend(ns)

        return ynums, ycats, n


    def test_interpolator(self, initial_values, sanitization_verbose = False):
        """ Method that can take in a 2-D numpy array for more efficient use of the interpolator

        Parameters
        ----------
        initial_values : numpy array
            A numpy array containing the in-key values of the binaries to be evolved

        Return
        ------
        Interpolated values
        """
        new_values = np.array(initial_values, copy = True)

        if self.interp_in_q:
            new_values.T[1] = new_values.T[1] / new_values.T[0]

        final_values = []

        for interpolator in self.interpolators:
            final_values.append(interpolator.test_interpolator(new_values).T)

        return np.concatenate(final_values).T

    def test_classifiers(self, initial_values):
        """ Method that can take in a 2-D numpy array for more efficient use of classifiers

        Parameters
        ----------
        initial_values : numpy array
            A numpy array containing the in-key values of the binaries to be classified

        Return
        ------
        Classified values
        """
        new_values = np.array(initial_values, copy = True)

        if self.interp_in_q:
            new_values.T[1] = new_values.T[1] / new_values.T[0]

        classes = {}

        for interpolator in self.interpolators:

            classes = {**classes, **interpolator.test_classifiers(new_values)}

        return classes



    def load(self, filename):
        """Load a saved IFInterpolator from a .pkl file.

        Parameters
        ----------
        filename : string
            Path to the .pkl file

        """
        with open(filename, 'rb') as f:
            self.interpolators = pickle.load(f)

        self.interp_in_q = self.interpolators[0].interp_in_q

    def save(self, filename):
        """Save the interpolator to a .pkl file.

        Parameters
        ----------
        filename : string
            Path where .pkl file will be saved

        """
        with open(filename, "wb") as f:
            pickle.dump(self.interpolators, f)

class BaseIFInterpolator:

    def __init__(self, grids, in_keys, out_keys, interpolation_class, max_k):

        if type(grids) != list:
            sys.exit("Please provide a list of PSyGrids containing both a training and validation grid to train the interpolator")
        else:

            self.in_keys = in_keys

            self.continous_out_keys = out_keys # keys to be interpolated which correspond to numerical quantities
            self.discrete_out_keys = out_keys # keys to be interpolated which correspond to discrete quantities
            
            self.interpolation_class = interpolation_class
            self.max_k = max_k

            self.training_grid = self.preprocess_grid(grids[0], training_grid = True)
            self.validation_grid = self.preprocess_grid(grids[1])

            self.triangulate(self.training_grid)

    def train(self):
        self.find_hyperparameters()
        # self.optimize_normalization()

    def evaluate(self, initial_values):

        if self.classifier is None:
            sys.exit("Please find classifier hyperparameters before using interpolator")

        if not self.training:
            initial_values = ( np.log10(initial_values + eps) - self.iv_min ) / (self.iv_max - self.iv_min + eps)

        classes = self.classifier.predict(initial_values)

        interpolated_values = []
        n = []

        for iv, klass in zip(initial_values, classes):
            if klass == "initial_MT":
                continue

            triangulation = self.training_grid["triangulations"][klass]

            simplex = triangulation.find_simplex(iv)

            if simplex == -1:
                interpolated_values.append(
                    self.get_nearest_neighbor(iv)
                )
                continue

            vertices = triangulation.simplices[simplex]

            class_inds = self.training_grid["class_inds"][klass]

            final_values = self.training_grid["final_values"][class_inds][vertices]

            if self.training:
                final_values = self.normalize_output(final_values, klass)

            barycentric_weights = self.compute_barycentric_coordinates(iv, triangulation.points[vertices])[..., np.newaxis]
            
            interpolated = np.sum(final_values * barycentric_weights, axis = 0)

            n.append({
                "weights": barycentric_weights, 
                "ics": triangulation.points[vertices],
                "ic": iv,
                "final_values": final_values,
                "interpolated": interpolated
            })

            # interpolated = self.denormalize_output(interpolated, klass)
            interpolated_values.append(interpolated)

        interpolated_values = np.array(interpolated_values)
        classes = np.array(classes)

        return interpolated_values, classes, n

    def find_hyperparameters(self):

        k_accuracies = []
        
        for k in range(1, self.max_k):

            validation_classifier = KNeighborsClassifier(n_neighbors = k, weights = "distance")
            validation_classifier.fit(self.training_grid["initial_values"], self.training_grid["classes"])

            predicted_classes = validation_classifier.predict(self.validation_grid["initial_values"])
            
            k_accuracies.append(
                balanced_accuracy_score(predicted_classes, self.validation_grid["classes"])
            )

        self.k_accuracies = np.array(k_accuracies)

        self.k_star = self.k_accuracies.argmax()

        self.classifier = KNeighborsClassifier(n_neighbors = self.k_star, weights = "distance")
        self.classifier.fit(self.training_grid["initial_values"], self.training_grid["classes"])

    def optimize_normalization(self):

        output_dim = len(self.continous_out_keys)

        all_scalers = []
        scaler_accuracies = np.zeros(shape = (len(SCALING_OPTIONS), len(self.classes), output_dim))

        self.training = True

        for i, option in enumerate(SCALING_OPTIONS):
            scalers = [[DataScaler()] * output_dim] * len(self.classes)

            for j, klass in enumerate(self.classes):
                for parameter in range(output_dim):
                    param_min = np.nanmin(self.training_grid["final_values"][:, parameter])
                    if param_min < 0 and "log" in option:
                        scalers[j][parameter].fit(self.training_grid["final_values"][:, parameter], "none")
                    else:
                        scalers[j][parameter].fit(self.training_grid["final_values"][:, parameter], option)
                
                self.out_scalers = scalers

                class_mask = self.validation_grid["class_inds"][klass]
                class_ivs = self.validation_grid["initial_values"][class_mask]

                interpolated_values, predicted_classes, _ = self.evaluate(class_ivs)
                predicted_class_mask = np.where(predicted_classes != "initial_MT") # taking out mispredictions of initial MT

                relative_errors = np.abs(
                    (interpolated_values - self.validation_grid["final_values"][class_mask][predicted_class_mask]) 
                    / (self.validation_grid["final_values"][class_mask][predicted_class_mask] + eps)
                ).mean(axis = 0)

                scaler_accuracies[i][j] = relative_errors

            all_scalers.append(scalers)

        all_scalers = np.array(all_scalers, dtype = object)

        star_scalers = []
        star_scaler_inds = scaler_accuracies.argmax(axis = 0)

        for c in range(len(self.classes)):
            star_scalers.append(
                all_scalers[:, c][star_scaler_inds[c], np.arange(len(self.continous_out_keys))]
            )

        self.out_scalers = np.array(star_scalers, dtype = object)
        self.normalize_triangulations()

        self.training = False

    # =================== helper methods below ===========================

    def preprocess_grid(self, grid, training_grid = False):

        final_values = np.array(grid.final_values[self.continous_out_keys].tolist())

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

        class_labels = np.unique(grid.final_values[valid_inds][self.interpolation_class])
        class_inds = dict(zip(
            class_labels, 
            [np.where(grid.final_values[valid_inds][self.interpolation_class] == label)[0] for label in class_labels]
        ))

        return {
            "initial_values": (initial_values - self.iv_min) / (self.iv_max - self.iv_min + eps),
            "final_values": np.array(grid.final_values[self.continous_out_keys][valid_inds].tolist()),
            "final_classes": grid.final_values[valid_inds][self.discrete_out_keys],
            "classes": grid.final_values[valid_inds][self.interpolation_class],
            "class_inds": class_inds,
        }
    
    def triangulate(self, grid_dict):

        self.classes = np.unique(grid_dict["classes"]).tolist()
        self.classes.remove("initial_MT")

        triangulations = {}

        for klass in self.classes:

            class_inds = grid_dict["class_inds"][klass]
            triangulations[klass] = Delaunay(grid_dict["initial_values"][class_inds])

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

    def get_nearest_neighbor(self, iv):

        dists = np.sqrt(np.square(self.training_grid["initial_values"] - iv).sum(axis = 1))
        sorted_inds = dists.argsort()

        return self.training_grid["final_values"][sorted_inds[0]]

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

    def normalize_triangulations(self):

        new_final_values = np.zeros_like(self.training_grid["final_values"])

        for i, (klass, fv) in enumerate(zip(self.training_grid["classes"], self.training_grid["final_values"])):
            if klass == "initial_MT": # should be taken out in preprocessing
                continue
            new_final_values[i] = self.normalize_output(np.array([fv]), klass)

        self.training_grid["final_values"] = new_final_values
        