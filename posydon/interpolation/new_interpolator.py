"""
Module implementing initial-final (IF) interpolation.

"""

__authors__ = [
    "Philipp Moura Srivastava <philipp.msrivastava@northwestern.edu>",
]

import os
import pickle
import sys
import time
from datetime import date

import numpy as np
from scipy.spatial import Delaunay
from sklearn.metrics import balanced_accuracy_score
from sklearn.model_selection import train_test_split

# ML Imports
from sklearn.neighbors import KNeighborsClassifier

# POSYDON
from posydon.grids.psygrid import PSyGrid
from posydon.interpolation.constraints import (
    find_constraints_to_apply,
    sanitize_interpolated_quantities,
)
from posydon.interpolation.data_scaling import SCALING_OPTIONS, DataScaler
from posydon.interpolation.preprocessing import (
    IN_SCALING_OPTIONS,
    OUT_SCALING_OPTIONS,
    Transformer,
    find_normalization_evaluation_matrix,
)
from posydon.utils.posydonwarning import Pwarn

eps = 1.0e-16


class IFInterpolator:
    """ Class used to train interpolator and carry out interpolation. Familiarity with the over all system, which can
    be gained by referencing section 3 of 2411.02376, is required to understand the documentation
    """

    def __init__(self, grids = None, in_keys = None, out_keys = None, max_k = None, load = False):
        """ Class constructor

            Parameters
            ----------
            grids : list of PSyGrid
                Contains both the training grid and validation grid in the first and second postions, respectively
            in_keys: list of strings
                Contains the names of the parameters that define the input space
            out_keys: dict
                Keys correspond to classifiers to be trained, and values are a list of strings specifying which parameters
                are to be interpolated using the respective classifier
            max_k: int
                The maximum number of k that is considered when optimizing k for each classifier
        """

        if type(grids) != list and not load:
            sys.exit("Please provide a list of PSyGrids containing both a training and validation grid to train the interpolator")
        elif load:
            print("Constructed in Loading Mode")
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

            # variable to control whether or not we are debugging
            self.debug_mode = False

    def stats(self, _print = False):
        """ Returns statistics regarding the number of samples that were interpolated with their initial condition outside
        of the convex hull for its respective predicted class

            Parameters
            ----------
            _print: boolean
                Controls whether statistics are printed as well as returned as a dictionary

            Return Values
            -------------
            dict: a dictionary that has the percentage of points outside the convex hull for
            each classification scheme
        """
        percentages = []

        for key in self.discrete_out_keys:
            percentages.append(
                self.outside_convex_hull[key] / (self.outside_convex_hull[key] + self.inside_convex_hull[key])
            )
        if _print:
            print(f"Total of {sum(percentages) / len(percentages):.2f} outside of hull")

        return dict(zip(self.discrete_out_keys, percentages))

    def train(self):
        """ method used to find optimal hyperparamters for classification (k) and which normalization
        schemes work best for interpolation
        """
        self.is_training = True
        self.classifiers = dict(# Finding classification hyperparameters
            zip(
                self.discrete_out_keys,
                [self.find_hyperparameters(key) for key in self.discrete_out_keys]
            )
        )
        self.out_scalers = dict(# Fincing interpolation normalization schemes
            zip(
                self.discrete_out_keys,
                [self.optimize_normalization(key) for key in self.discrete_out_keys]
            )
        )
        self.is_training = False


    def interpolate(self, iv, klass, sn_model):
        """ a method which performs interpolation for a respective initial value and its
        predicted class (convex hull)

            Parameters
            ----------
            iv: np.ndarray (3,)
                an initial value containing the primary and secondary mass as well as the orbital period
            klass: string
                a predicted class label
            sn_model: string
                the specified super nova model being used
        """

        interpolated = []
        ics = {}
        weights = {}

        interpolation_class_ind = self.discrete_out_keys.index("interpolation_class")
        sn_class_ind = self.discrete_out_keys.index(sn_model)
        klass = [klass[interpolation_class_ind], klass[sn_class_ind]]
        classification_schemes = ["interpolation_class", sn_model]

        for key, c in zip(classification_schemes, klass): # interpolating based in mass transfer type and supernova outcome separately

            triangulation = self.training_grid["triangulations"][key][c]

            simplex = -1 if triangulation == "1NN" else triangulation.find_simplex(iv)

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

            # ============= handling cases where nans exist in final values ==========
            nans = np.isnan(final_values)
            num_nans = nans.sum(axis = 0)
            not_nans = np.where((num_nans < 2) & (num_nans > 0))[0] # indices of dimensions that should not be nan despite containing nans in neighbors
            nans_tie = np.where(num_nans == 2)[0] # if there is a tie in the number of nans, we flip a coin to decide if output will be nan or not

            if self.debug_mode:
                print(final_values)

            for tie in nans_tie:
                coin = np.random.rand() # draws from uniform distribution 0-1
                if coin > 0.5:
                    np.append(not_nans, tie)

            # =============== performing general interpolation ==================
            if not self.is_training:
                label_dict = self.out_scalers[key]["label_dict"]

            # final_values = self._transform.normalize(final_values) if self.is_training else self.out_scalers[key]["transform"][label_dict[c]].normalize(final_values)
            final_values = final_values
            barycentric_weights = self.compute_barycentric_coordinates(iv, triangulation.points[vertices])[..., np.newaxis]

            weights[key] = barycentric_weights

            i_values = np.sum(final_values * barycentric_weights, axis = 0)

            # =========== fixing those values that were interpolated to nan because they contained an acceptable amount of nans in neighbors =========
            for dim in not_nans:
                where_not_nan = np.where(nans[:, dim] == 0)
                nan_weights = barycentric_weights[where_not_nan]
                nan_weights /= nan_weights.sum()

                i_values[dim] = (final_values[:, dim][where_not_nan] * nan_weights[:, 0]).sum()


            # i_values = self._transform.unnormalize(i_values[np.newaxis, ...]) if self.is_training else self.out_scalers[key]["transform"][label_dict[c]].unnormalize(i_values[np.newaxis, ...])
            i_values = [i_values]

            interpolated.extend(i_values[0])


        meta_data = {
            "weights": weights,
            "ics": ics,
            "ic": iv,
            "interpolated": interpolated
        }

        return interpolated, meta_data

    def evaluate(self, initial_values, sn_model = "S1_SN_MODEL_v2_01_SN_type"):
        """ The main method of the class used to classify and interpolate

            Parameters
            ----------
            initial_values: np.ndaarray
                starting stellar masses and orbital period in numpy array
            sn_model: string
                specifies supernova model to be used during interpolation

            Return Values
            -------------
            interpolated_values: np.ndarray
                contains all interpolated values, all classification schemes are concatenated into one array
            classes: np.ndarray
                a list of lists, each list has two classes. The first is the mass transfer type and the second
                is the compact object type which is used for the supernova
            n: list of dicts
                each dict containing meta data information about the interpolation such as
                neighbors used and distances found

        """

        if self.classifiers is None:
            sys.exit("Please find classifier hyperparameters before using interpolator")

        interpolation_class_ind = self.discrete_out_keys.index("interpolation_class")

        classes = np.array([
            cl["classifier"].predict(cl["transform"].normalize(initial_values))
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
        """ finds optimal k for a specified classifier

            Parameters
            ----------
            klass: string
                classifier to consider

            Return Values
            -------------
                dict: dict
                    contains classifier information and more
        """

        input_matrix = []
        """
        matrix that considers different number of neighbors with different
        in_scaling options
        """

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
            """ helper function handling passing of arguments for functions given
            to the preprocessing modules
            """

            kwargs = {
                "self": kwargs["kwargs"]["self"],
                "k": kwargs["item"][0],
                "scaling": kwargs["item"][1]
            }

            return kwargs

        def eval_fnc(self, k, scaling):
            """ the preprocessing module evaluates every point in input_matrix (specified above)
            which considers different input_scalings and numbers of neighbors. This function gives
            a score for each value of k paired with a normalization

                Parameters
                ----------
                k: int
                    number of neighbors used
                scaling: string
                    specifies normalization

                Return Values
                -------------
                bacc: float
                    an accuracy score
                stats: statistics used
            """

            validation_classifier = KNeighborsClassifier(n_neighbors = k, weights = "distance")

            training_initial_values = self.training_grid["initial_values"]

            transform = Transformer(training_initial_values, scaling)
            training_initial_values = transform.normalize(training_initial_values)

            validation_classifier.fit(
                training_initial_values,
                self.training_grid["final_classes"][klass]
            )

            validation_initial_values = self.validation_grid["initial_values"]
            validation_initial_values = transform.normalize(validation_initial_values)
            predicted_classes = validation_classifier.predict(validation_initial_values)

            bacc = balanced_accuracy_score(
                self.validation_grid["final_classes"][klass],
                predicted_classes
            )

            return bacc, transform

        eval_matrix, stat_matrix = find_normalization_evaluation_matrix(eval_fnc, kwargs_fnc, kwargs) # getting matrix

        k_star = list(np.unravel_index(eval_matrix.argmax(), eval_matrix.shape)) # optimal number of neighbors

        classifier = KNeighborsClassifier(n_neighbors = k_star[0] + 1, weights = "distance") # defining classifier

        training_initial_values = self.training_grid["initial_values"]

        scaling = IN_SCALING_OPTIONS[k_star[1]]

        transform = Transformer(training_initial_values, scaling)
        training_initial_values = transform.normalize(training_initial_values) # taking care of normalization

        classifier.fit(
            training_initial_values,
            self.training_grid["final_classes"][klass]
        ) # training classifier

        return {
            "classifier": classifier,
            "transform": stat_matrix[*k_star],
            "log": "log" in IN_SCALING_OPTIONS[k_star[1]],
            "k_star": k_star,
            "eval_matrix": eval_matrix
        }

    def optimize_normalization(self, key):
        """ method to find optimal normalization per class

            Parameters
            ----------
                key: string
                specifies parameter

            Return Values
            -------------
                dict: dict
                    interpolator information and more
        """

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
            """ helper function handling passing of arguments for functions given
            to the preprocessing modules
            """

            kwargs = {
                "self": kwargs["kwargs"]["self"],
                "key": kwargs["kwargs"]["key"],
                "klass": kwargs["item"][0],
                "scaling": kwargs["item"][1]
            }

            return kwargs

        def eval_fnc(self, key, klass, scaling):
            """ the preprocessing module evaluates every point in input_matrix (specified above)
            which considers different classes and output scalings. This function gives
            a score for each value of class label paired with a normalization

                Parameters
                ----------
                key: string
                    parameter considered
                klass: string
                    class label
                scaling: string
                    specifies normalization

                Return Values
                -------------
                errors: float
                    an accuracy score
                stats: statistics used
            """
            self.training = True
            self.scaling = scaling

            klass_inds = np.where(self.validation_grid["final_classes"][key] == klass)[0]

            training_final_values = self.training_grid["final_values"][key][self.training_grid["class_inds"][key][klass]]

            self._transform = Transformer(training_final_values, scaling)

            interpolated, classes, _ = self.evaluate(self.validation_grid["initial_values"][klass_inds])

            classes = classes[np.where(classes[:, 0] != "initial_MT")[0]]
            predicted_klass_inds = np.where((classes[:, 0] == klass) | (classes[:, 1] == klass))[0]

            # needs to be fixed to include any arbitrary SN model but this will do for now
            ground_truth = np.concatenate(
                [self.validation_grid["final_values"]["interpolation_class"], self.validation_grid["final_values"]["S1_SN_MODEL_v2_01_SN_type"]], axis = 1
            )

            # some interpolated and ground truth values will be NaN, e.g., for keys describing CE phenomnea are NaN if no CE is present, need to filter these NaNs out
            nans_mask = ~(np.isnan(interpolated[predicted_klass_inds]) + np.isnan(ground_truth[klass_inds][predicted_klass_inds]))

            errors = np.abs(
                (interpolated[predicted_klass_inds][nans_mask] - ground_truth[klass_inds][predicted_klass_inds][nans_mask]) /
                (ground_truth[klass_inds][predicted_klass_inds][nans_mask] + eps)
            )

            self.training = False

            error_mean = errors.mean() if len(errors) > 0 else np.inf
            return error_mean, Transformer(training_final_values, scaling)

        eval_matrix, stat_matrix = find_normalization_evaluation_matrix(eval_fnc, kwargs_fnc, kwargs) # finding normalization

        # opt = tuple(np.unravel_index(eval_matrix.argmin(axis = 0), eval_matrix.shape))
        opt = [eval_matrix.argmin(axis = 0), np.arange(eval_matrix.shape[1])]

        return {
            "transform": stat_matrix[opt[0], opt[1]],
            "eval_matrix": eval_matrix,
            "label_dict": dict(zip(labels, np.arange(len(labels))))
        }

    # =================== helper methods below ===========================

    def preprocess_grid(self, grid, training_grid = False):
        """ method that takes PSyGrid object and processes it into nice
        numpy arrays and dictionaries

            Parameters
            ---------
            grid: PSyGrid
                grid
            training_grid: bool
                specifies whether this is a training grid (regularly sampled)

            Return Values
            -------------
            dict
                PSyGrid processed into dict such that it contains only information
                important to IF interpolation
        """

        final_values = np.array(grid.final_values[self.continuous_out_keys].tolist())

        valid_inds = np.where(
            (grid.final_values["interpolation_class"] != "not_converged") &
            (grid.final_values["interpolation_class"] != "ignored_no_RLO") &
            (grid.final_values["interpolation_class"] != "ignored_no_binary_history")
        )[0]

        initial_values = np.array(grid.initial_values[self.in_keys][valid_inds].tolist())
        # determining if should interp in q
        if training_grid:
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
        """ method that constructs Delaunay triangulations stored in class memory
        when given a grid

            Parameters
            ----------
            grid_dict: dict
                a dictionary containing grid information created with preprocess grid
                method

        """

        triangulations = {}

        for label_name in self.discrete_out_keys:
            classes = np.unique(grid_dict["final_classes"][label_name]).tolist()
            if "initial_MT" in classes:
                classes.remove("initial_MT")

            class_triangulations = {}

            for klass in classes:

                class_inds = grid_dict["class_inds"][label_name][klass]

                if class_inds.shape[0] < 5:
                    print(f"too few training samples for {klass}")
                    class_triangulations[klass] = "1NN"
                else:

                    try:
                        class_triangulations[klass] = Delaunay(grid_dict["initial_values"][class_inds])
                    except:
                        print(f"Geometry wrong for {klass}, using 1NN")
                        class_triangulations[klass] = "1NN"

            triangulations[label_name] = class_triangulations

        grid_dict["triangulations"] = triangulations

    def compute_barycentric_coordinates(self, point, coords):
        """ helper method that computes barycentric coordinates which are the weights to use
        for the nearest neighbors

            Parameters
            ----------
            point: np.ndarray
                initial values inside tetrahedral (we are trying to predict)
            coords: np.ndarray
                coordinates of neighbors which are vertices of tetrahedra

            Return Values
            -------------
            np.ndarray: weight for each coord
        """

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
        """ finds the nearest neighbor in the training grid

            Parameters
            ----------
            iv: np.ndarray
                contains initial stellar masses and orbital period
            key: string
                specifies parameter

            Return Values
            -------------
                np.ndarray: output parameters of nearest neighbors
        """

        dists = np.sqrt(np.square(self.training_grid["initial_values"] - iv).sum(axis = 1))
        sorted_inds = dists.argsort()

        return np.array(self.training_grid["final_values"][key][sorted_inds[0]].tolist())

    def apply_continuous_constraints(self, interpolated, sn_model):
        """ method that applies constraints to our outputs

            Parameters
            ----------
            interpolated: np.ndarray
                interpolated values
            sn_model: string
                specifies super nova model used

            Return Values
            -------------
            np.ndarray: contains interpolated values with constraints applied

        """
        keys = self.out_key_dict["interpolation_class"] + self.out_key_dict[sn_model]

        sanitized = sanitize_interpolated_quantities(
            dict(zip(keys, interpolated)),
            self.constraints, verbose=False
        )
        return np.array([sanitized[key] for key in keys])

    def save(self, filename):
            """
            Saves the IFInterpolator instance to a pickle file.

            Parameters
            ----------
            filename : str
                Path or filename where the object should be saved (e.g., 'interpolator.pkl').
            """
            try:
                with open(filename, 'wb') as f:
                    pickle.dump(self, f, protocol=pickle.HIGHEST_PROTOCOL)
                print(f"Successfully saved interpolator to {filename}")
            except Exception as e:
                print(f"Error saving interpolator: {e}")

    def load(self, filename):
        """
        Loads an IFInterpolator instance from a pickle file.
        """
        with open(filename, 'rb') as f:
            return pickle.load(f)



