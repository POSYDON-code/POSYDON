"""The PSY-CRIS regression module."""


__authors__ = [
    "Kyle Akira Rocha <kylerocha2024@u.northwestern.edu>",
    "Scott Coughlin <scottcoughlin2014@u.northwestern.edu>",
]


import collections
import time
from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sklearn.gaussian_process as gp

# -------- regressors --------
from scipy.interpolate import LinearNDInterpolator, Rbf
from scipy.spatial.qhull import QhullError

# from sklearn.gaussian_process import GaussianProcessRegressor
# from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C

# -----------------------------


LinearNDInterpolator_names = [
    "linear",
    "lin",
    "linearndinterpolator",
    "linear nd interpolator",
]
RBF_names = ["rbf", "radialbasisfunction", "radial basis function"]
GaussianProcessRegressor_names = ["gp", "gpr", "gaussianprocessregressor"]


def makehash():
    """Manage nested dictionaries."""
    return collections.defaultdict(makehash)


class Regressor:
    """Perform regression/interpolation with different regression algorithms.

    Regression algorithms are trained by class and by output column in the data
    set and stored as instance variables in nested dictionaries.

    This class inlcudes a 'cross validation' method that trains with the
    holdout method but calculates differences instead of a single accuracy.
    """

    def __init__(self, TableData_object):
        """Initialize the Regressor instance.

        Parameters
        ----------
        TableData_object : instance of <class, TableData>
            An instance of the TableData class.

        """
        self._TableData_ = TableData_object

        holder = self._TableData_.get_regr_data(what_data="full")
        self.input_dict = holder[0]  # _regr_inputs_
        self.output_dict = holder[1]  # _regr_outputs_
        self.regr_dfs_per_class = holder[2]  # _regr_dfs_per_class_

        max_apc_vals = []
        apc_dfs_per_class = self.regr_dfs_per_class.copy()
        for key, val in apc_dfs_per_class.items():
            non_APC_cols = [i for i in val.columns if "APC" not in i]
            if len(non_APC_cols) == len(val.columns):
                continue  # no APC cols to work on
            apc_dfs_per_class[key] = val.drop(columns=non_APC_cols)
            abs_max_val = np.nanmax(abs(apc_dfs_per_class[key].to_numpy()))
            max_apc_vals.append(abs_max_val)
        if len(max_apc_vals) != 0:
            self.abs_max_APC = np.max(max_apc_vals)
        else:
            self.abs_max_APC = None

        self._undefined_p_change_val_ = self._TableData_._return_data_(
            "undefined_p_change_val"
        )

        self._regressors_ = makehash()
        self._cv_regressors_ = makehash()

        self._log_history_ = makehash()
        self._cv_log_history = makehash()

        self.__train_cross_val = False

    def train_everything(self, regressor_names, verbose=False):
        """Train all classes and columns with the specified list of regressors.

        Parameters
        ----------
        regressor_names : list
            List of strings specifying all the regressors to train.
        verbose : optional, bool
            Print useful information.

        Returns
        -------
        None

        """
        for regr_name in regressor_names:
            if verbose:
                print("Regressor: {0}".format(regr_name))
            class_keys = list(self.regr_dfs_per_class.keys())
            for class_name in class_keys:
                self.train(regr_name, [class_name], None, verbose=verbose)
        if verbose:
            print("\nDone Regressor train_everything.")
        return None

    def train(self, regressor_name, class_keys, col_keys, di=None,
              verbose=False):
        """Train a regression algorithm.

        Implemented regressors:
            LinearNDInterpolator ('linear', ...)
            Radial Basis Function ('rbf', ...)
            GaussianProcessRegressor ('gp', ...)

        >>> rg = Regressor( TableData_object )
        >>> rg.train('linear', di = np.arange(0, Ndatapoints, 5), verbose=True)

        Trained regressor objects are uniquely defined by the algorithm used to
        train, the data set used to train (grouped by class), and finally the
        output column (there could be more than one). This motivates the data
        structure for storing the regressor objects as follows:
            Algorithm -> Class -> Output Column -> Object
        Here is more realistic example of what it could look like:
          {RBF: {"class_1": {"output_1": {instance of scipy.interpolate.rbf}}}}

        Parameters
        ----------
        regressor_name : string
            Name of regressor to train.
        class_keys : list
            List of class(es) to train on.
        col_keys : list or None
            For a given class, what columns to train on.
            If None, it trains on all columns in one class.
        di : optional, array
            Array indicies of data used to train (training on a subset).
            If None (default) - train on whole data set
        verbose : optional, bool
            Print statements with more information while training.

        Returns
        -------
        None


        Note: You can train mutliple classes at once as long as they have the
        same columns specified in col_keys.
        """
        regressor_key = self.get_regressor_name_to_key(regressor_name)

        if col_keys is None:
            first_class_data = self.regr_dfs_per_class[class_keys[0]]
            if isinstance(first_class_data, pd.DataFrame):
                col_keys = np.array(first_class_data.keys())
                if verbose:
                    print(
                        "\t Training on all {0} columns in '{1}'...".format(
                            len(col_keys), class_keys[0]
                        )
                    )
            else:
                if verbose:
                    print("No regression data for {0}.".format(class_keys[0]))
                return

        if regressor_key == "LinearNDInterpolator":
            regr_holder = self.fit_linear_ND_interpolator(
                class_keys, col_keys, data_interval=di, verbose=verbose
            )
        elif regressor_key == "RBF":
            regr_holder = self.fit_rbf_interpolator(
                class_keys, col_keys, data_interval=di, verbose=verbose
            )
        elif regressor_key == "GaussianProcessRegressor":
            regr_holder = self.fit_gaussian_process_regressor(
                class_keys, col_keys, data_interval=di, verbose=verbose
            )
        else:
            print("No trainers with name {0}".format(regressor_name))
            return

        for class_key, class_dict in regr_holder.items():
            for col_key, interpolated_obj in class_dict.items():
                if verbose:
                    print(
                        "\tdict loc: {0}, {1}, {2},".format(
                            regressor_key, class_key, col_key
                        )
                    )
                if self.__train_cross_val:
                    self._cv_regressors_[regressor_key][class_key][
                        col_key
                    ] = interpolated_obj
                else:
                    self._regressors_[regressor_key][class_key][
                        col_key
                    ] = interpolated_obj

        if verbose:
            print("\tEXIT TRAIN\n")
        return None

    def _get_cleaned_regression_data_(self, training_x, training_y,
                                      class_key, col_key):
        """Check for NaNs and user-specified `undefined_p_change_val`.

        Given a set of training data, the output is checked for nans and
        user-specified undefined_p_change_val. All instances are removed
        before training. Returns the new training input and output data:
        training_x, training_y.

        Parameters
        ----------
        training_x : ndarray
            Input data to clean.
        training_y : array
            Ouptut data to clean.
        class_key : str
            Which class is being cleaned.
        col_key : str
            Which column is being cleaned.

        Returns
        -------
        training_x : ndarray
            Cleaned input data free of undefined values.
        training_y : array
            Cleaned output data free of undefined values.

        """
        if np.sum(pd.isna(training_y)) > 0:
            where_undef = np.where(pd.isna(training_y))[0]
            where_def = np.where(pd.notna(training_y))[0]
            need_to_clean = True
        elif self._undefined_p_change_val_ in training_y:
            where_undef = np.where(self._undefined_p_change_val_
                                   == training_y)[0]
            where_def = np.where(self._undefined_p_change_val_
                                 != training_y)[0]
            need_to_clean = True
        else:
            need_to_clean = False
        if need_to_clean:
            print("Not training on {0} value(s) in {1}, {2}.".
                  format(len(where_undef), class_key, col_key))
            training_x = training_x[where_def]
            training_y = training_y[where_def]
        return training_x, training_y

    def fit_linear_ND_interpolator(self, class_keys, col_keys,
                                   data_interval=None, verbose=False):
        """Fit linear ND interpolator.

        Implementation from: scipy.interpolate.LinearNDInterpolator
        (https://docs.scipy.org/doc/scipy/reference/interpolate.html)

        Parameters
        ----------
        class_keys : list
            List of classes to train on.
        col_keys : list
            List of columns in the class to train on.
            If multiple classes are given, it is assumed they all contain
            the supplied columns.
        data_interval : array, optional
            Array indicies of data used to train (training on a subset).
            If None (default) train on whole data set
        verbose : bool, optional
            Print statements with more information while training.

        Returns
        -------
        regressor_holder : dict
            Ordered by class specific data and then by column. Nested
            dictionary maps to a trained linearNDinterpolator object.
        """
        if verbose:
            print("--- Fit LinearNDInterpolator ---")

        start_time = time.time()
        regressor_holder = OrderedDict()

        for class_key in class_keys:
            this_class_dict = OrderedDict()  # will hold columns

            # extract the output data associated with class_key
            which_class_data = self.regr_dfs_per_class[class_key]

            for col_key in col_keys:

                if data_interval is None:
                    training_x = self.input_dict[class_key].to_numpy(float)
                    training_y = which_class_data[col_key].to_numpy(float)
                else:
                    di = np.array(data_interval)
                    training_x = self.input_dict[class_key].to_numpy(float)[di]
                    training_y = which_class_data[col_key].to_numpy(float)[di]

                # if any undefined_p_change_val in regression data, remove it
                training_x, training_y = self._get_cleaned_regression_data_(
                    training_x, training_y, class_key, col_key)

                if verbose:
                    print(
                        "%s: %s - %.0f training points"
                        % (class_key, col_key, len(training_x))
                    )

                try:
                    line = LinearNDInterpolator(training_x, training_y)
                except QhullError as err:
                    if verbose:
                        print("Error: {}".format(err))
                    print("Skipping linearNDinterpolator training")
                    line = None

                this_class_dict[col_key] = line
            regressor_holder[class_key] = this_class_dict

        if verbose:
            print("--- Done in {0:.2f} seconds. ---".
                  format(time.time() - start_time))
        return regressor_holder

    def fit_rbf_interpolator(self, class_keys, col_keys, data_interval=None,
                             verbose=False):
        """Fit RBF interpolator - binary classification (one against all).

        Implementation from: scipy.interpolate.Rbf
        (https://docs.scipy.org/doc/scipy/reference/interpolate.html)

        Parameters
        ----------
        class_keys :
            List of classes to train on.
        col_keys :
            If multiple classes are given, it is assumed they all contain
            the supplied columns.
        data_interval : array, optional
            Array indicies of data used to train (training on a subset).
            if None (default) train on whole data set
        verbose : bool, optional
            Print statements with more information while training.

        Returns
        -------
        regressor_holder : dict
            Ordered by class specific data and then by column. Nested
            dictionary maps to a trained RBF object.

        """
        if verbose:
            print("--- Fit RBF ---")

        start_time = time.time()
        regressor_holder = OrderedDict()

        for class_key in class_keys:
            this_class_dict = OrderedDict()  # will hold columns

            # extract the output data associated with class_key
            which_class_data = self.regr_dfs_per_class[class_key]

            for col_key in col_keys:

                if data_interval is None:
                    training_x = self.input_dict[class_key].to_numpy(float)
                    training_y = which_class_data[col_key].to_numpy(float)
                else:
                    di = np.array(data_interval)
                    training_x = self.input_dict[class_key].to_numpy(float)[di]
                    training_y = which_class_data[col_key].to_numpy(float)[di]

                # if any undefined_p_change_val in regression data, remove it
                training_x, training_y = self._get_cleaned_regression_data_(
                    training_x, training_y, class_key, col_key)

                argList = []
                for col in range(len(training_x[0])):
                    argList.append(training_x.T[col])
                argList.append(training_y)

                if verbose:
                    print(
                        "%s: %s - %.0f training points"
                        % (class_key, col_key, len(training_x))
                    )
                if len(training_x) <= 1:
                    print("Skipping training... not enough points for Rbf")
                    # Rbf will fail for training with one point.
                    # So we put None here.
                    line = None
                else:
                    line = Rbf(*argList)

                this_class_dict[col_key] = line
            regressor_holder[class_key] = this_class_dict

        if verbose:
            print("--- Done in {0:.2f} seconds. ---".
                  format(time.time() - start_time))

        return regressor_holder

    def fit_gaussian_process_regressor(self, class_keys, col_keys,
                                       data_interval=None, verbose=False):
        """Fit a Gaussian Process regressor.

        Implementation from: sklearn.gaussian_process
        (https://scikit-learn.org/stable/modules/gaussian_process.html)

        Parameters
        ----------
        class_keys :
            List of classes to train on.
        col_keys :
            If multiple classes are given, it is assumed they all contain
            the supplied columns.
        data_interval : array, optional
            Array indicies of data used to train (training on a subset).
            if None (default) train on whole data set
        verbose : bool, optional
            Print statements with more information while training.

        Returns
        -------
        regressor_holder : dict
            Ordered by class specific data and then by column. Nested
            dictionary maps to a trained GaussianProcessRegressor object.
        """
        if verbose:
            print("--- Fit GaussianProcessRegressor ---")

        start_time = time.time()
        n_restarts = 3
        regressor_holder = OrderedDict()

        for class_key in class_keys:
            this_class_dict = OrderedDict()  # will hold columns
            # extract the output data associated with class_key
            which_class_data = self.regr_dfs_per_class[class_key]

            for col_key in col_keys:

                if data_interval is None:
                    training_x = self.input_dict[class_key].to_numpy(float)
                    training_y = which_class_data[col_key].to_numpy(float)
                else:
                    di = np.array(data_interval)
                    training_x = self.input_dict[class_key].to_numpy(float)[di]
                    training_y = which_class_data[col_key].to_numpy(float)[di]

                # if any undefined_p_change_val in regression data, remove it
                training_x, training_y = self._get_cleaned_regression_data_(
                    training_x, training_y, class_key, col_key)

                if verbose:
                    print(
                        "%s: %s - %.0f training points"
                        % (class_key, col_key, len(training_x))
                    )

                num_dim = len(training_x[0])
                starting_loc = [1 for i in range(num_dim)]
                axis_ranges = [(1e-3, 1e3) for i in range(num_dim)]
                # kernel = C( 1e3, (1e2, 5e4) ) * RBF(
                #     [10, 500, 300.], [(1e0, 1e3), (1e0, 1e3), (1e-1, 5e3)])
                kernel = gp.kernels.RBF(starting_loc, axis_ranges)
                gpr = gp.GaussianProcessRegressor(
                    kernel=kernel, n_restarts_optimizer=n_restarts
                )

                if verbose:
                    print(
                        " PRE-fit params:\n{0}".format(gpr.kernel.get_params())
                    )  # helpful for kernel things
                gpr.fit(training_x, training_y)
                if verbose:
                    print("POST-fit params:\n{0}".
                          format(gpr.kernel_.get_params()))

                this_class_dict[col_key] = gpr
            regressor_holder[class_key] = this_class_dict

        if verbose:
            print("--- Done in {0:.2f} seconds. ---".
                  format(time.time() - start_time))
        return regressor_holder

    def get_predictions(self, regressor_names, class_keys, col_keys,
                        test_input, return_std=False):
        """Get predictions from trained regressors for a set of inputs.

        Parameters
        ----------
        regressor_names : list
            List of regressor algorithm names to use to predict.
        class_keys : list
            List of classes to get predictions for.
        col_keys : list
            List of columns to get predictions for.
        test_input : ndarray
            Array of input points for which predictions will be found.
        return_std : optional, bool
            Return the STD is when using GaussianProcessRegressor.

        Returns
        -------
        predictions : dict
            Dictionary ordered by algorithm, class, and output column mapping
            to an array of predictions for the test input points.

        """
        predictions = OrderedDict()
        for regr_name in regressor_names:
            regr_key = self.get_regressor_name_to_key(regr_name)
            this_class_dict = OrderedDict()
            for class_key in class_keys:
                these_cols_dict = OrderedDict()

                for col_key in col_keys:
                    # will return None for failed Rbf, otherwise ndarray
                    pred_vals = self._predict(regr_key, class_key,
                                              col_key, test_input,
                                              return_std=return_std)
                    these_cols_dict[col_key] = pred_vals

                this_class_dict[class_key] = these_cols_dict
            predictions[regr_key] = this_class_dict

        return predictions

    def _predict(self, regressor_name, class_key, col_key, test_input,
                 return_std=False):
        """Evaluate the trained regressor at test_input and return predictions.

        If using GaussianProcessRegressor, the std is optionally returned.

        """
        if isinstance(test_input, list):
            test_input = np.array(test_input)
        if test_input.ndim == 1:
            test_input = np.array([test_input])
        if len(test_input) == 0:
            # given bad data
            return None

        sigma = None  # default

        # if empty
        if not bool(self._regressors_) and not bool(self._cv_regressors_):
            raise Exception("\n\nNo trained interpolators exist.")

        regressor_key = self.get_regressor_name_to_key(regressor_name)

        if self.__train_cross_val:
            interpolators = self._cv_regressors_[regressor_key]
        else:
            interpolators = self._regressors_[regressor_key]

        interp = interpolators[class_key][col_key]

        if regressor_key == "RBF":
            # When Rbf training fails for small classes, interpolator is None
            if interp is None:
                return None

            argList = []
            for col in range(len(test_input[0])):
                argList.append(test_input.T[col])
            pred = interp(*argList)
        elif regressor_key == "GaussianProcessRegressor":
            if return_std:
                pred, sigma = interp.predict(test_input, return_std=True)
            else:
                pred = interp.predict(test_input)
        elif regressor_key == "LinearNDInterpolator":
            pred = interp(test_input)
        else:
            print("Name not recognized: {0}".format(regressor_name))

        if return_std:
            return np.array(pred), np.array(sigma)
        else:
            return np.array(pred)

    def get_regressor_name_to_key(self, name):
        """Return the standard key (str) of a classifier."""
        if name.lower() in LinearNDInterpolator_names:
            key = "LinearNDInterpolator"
        elif name.lower() in RBF_names:
            key = "RBF"
        elif name.lower() in GaussianProcessRegressor_names:
            key = "GaussianProcessRegressor"
        else:
            print("No regressor with name '%s'." % name)
            return None
        return key

    def show_structure(self):
        """Show (print) the structure of the regression data."""
        for outer_key, outer_val in self.regr_dfs_per_class.items():
            print("CLASS: {0}".format(outer_key))
            if isinstance(outer_val, pd.DataFrame):
                print("\tCOLS:")
                for mid_key, mid_val in outer_val.items():
                    print("\t" + mid_key)
        print("")
        return None

    def get_cross_val_data(self, class_key, col_key, alpha):
        """Randomly sample the data set and seperate training and test data.

        Parameters
        ----------
        class_key : str, class_dtype(int or other)
            Class key specifying the class to get data from.
        col_key : str
            Column key specifying the output column to get data.
        alpha : float
            Fraction of data set to use for training. (0.05 = 5% of data set)

        Returns
        -------
        cross_val_test_input_data : ndarray
            Input data used to test after training on a subset.
        cross_val_test_output_data : ndarray
            Output data used to test after training on a subset.
        sorted_rnd_int_vals : array
            Indicies of original data that were used as training points.
        """
        num_points = int(len(self.input_dict[class_key]) * alpha)
        rnd_input_train = []
        rnd_outout_train = []
        rnd_int_vals = []
        rnd_int_set = set()

        # print("Num points", num_points)
        if alpha > 1 or alpha <= 0:
            raise ValueError("Alpha must be in the range (0,1].")

        ct = 0
        while len(rnd_int_vals) < num_points and ct < 1e7:
            rnd_int = int(np.random.random() * len(self.input_dict[class_key]))

            if rnd_int not in rnd_int_set:
                rnd_int_vals.append(rnd_int)
                rnd_int_set.add(rnd_int)
            ct += 1

        train_rnd_int_vals = np.array(sorted(rnd_int_vals))

        # Random training data
        # cross_val_train_input_data = (self.input_dict[class_key].
        #                               to_numpy(float))[train_rnd_int_vals, :]
        # cross_val_train_class_data = (
        #     self.regr_dfs_per_class[class_key][col_key].to_numpy(float)
        # )[train_rnd_int_vals]

        test_int_vals = []
        for i in range(len(self.input_dict[class_key])):
            if i in train_rnd_int_vals:
                pass
            else:
                test_int_vals.append(i)

        # The remainder which will be used to test fits
        cross_val_test_input_data = (self.input_dict[class_key].
                                     to_numpy(float))[test_int_vals, :]
        cross_val_test_output_data = (
            self.regr_dfs_per_class[class_key][col_key].to_numpy(float))[
                test_int_vals]

        return (cross_val_test_input_data,
                cross_val_test_output_data,
                train_rnd_int_vals)

    def cross_validate(self, regressor_name, class_key, col_key, alpha,
                       verbose=False):
        """Our method of cross validation for regression.

        Train on a subset of the data and predict values for the rest.
        Then calculate the difference between the true and predicted value.

        Parameters
        ----------
        regressor_name :
            Regressor name to use for analysis.
        class_key :
            Class key to take differences.
        col_key :
            Column key to take differences.
        alpha : float
            Fraction of data set used to find differences.
        verbose : bool, optional
            Print useful information.

        Returns
        -------
        percent_diffs : array
            Percent difference.
        diffs : array
            Absolute difference.

        """
        (
            cross_val_test_input,
            cross_val_test_output,
            train_data_indicies,
        ) = self.get_cross_val_data(class_key, col_key, alpha)

        if verbose:
            print(
                "alpha: %f, num_training_points %.0f"
                % (alpha, len(train_data_indicies))
            )

        regressor_key = self.get_regressor_name_to_key(regressor_name)

        # Train classifier
        start_time = time.time()
        try:
            self.__train_cross_val = True
            if regressor_key == "LinearNDInterpolator":
                # if linear - train rbf to use if linear predicts nan
                self.train(
                    regressor_key,
                    [class_key],
                    [col_key],
                    di=train_data_indicies,
                    verbose=verbose,
                )
                self.train(
                    "RBF",
                    [class_key],
                    [col_key],
                    di=train_data_indicies,
                    verbose=verbose,
                )
            else:
                self.train(
                    regressor_key,
                    [class_key],
                    [col_key],
                    di=train_data_indicies,
                    verbose=verbose,
                )
            time_to_train = time.time() - start_time

            # Make Predictions
            if regressor_key == "LinearNDInterpolator":
                predicted_values_linear = self._predict(
                    regressor_key, class_key, col_key, cross_val_test_input
                )
                predicted_values_rbf = self._predict(
                    "RBF", class_key, col_key, cross_val_test_input
                )
                where_nan = np.where(pd.isna(predicted_values_linear))[0]
                if len(where_nan) > 0:
                    print("{0}: {1} nan points out of {2}. Used rbf instead.".
                          format(regressor_key, len(where_nan),
                                 len(predicted_values_linear)))
                    predicted_values_linear[where_nan] = predicted_values_rbf[
                        where_nan]
                predicted_values = predicted_values_linear
            else:
                predicted_values = self._predict(regressor_key, class_key,
                                                 col_key, cross_val_test_input)
        except Exception:
            self.__train_cross_val = False
            print("FAILED DURING CROSS VAL PREDICT")
            raise
        self.__train_cross_val = False

        # Calculate the difference
        diffs = predicted_values - cross_val_test_output

        where_zero = np.where(cross_val_test_output == 0)[0]  # 1d array
        where_not_zero = np.where(cross_val_test_output != 0)[0]  # 1d array

        if len(where_zero) > 0:
            percent_diffs = (
                diffs[where_not_zero] / cross_val_test_output[where_not_zero]
            ) * 100
            print("{0} output(s) with value zero. Omitting for percent change "
                  "calculation.".format(len(where_zero)))
        else:
            percent_diffs = (diffs / cross_val_test_output) * 100

        return percent_diffs, diffs

    def get_max_APC_val(self, regressor_name, class_key, args):
        """Return the maximum interpolated average percent change for a class.

        For a given class, and regression method. Return the maximum
        interpolated average percent change value across all APC columns in
        the class sorted data set. Helper method for constructing target
        distributions for the Sampler.

        Parameters
        ----------
        regressor_name : str
            Name of regression algorithm to use.
        class_key : str
            Class key to use for data.
        args : array
            Locations for the APC value to be predicted.

        Returns
        -------
        max_APC : float
            Maximum average percent change (APC) value.
        which_col_max : int
            Index of which column had the maximum APC.

        """
        regr_column_names = self.regr_dfs_per_class[class_key].keys()
        good_col_keys = [
            i for i in regr_column_names if "APC" in i
        ]  # columns with average percent change data

        # No APC for this class
        if not good_col_keys:
            return 0, None

        regr_key = self.get_regressor_name_to_key(regressor_name)
        predictions = self.get_predictions([regr_key], [class_key],
                                           good_col_keys, args)
        dict_with_APC_data = predictions[regr_key][class_key]
        max_APC_vals = [i[0] for i in dict_with_APC_data.values()]
        max_APC = np.max(max_APC_vals)
        which_col_max = list(
            dict_with_APC_data.keys())[np.argmax(max_APC_vals)]
        return max_APC, which_col_max

    def mult_diffs(self, regressor_name, class_key, col_keys, alpha, cutoff,
                   verbose=False):
        """For multiple calls to cross_validate.

        Parameters
        ----------
        regressor_name : str
            Name of regression algorithm to use.
        class_key : str, class_dtype(int or other)
            Name of class data to use.
        col_keys : str
            Column keys to cross validate on.
        alpha : float
            Fraction of data set to cross validate on.
        cutoff : float
            Sets the cutoff percentage at which to calculate
            the fraction of the data set above or below.
        vebose : bool, optional
            Print useful diagnostic information.

        Returns
        -------
        p_diffs_holder : ndarray
            Percent differencs per column.
        attr_holder : ndarray
            Contains the number of points outside the cutoff, mean,
            and standard deviation of the percent difference calculations.

        """
        # col_keys = self.regr_dfs_per_class[class_key].keys()
        if verbose:
            print("MULT DIFFS:", regressor_name, col_keys)

        p_diffs_holder = []
        for col_key in col_keys:
            p_diffs, diffs = self.cross_validate(
                regressor_name, class_key, col_key, alpha, verbose=verbose
            )
            where_not_nan = np.where(pd.notna(p_diffs))[0]

            p_diffs_holder.append(p_diffs[where_not_nan])

        attr_holder = []
        for p_diff in p_diffs_holder:
            holder = []

            outside_cutoff = abs(p_diff) >= cutoff * 100
            num_outside = np.sum(outside_cutoff)

            holder.append(num_outside / len(p_diff) * 100)  # percent outside
            holder.append(np.mean(p_diff))  # mean
            holder.append(np.std(p_diff))  # standard deviation

            attr_holder.append(holder)

        return np.array(p_diffs_holder), np.array(attr_holder)

    def plot_regr_data(self, class_name):
        """Plot all regression data from the chosen class.

        Parameters
        ----------
        class_name : str
            Specify what class data will plotted.

        Returns
        -------
        matplotlib figure
            Plots with all regression data for a given class.

        """
        data_out = self.regr_dfs_per_class[class_name]
        data_in = self.input_dict[class_name]

        if isinstance(data_out, pd.DataFrame):
            pass
        else:
            print(
                "Output for class '{0}': {1} \nNo valid data to plot.".format(
                    class_name, str(data_out)
                )
            )
            return None

        key_in = np.array(data_in.columns)
        key_out = np.array(data_out.columns)

        # note they are still data frames until this point
        num_x_axis = len(data_in.keys())
        num_y_axis = len(data_out.keys())

        # inches per subplot - these ratios can be changed
        fig_x_ratio = 4 + 1 / 3
        fig_y_ratio = 3 + 1 / 3

        fig, subs = plt.subplots(
            nrows=num_y_axis,
            ncols=num_x_axis,
            dpi=100,
            figsize=(fig_x_ratio * num_x_axis, fig_y_ratio * num_y_axis),
        )

        # so that the indexing below works
        if num_y_axis == 1:
            subs = np.array([subs])

        print("Plotting all regression data from class '{0}'. "
              "This could take some time...".format(class_name))

        for i in range(num_x_axis):
            for k in range(num_y_axis):
                data_x = np.array(data_in[key_in[i]]).astype(float)
                data_y = np.array(data_out[key_out[k]]).astype(float)

                subs[k, i].plot(data_x, data_y, ".")
                subs[k, i].set_xlabel(key_in[i])
                subs[k, i].set_ylabel(key_out[k])
        fig.tight_layout()
        return fig

    def get_rnd_test_inputs(self, class_name, N, other_rng={}, verbose=False):
        """Produce randomly sampled 'test' inputs inside domain of input_data.

        Input data is seperated by class.

        Parameters
        ----------
        class_name : str
            Class name to specify which input data you want to look at.
        N : int
            Number of test inputs to return.
        other_rng: dict, optional
            Change the range of random sampling in desired axis. By default,
            the sampling is done in the range of the training data.
            The axis is specified with an integer key and the value
            is a list specifying the range. {1:[min, max]}
        verbose: bool, optional
            Print diagnostic information. (default False)

        Returns
        -------
        rnd_test_points : ndarray
            Test points randomly sampled in the range of the training data
            in each axis unless otherwise specified in 'other_rng'.
            Has the same shape as input data from TableData.

        """
        num_axis = len(self.input_dict[class_name].values[0])

        # find max and min in each axis
        a_max = []
        a_min = []
        for i in range(num_axis):
            a_max.append(max(self.input_dict[class_name].values.T[i]))
            a_min.append(min(self.input_dict[class_name].values.T[i]))
        # sample N points between max & min in each axis
        axis_rnd_points = []
        for i in range(num_axis):
            if i in other_rng:
                b_min, b_max = other_rng[i]
            else:
                b_min, b_max = a_min[i], a_max[i]
            if verbose:
                print("{0} - min: {1}, max: {2}".format(i, b_min, b_max))

            r = np.random.uniform(low=b_min, high=b_max, size=N)

            # this reshape is necessary to concatenate
            axis_rnd_points.append(r[:, np.newaxis])

        # now put the random points back together with same shape as input_data
        rnd_test_points = np.concatenate(axis_rnd_points, axis=1)
        return rnd_test_points
