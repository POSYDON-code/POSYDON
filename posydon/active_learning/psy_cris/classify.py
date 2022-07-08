"""The PSY-CRIS classification module."""


__authors__ = [
    "Kyle Akira Rocha <kylerocha2024@u.northwestern.edu>",
    "Scott Coughlin <scottcoughlin2014@u.northwestern.edu>",
]


import collections
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import time
from collections import OrderedDict

# -------- classifiers --------
from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import Rbf
import sklearn.gaussian_process as gp

# -----------------------------

LinearNDInterpolator_names = [
    "linear",
    "lin",
    "linearndinterpolator",
    "linear nd interpolator",
]

RBF_names = ["rbf", "radialbasisfunction", "radial basis function"]
GaussianProcessClassifier_names = ["gp", "gpc", "gaussianprocessclassifier"]


def makehash():
    """Manage nested dictionaries."""
    return collections.defaultdict(makehash)


class Classifier:
    """Classifier class.

    Perform one against all classification with a variety of different
    classification algorithms (interpolators). Different classifcation
    algorithms are trainined and stored for recall as instance variables
    inside nested dictionaries.
    This class also supports model validation through cross validation
    using the holdout method.

    """

    def __init__(self, TableData_object):
        """Initialize the classifier.

        Parameters
        ----------
        TableData_object : instance of <class, TableData>
            An instance of the TableData class with training data.

        """
        self._TableData_ = TableData_object

        holder = self._TableData_.get_class_data(what_data="full")
        self.class_names = holder[1]  # _unique_class_keys_
        self.classes_to_ids = holder[2]  # _class_col_to_ids_
        self.class_id_mapping = holder[3]  # _class_id_mapping_
        self.binary_class_data = holder[4]  # _binary_data_
        self.input_data = self._TableData_.get_data(what_data="input")

        self._interpolators_ = makehash()
        self._cv_interpolators_ = makehash()

        self.__train_cross_val = False

    def train_everything(self, classifier_names, verbose=False):
        """Trains multiple classifiers at once.

        Parameters
        ----------
        classifier_names : list
            List of strings specifying classification algorithms to use.

        Returns
        -------
        None
        """
        for cls_name in classifier_names:
            self.train(cls_name, di=None, verbose=verbose)
        return None

    def train(self, classifier_name, di=None, verbose=False, **kwargs):
        """Train a classifier.

        Implemented classifiers:
            LinearNDInterpolator ('linear', ...)
            Radial Basis Function ('rbf', ...)
            GaussianProcessClassifier ('gp', ...)

        >>> cl = Classifier( TableData_object )
        >>> cl.train('linear', di = np.arange(0, Ndatapoints, 5), verbose=True)

        Parameters
        ----------
        classifier_name : str
            Name of classifier to train.
        di : array_int, optional
            Array indicies of data used to train (training on a subset).
            if None - train on whole data set
        train_cross_val : bool, optional
            For storing regular trained interpolators and cross val
            interpolators. Used in the cross_validate() method.
            if False - save normal interpolators
            if True  - save cross validation interpolators
        verbose: bool, optional
            Print statements with more information while training.

        Returns
        -------
        None
        """
        classifier_key = self.get_classifier_name_to_key(classifier_name)

        if classifier_key == "LinearNDInterpolator":
            bi_cls_holder = self.fit_linear_ND_interpolator(data_interval=di,
                                                            verbose=verbose)
        elif classifier_key in "RBF":
            bi_cls_holder = self.fit_rbf_interpolator(data_interval=di,
                                                      verbose=verbose)
        elif classifier_key in "GaussianProcessClassifier":
            bi_cls_holder = self.fit_gaussian_process_classifier(
                data_interval=di, verbose=verbose, **kwargs)
        else:
            print("No classifiers with name {0}.".format(classifier_name))
            return

        for col_key, cls_obj in bi_cls_holder.items():
            if verbose:
                print("\tdict loc: {0} {1}".format(classifier_key, col_key))
            if self.__train_cross_val:
                self._cv_interpolators_[classifier_key][col_key] = cls_obj
            else:
                self._interpolators_[classifier_key][col_key] = cls_obj

        if verbose:
            print("Done training {0}.".format(classifier_key))
        return None

    def fit_linear_ND_interpolator(self, data_interval=None, verbose=False):
        """Fit linear ND interpolator - binary one-against-all classification.

        Implementation from: scipy.interpolate.LinearNDInterpolator
        (https://docs.scipy.org/doc/scipy/reference/interpolate.html)

        Parameters
        ----------
        data_interval : array_int, optional
            Array indicies of data used to train (training on a subset).
            if None (default) train on whole data set
        verbose : bool, optional
            Print statements with more information while training.

        Returns
        -------
        binary_classifier_holder : dict
            Sorted by class, each key maps to a trained linearNDinterpolator
            object.
        """
        if verbose:
            if data_interval is None:
                print("N points: %i" % (len(self.input_data)))
            else:
                print("N points: %i" % (len(data_interval)))

        start_time = time.time()

        binary_classifier_holder = OrderedDict()

        for i, cls_data in enumerate(self.binary_class_data):
            # iter_time = time.time()

            # for running with a subset of the data
            if data_interval is None:
                line = LinearNDInterpolator(self.input_data, cls_data)
            else:
                di = np.array(data_interval)
                line = LinearNDInterpolator(self.input_data[di], cls_data[di])

            binary_classifier_holder[self.class_names[i]] = line

            time_print = time.time() - start_time
            if verbose:
                if i == 0:
                    len_classes = len(self._TableData_.class_ids)
                    print(
                        "Time to fit %s classifiers ~ %.3f\n"
                        % (len_classes, time_print * len_classes)
                    )
                print(
                    "LinearNDInterpolator class %s -- current time: %.3f"
                    % (i, time_print)
                )

        return binary_classifier_holder

    def fit_rbf_interpolator(self, data_interval=None, verbose=False):
        """Fit RBF interpolator - binary classification (one against all).

        Implementation from: scipy.interpolate.Rbf
        (https://docs.scipy.org/doc/scipy/reference/interpolate.html)

        Parameters
        ----------
        data_interval : array_int, optional
            Array indicies of data used to train (training on a subset).
            if None (default) train on whole data set
        verbose : bool, optional
            Print statements with more information while training.

        Returns
        -------
        binary_classifier_holder : dict
            Sorted by class, each key maps to a trained RBF object.
        """
        if verbose:
            if data_interval is None:
                print("N points: %i" % (len(self.input_data)))
            else:
                print("N points: %i" % (len(data_interval)))

        start_time = time.time()

        binary_classifier_holder = OrderedDict()

        for i, cls_data in enumerate(self.binary_class_data):
            # iter_time = time.time()

            # for running with a subset of the data
            if data_interval is None:
                argList = []
                for col in range(len(self.input_data[0])):
                    argList.append(self.input_data.T[col])
                argList.append(cls_data)

                line = Rbf(*argList)
            else:
                di = np.array(data_interval)
                argList = []
                for col in range(len(self.input_data[0])):
                    argList.append(self.input_data.T[col][di])
                argList.append(cls_data[di])

                line = Rbf(*argList)

            binary_classifier_holder[self.class_names[i]] = line

            time_print = time.time() - start_time
            if verbose:
                if i == 0:
                    len_classes = len(self._TableData_.class_ids)
                    print(
                        "Time to fit %s classifiers ~ %.3f\n"
                        % (len_classes, time_print * len_classes)
                    )
                print("RBF class %s -- current time: %.3f" % (i, time_print))

        return binary_classifier_holder

    def fit_gaussian_process_classifier(self, data_interval=None,
                                        my_kernel=None, n_restarts=5,
                                        verbose=False):
        """Fit a Gaussian Process classifier.

        Implementation from: sklearn.gaussian_process
        (https://scikit-learn.org/stable/modules/gaussian_process.html)

        Parameters
        ----------
        data_interval : array_int, optional
            Array indicies of data used to train (training on a subset).
            if None (default) train on whole data set
        my_kernel : kernel
            Set the kernel for the GPC.
        n_restarts : int
            Number of restarts for the GPC.
        verbose : bool, optional
            Print statements with more information while training.

        Returns
        -------
        binary_classifier_holder : array_like
            Sorted by class, each key maps to a trained
            GaussianProcessClassifier object.
        """
        if verbose:
            if data_interval is None:
                print("N points: %i" % (len(self.input_data)))
            else:
                print("N points: %i" % (len(data_interval)))

        start_time = time.time()

        binary_classifier_holder = OrderedDict()

        for i, cls_data in enumerate(self.binary_class_data):
            # iter_time = time.time()

            if my_kernel is None:
                num_dim = len(self.input_data[0])
                starting_loc = [1 for i in range(num_dim)]
                axis_ranges = [(1e-3, 1e3) for i in range(num_dim)]
                # kernel = gp.kernels.RBF([1,1,1], [(1e-3,1e3), (1e-3,1e3),
                #                                   (1e-3, 1e3)])
                kernel = gp.kernels.RBF(starting_loc, axis_ranges)
            else:
                kernel = my_kernel

            gpc = gp.GaussianProcessClassifier(kernel=kernel,
                                               n_restarts_optimizer=n_restarts)

            # for running with a subset of the data
            if data_interval is None:
                line = gpc.fit(self.input_data, cls_data)
            else:
                di = np.array(data_interval)
                line = gpc.fit(self.input_data[di], cls_data[di])

            if verbose:
                print("\t kernel:\n{0}".format(kernel))

            binary_classifier_holder[self.class_names[i]] = line

            time_print = time.time() - start_time

            if verbose:
                if i == 0:
                    len_classes = len(self._TableData_.class_ids)
                    print(
                        "Time to fit %s classifiers ~ %.3f\n"
                        % (len_classes, time_print * len_classes)
                    )
                print(
                    "GaussianProcessClassifier class %s -- current time: %.3f"
                    % (i, time_print)
                )

        return binary_classifier_holder

    def get_classifier_name_to_key(self, classifier_name):
        """Return the standard key (str) of a classifier.

        Parameters
        ----------
        classifier_name : str
            Name of classification algorithm to use.

        Returns
        -------
        key : str
            Key to access trained classifier objects.
        """
        if classifier_name.lower() in LinearNDInterpolator_names:
            key = "LinearNDInterpolator"
        elif classifier_name.lower() in RBF_names:
            key = "RBF"
        elif classifier_name.lower() in GaussianProcessClassifier_names:
            key = "GaussianProcessClassifier"
        else:
            print("No classifiers with name '{0}'.".format(classifier_name))
            return None
        return key

    def __remove_nans(self, trans_probs, verbose=False):
        """Given an array of probabilities, remove the nans.

        Return where not nan
        which is useful for finding good or bad values.

        """
        bool_row_index_where_nan = np.isnan(trans_probs.T[0])
        # where_nan = np.where(bool_row_index_where_nan)[0]
        where_not_nan = np.where(~bool_row_index_where_nan)[0]
        # how_many_rows = len(trans_probs.T[0])
        how_many_nans = np.sum(bool_row_index_where_nan)

        clean_trans_probs = np.array([trans_probs[i] for i in where_not_nan])
        if verbose:
            print("Nans omitted: {0}".format(how_many_nans))
        return clean_trans_probs, where_not_nan

    def return_probs(self, classifier_name, test_input, verbose=False):
        """Return probability that a given input corresponds to a class.

        The probability is calculated using trained classifiers.

        Parameters
        ----------
        classifier_name : str
            Name of classifier to train.
        test_input : ndarray
            N dimensional inputs to be classified.
            The shape should be N_points x N_dimensions.
        verbose : optional, bool
            Print useful information.

        Returns
        -------
        normalized_probs : ndarray
            Array holding the normalized probability for a point to be in any
            of the possible classes. Shape is N_points x N_classes.
        where_not_nan: ndarray
            Indicies of the test inputs that did not result in nans.

        """
        if isinstance(test_input, list):
            test_input = np.array(test_input)
        if test_input.ndim == 1:
            test_input = np.array([test_input])

        probs = []
        if not bool(self._interpolators_) and not bool(
            self._cv_interpolators_
        ):  # if empty
            raise Exception("\n\nNo trained interpolators exist.")

        # convert the user input shorthand into a valid key for dict
        classifier_name = self.get_classifier_name_to_key(classifier_name)

        if self.__train_cross_val:
            interpolators = self._cv_interpolators_[classifier_name].items()
        else:
            interpolators = self._interpolators_[classifier_name].items()

        for key, interp in interpolators:
            # - each interpolator is on a different class
            if classifier_name == "RBF":
                argList = []
                for col in range(len(test_input[0])):
                    argList.append(test_input.T[col])
                uncleaned_data = interp(
                    *argList
                )  # inherent overshoot occurs in rbf interpolation
                step_1 = np.where(
                    uncleaned_data > 1, 1, uncleaned_data
                )  # values above 1 replace with 1
                cleaned_data = np.where(
                    step_1 < 0, 0, step_1
                )  # negative values set to 0
                probs.append(cleaned_data)
            elif classifier_name == "GaussianProcessClassifier":
                # print(key, interp.predict(test_input),
                #       interp.predict_proba(test_input ).T[1])
                probs.append(interp.predict_proba(test_input).T[1])
                # The [1] is selecting the second output from predict_proba for
                # all points. The second output being the probability that it
                # is the current class. This is a similar form to all the other
                # classifier output.
            elif classifier_name == "LinearNDInterpolator":
                unfiltered_output = interp(
                    test_input
                )  # can return nan if out of bounds
                probs.append(unfiltered_output)
            else:
                print("Name not recognized: '{0}'".format(classifier_name))
        # for loop - all test points do one interpolator
        # 1st array in probs = [ <class 0 prob on 1st test point>,
        #                        <class 0 prob on 2nd test point>,... ]
        # to get a row where each element is P for a diff class -> transpose
        # probs

        trans_probs = np.array(probs).T
        clean_trans_probs, where_not_nan = self.__remove_nans(
            trans_probs, verbose=verbose
        )  # remove nans

        try:
            totals_per_input = np.sum(clean_trans_probs, axis=1)
            normalized_probs = clean_trans_probs / totals_per_input[:,
                                                                    np.newaxis]
        except Exception:
            normalized_probs = [[0]]

        return normalized_probs, where_not_nan

    def get_class_predictions(self, classifier_name, test_input,
                              return_ids=True):
        """Return the class predictions.

        The predictions are in the form of class IDs or the original
        classification key. This method also returns the probability
        of the class that was predicted.

        Parameters
        ----------
        classifier_name : str
            Name of classification algorithm to use.
        test_input : ndarray
            Input values to predict. Same shape as input data.
        return_ids : bool, optional
            If True (default), return class IDs.
            Else, return the original classification keys.

        Returns
        -------
        pred_class_ids : array
            Predicted class IDs given test input.
        max_probs : array
            Probability the classifier gives for the chosen class.
        where_not_nan : array
            Inidices where there are no nans (from LinearNDInterpolator).
            You may use this to pick out which input data gives a valid
            classification.

        """
        all_probs, where_not_nan = self.return_probs(classifier_name,
                                                     test_input)
        max_probs = np.max(all_probs, axis=1)
        pred_class_ids = np.argmax(all_probs, axis=1)

        if return_ids:
            return pred_class_ids, max_probs, where_not_nan
        else:
            pred_classes = [self.class_id_mapping[i] for i in pred_class_ids]
            return pred_classes, max_probs, where_not_nan

    def get_cross_val_data(self, alpha):
        """Randomly sample the data set and seperate training and test data.

        Parameters
        ----------
        alpha : float
            Fraction of data set to use for training. (0.05 = 5% of data set)

        Returns
        -------
        sorted_rnd_int_vals : array
            Array indicies for data used to train interpolators.
        cv_test_input_data : array
            Input test data to perform cross validation.
        cv_test_output_data : array
            Output test data to perform cross validation.
        """
        num_points = int(len(self.input_data) * alpha)
        # rnd_input_train = []
        # rnd_outout_train = []
        rnd_int_vals = []
        rnd_int_set = set()

        ct = 0
        while len(rnd_int_vals) < num_points and ct < 1e7:
            rnd_int = int(np.random.random() * len(self.input_data))

            if rnd_int not in rnd_int_set:
                rnd_int_vals.append(rnd_int)
                rnd_int_set.add(rnd_int)
            ct += 1

        sorted_rnd_int_vals = sorted(rnd_int_vals)

        # Random training data
        # self.cross_val_train_input_data = \
        #    self.input_data[sorted_rnd_int_vals,:]
        # self.cross_val_train_class_data = \
        #    np.argmax(self.binary_class_data.T[sorted_rnd_int_vals,:], axis=1)

        test_int_vals = []
        for i in range(len(self.input_data)):
            if i in sorted_rnd_int_vals:
                pass
            else:
                test_int_vals.append(i)

        # The remainder which will be used to test fits
        cv_test_input_data = self.input_data[test_int_vals, :]
        cv_test_output_data = np.argmax(
            self.binary_class_data.T[test_int_vals, :], axis=1
        )

        return sorted_rnd_int_vals, cv_test_input_data, cv_test_output_data

    def cross_validate(self, classifier_names, alpha, verbose=False):
        """Cross validate classifiers on data from TableData object.

        For each iteration, the classifiers specified are all trained
        and tested on the same random subset of data.

        Parameters
        ----------
        classifier_names : array
            Names of classifiers to train.
        alpha : float
            Fraction of data set to use for training. (0.05 = 5% of data set)
        verbose : bool, optional
            Print statements with more information while training.

        Returns
        -------
        percent_correct: ndarray
            Percent correct classification on (1-alpha)% of the data set.
            Element order matches the order of classifier_names.
        time_to_train : ndarray
            Time to train classifiers on a data set. Element order matches
            the order of classifier_names.

        """
        self.__train_cross_val = True

        (
            train_data_indicies,
            cv_test_input_data,
            cv_test_output_data,
        ) = self.get_cross_val_data(alpha)

        classifier_names = [
            self.get_classifier_name_to_key(x) for x in classifier_names
        ]

        if verbose:
            print(
                "alpha: {0}, num_training_points: {1}".format(
                    alpha, len(train_data_indicies)
                )
            )

        time_to_train = []
        # Train classifiers
        start_time = time.time()
        for name in classifier_names:
            self.train(name, di=train_data_indicies)
            time_to_train.append(time.time() - start_time)

        predicted_class_ids = []
        where_not_nans_holder = []
        num_corr_holder = []
        # Test classification
        for name in classifier_names:
            pred_ids, probs, where_not_nan = self.get_class_predictions(
                name, cv_test_input_data
            )

            predicted_class_ids.append(pred_ids)
            where_not_nans_holder.append(where_not_nan)
            # pred_ids already has nans omitted...
            # the output data needs to have those removed...
            num_correct = np.sum(
                pred_ids == cv_test_output_data[where_not_nan])
            num_corr_holder.append(num_correct)

        percent_correct = []
        for j in range(len(num_corr_holder)):
            top = num_corr_holder[j]
            bot = len(cv_test_input_data[where_not_nans_holder[j]])
            percent_correct.append(top / bot * 100.0)

        if verbose:
            print("% Correct      Interpolator")
            print("---------    ----------------")
            for i in range(len(classifier_names)):
                offset = 5 + len(classifier_names[i])
                print(
                    " {1:.3f} {0}".format(
                        classifier_names[i].rjust(offset), percent_correct[i]
                    )
                )
            print("")

        self.__train_cross_val = False
        return np.array(percent_correct), np.array(time_to_train)

    def make_cv_plot_data(self, interp_type, alphas, N_iterations,
                          folder_path="cv_data/"):
        """Script for running many instances of the method `cross_validate()`.

        Cross validation score and timing data produced are saved locally.

        ! Time to train GaussianProcessClassifier becomes large for num
        training points > 1000. !

        Files saved every 5 iterations to prevent loss of data for large
        N_iterations.
        Known expection occurs in GP classifier for low alpha due to data set
        with only one class.

        Parameters
        ----------
        interp_type : array_str
            Names of classifiers to train.
        alphas : array_floats
            Fractions of data set to use for training. (0.05 = 5% of data set)
            (ex. [0.01, 0.02, ...])
        N_iterations : int
            Number of iterations to run cross validation at a given alpha.
        folder_path : str
            Folder path where to save cross validation and timing data
            ("your_folder_path/").

        Returns
        -------
        None

        """
        N_total_points = len(self.input_data)

        for frac in alphas:
            print("alpha: {0}".format(frac))
            print("N training points: {0}".format(int(N_total_points * frac)))

            cross_df = pd.DataFrame(columns=interp_type)
            time_df = pd.DataFrame(columns=interp_type)

            good_cv_count = 0
            # iteration_counter = 0
            while good_cv_count < (N_iterations):
                try:
                    percent_correct, time_to_train = self.cross_validate(
                        interp_type, frac
                    )

                    # .loc[i] is setting a row
                    cross_df.loc[good_cv_count] = percent_correct
                    time_df.loc[good_cv_count] = time_to_train

                    if good_cv_count % 5 == 0:
                        print("\t iter", good_cv_count)
                        cross_df.to_csv(folder_path
                                        + "cross_val_data_f%s" % (str(frac)))
                        time_df.to_csv(folder_path
                                       + "timing_data_f%s" % (str(frac)))

                    good_cv_count += 1

                except Exception as ex:
                    print("!expection!", ex)

                if good_cv_count == 1:
                    print("Time to run %i at alpha %.3f" % (N_iterations,
                                                            frac))
                    print("~%.2f seconds" % (np.sum(time_to_train)
                                             * N_iterations))

            cross_df.to_csv(folder_path + "cross_val_data_f%s" % (str(frac)))
            time_df.to_csv(folder_path + "timing_data_f%s" % (str(frac)))
        print(" ------------ DONE ------------ ")
        return

    def get_rnd_test_inputs(self, N, other_rng=dict(), verbose=False):
        """Produce randomly sampled 'test' inputs inside domain of input_data.

        Parameters
        ----------
        N : int
            Number of test inputs to return.
        other_rng : dict, optional
            Change the range of random sampling in desired axis. By default,
            the sampling is done in the range of the training data. The axis is
            specified with an integer key in [0,N-1] mapping to a list
            specifying the range. (e.g. {1:[my_min, my_max]})
        verbose : bool, optional
            Print diagnostic information.

        Returns
        -------
        rnd_test_points : ndarray
            Test points randomly sampled in the range of the training data
            in each axis unless otherwise specified in 'other_rng'.
            Has the same shape as input data from TableData.

        """
        num_axis = len(self.input_data[0])

        # find max and min in each axis
        a_max = []
        a_min = []
        for i in range(num_axis):
            a_max.append(max(self.input_data.T[i]))
            a_min.append(min(self.input_data.T[i]))
        # sample N points between max & min in each axis
        axis_rnd_points = []
        for i in range(num_axis):
            if i in other_rng:
                b_min, b_max = other_rng[i]
            else:
                b_min, b_max = a_min[i], a_max[i]
            if verbose:
                print("{0} - min: {1}, max: {2}".format(i, b_min, b_max))

            r = np.random.uniform(low=b_min, high=b_max, size=int(N))

            # this reshape is necessary to concatenate
            axis_rnd_points.append(r[:, np.newaxis])

        # now put the random points back together with same shape as input_data
        rnd_test_points = np.concatenate(axis_rnd_points, axis=1)
        return rnd_test_points

    def make_max_cls_plot(self, classifier_name, axes_keys, other_rng=dict(),
                          N=4000, **kwargs):
        """Make the maximum classification probablity plot.

        Not generalized yet to slice along redundant axes.

        """
        vals = self.get_rnd_test_inputs(N, other_rng=other_rng)
        class_vals, probs, where_not_nan = self.get_class_predictions(
            classifier_name, vals, return_ids=False
        )

        fig, (cls_data, prob) = plt.subplots(
            1, 2, figsize=(12, 4), gridspec_kw={"width_ratios": [1.0, 1.2]}
        )

        input_data = self._TableData_.get_data(
            what_data="input", return_df=True)
        output_cls_data = self._TableData_.get_class_data(
            what_data="class_col")
        x_data = input_data[axes_keys[0]].values
        y_data = input_data[axes_keys[1]].values
        cls_data.set_title("Input data from TableData")
        cls_data.set_xlabel(axes_keys[0])
        cls_data.set_ylabel(axes_keys[1])
        color_dict = OrderedDict()
        for j, color_str in enumerate(self._TableData_._class_colors_):
            color_dict[j] = color_str

        cls_data.scatter(x_data, y_data, c=[color_dict[i]
                                            for i in output_cls_data])
        # You need to use a mapping to get classes to colors
        # I should add this to data.py
        prob.set_xlabel(axes_keys[0])
        prob_plt = prob.scatter(vals.T[0], vals.T[1], c=probs, **kwargs)
        cb = fig.colorbar(prob_plt)
        cb.ax.set_ylabel("Maximum classification probability")
        prob.set_title("{}".format(
            self.get_classifier_name_to_key(classifier_name)))
        return fig, (cls_data, prob)
