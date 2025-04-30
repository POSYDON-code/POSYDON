"""Module for handling data for PSY-CRIS."""


__authors__ = [
    "Kyle Akira Rocha <kylerocha2024@u.northwestern.edu>",
    "Scott Coughlin <scottcoughlin2014@u.northwestern.edu>",
]


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd
import math
import time

from collections import OrderedDict
from sklearn.neighbors import NearestNeighbors


def calc_avg_dist(data, n_neighbors, neighbor=None):
    """Get the average distance to the nearest neighbors in the data set.

    (NearestNeighbors from sklearn.neighbors)

    Parameters
    ----------
    data : ndarray
        Data to train the NearestNeighbors class on.
    n_neighbors : int
        Number of neighbors to use when finding average distance.
    neighbor : instance of NearestNeightbors class
        For passing your own object.

    Returns
    -------
    avg_dist : array
        The average distance between nearest neighbors.
    g_indi : array
        Indicies that correspond to the nearest neighbors.
    """
    if neighbor:
        neigh = neighbor
    else:
        neigh = NearestNeighbors()
    neigh.fit(data)

    # Since we pass the original data, the first nearest point is itself. If we
    # want 2 nearest neighbors, we ask for 3 because the first will always be 0
    dist, indicies = neigh.kneighbors(data, n_neighbors=(n_neighbors + 1))

    g_dist = (dist.T[1:len(dist)]).T
    g_indi = (indicies.T[1:len(indicies)]).T
    avg_dist = np.mean(g_dist, axis=1)
    return avg_dist, g_indi


def calc_avg_p_change(data, where_nearest_neighbors,
                      undefined_p_change_val=None):
    """Calculate the average fractional change in a given data set.

    The method uses the N nearest neighbors (calculated beforehand).

    Parameters
    ----------
    data : ndarray
        Data set to calculate percent change.
    where_nearest_neighbors : dict
        Indicies in data for the n nearest neighbors in input space.
    undefined_p_change_val : optional
        For output with an undefined percent change (zero value), this kwarg
        defines what is put in its place. Defaults to nan.

    Returns
    -------
    avg_p_change_holder : ndarray
        Each element conatins the average percent change for a given number of
        neighbors. If any output values are 0, then a nan is put in place of a
        percent change.

    """
    where_zero = np.where(data == 0)[0]
    where_not_zero = np.where(data != 0)[0]
    if len(where_zero) == len(data):
        return None

    avg_p_change_holder = []
    for n_, indicies in where_nearest_neighbors.items():
        diff_holder = []
        for i in range(n_):
            nearest_data = data[indicies.T[i]]
            diff_holder.append(data - nearest_data)

        diffs = np.array(diff_holder).T
        avg_diffs = np.mean(diffs, axis=1)
        avg_p_change = np.zeros(data.shape)
        avg_p_change[where_not_zero] = (
            abs(avg_diffs[where_not_zero]) / data[where_not_zero]
        )
        if undefined_p_change_val is None:
            avg_p_change[where_zero] = np.nan
        else:
            avg_p_change[where_zero] = undefined_p_change_val

        avg_p_change_holder.append(avg_p_change)

    return np.array(avg_p_change_holder)


class TableData:
    """
    For managing data sets used for classification and regression.

        Reads tables of simulation data where a single row represents one
        simulation. Each column in a row represents different inputs
        (initial conditions) and outputs (result, continuous variables).
        If using multiple files, each file is assumed to have the same columns.
        You may also directly load a pandas DataFrame instead of reading in
        files.

        Example data structure expected in files or pandas DataFrame:

        0  input_1  input_2  outcome  output_1   output_2  output_3 ...
        1    1.5      2.6      "A"       100       0.19       -     ...
        2    1.5      3.0      "B"        -          -        -     ...
        3    2.0      2.6      "C"        -          -        6     ...
        ...

        The above table has dashes '-' in output columns to indicate NaN
        values. You may have a similar structure if different classes have
        fundamentally different outputs.
    """

    def __vb_helper(self, verbose_bool, info_string):
        """Help clean up verbose print statements and storing info.

        By default, store the passed info_string regardless of verbose_bool.

        """
        if verbose_bool:
            print(info_string)
        self._for_info_.append(info_string)

    def __init__(
        self,
        table_paths,
        input_cols,
        output_cols,
        class_col_name,
        my_DataFrame=None,
        omit_vals=None,
        omit_cols=None,
        subset_interval=None,
        verbose=False,
        my_colors=None,
        neighbor=None,
        n_neighbors=None,
        undefined_p_change_val=None,
        read_csv_kwargs={},
        **kwargs
    ):
        """Initialize the TableData instance.

        Parameters
        ----------
        table_paths : list
            List of file paths to read in as data.
            if None, a pandas DataFrame is used instead
        input_cols : list
            List of names of the columns which will be considered 'input'.
        output_cols : list
            List of names of the columns which will be considered 'output'.
            This should include the class column name.
        class_col_name : str
            Name of column which contains classification data.
        my_DataFrame : pandas DataFrame, optional
            If given, use this instead of reading files.
        omit_vals : list, optional
            Numerical values that you wish to omit from the entire data set.
            If a row contains the value, the entire row is removed.
            (For example you may want to omit all rows if they contain "-1" or
            "failed".)
        omit_cols : list, optional
            Column names that you wish to omit from the data set.
        subset_interval : array, optional
            Use some subset of the data files being loaded in.
            An array with integers indicating the rows that will be kept.
        my_colors : list, optional
            Colors to use for classification plots.
        n_neighbors : list, optional
            List of integers that set the number of neighbors to use to
            calculate average distances. (default None)
        neighbor : instance of sklearn.neighbors.NearestNeighbors, optional
            To use for average distances. See function 'calc_avg_dist()'.
        undefined_p_change_val : optional, float
            Sets the undefined value used when calculating percent change fails
            due to zero values in the output data. Default uses nan.
        verbose : bool, optional
            Print statements with extra info.
        read_csv_kwargs : dict, optional
            Kwargs passed to the pandas function 'read_csv()'.
        **kwargs
            Extra kwargs

        """
        start_time = time.time()
        self._for_info_ = []  # storing info strings

        # --------- data pre-processing ---------
        self._df_list_ = []  # data frame list
        self._df_index_keys_ = []  # index keys (for rows in _full_data_)
        self._files_ = table_paths
        # assumed to be one column name + string
        self.class_col_name = class_col_name

        if isinstance(my_DataFrame, pd.DataFrame):
            self.__vb_helper(verbose, "Using loaded Pandas DataFrame")
            self._full_data_ = my_DataFrame
        else:
            # Read in all data files and add to _df_list_
            info_str_01 = "Reading in data from {0} file(s).".format(
                len(table_paths))
            self.__vb_helper(verbose, info_str_01)

            for num, path in enumerate(table_paths):
                info_str_02 = "\t'{0}'".format(path)
                self.__vb_helper(verbose, info_str_02)

                df = pd.read_csv(path, **read_csv_kwargs)
                self._df_list_.append(df)
                self._df_index_keys_.append("df" + str(num))

            info_str_03 = "Finished reading data.\n"
            self.__vb_helper(verbose, info_str_03)

            # _df_index_keys_ setting index of Nth file
            # the data is from with 'dfN'.
            self._full_data_ = pd.concat(self._df_list_, join="outer",
                                         ignore_index=False,
                                         keys=self._df_index_keys_, sort=True)
        # remove rows and columns with unwanted data
        if omit_vals is not None:
            ct = 0  # counter
            for j, val in enumerate(omit_vals):
                row_where_omit_val = np.unique(np.where(
                    self._full_data_.values == val)[0])
                df_keys_for_omit = self._full_data_.index[row_where_omit_val]
                self._full_data_ = self._full_data_.drop(df_keys_for_omit,
                                                         axis=0)
                ct += len(row_where_omit_val)

                info_str_04 = " - Removed {0} rows containing: {1}".format(
                    len(row_where_omit_val), omit_vals[j]
                )
                self.__vb_helper(verbose, info_str_04)
            info_str_05 = "Removed a total of {0} rows.".format(ct)
            self.__vb_helper(verbose, info_str_05)

        # remove entire columns
        if omit_cols is not None:
            self._full_data_ = self._full_data_.drop(columns=omit_cols)
            info_str_06 = "Removed columns: {0}".format(omit_cols)
            self.__vb_helper(verbose, info_str_06)

        # use a subset of the original data table
        if subset_interval is not None:
            len_original_data = len(self._full_data_)
            self._full_data_ = self._full_data_.iloc[subset_interval, :]
            info_str_07 = ("--Using Subset--\n{0} percent of total data set.".
                           format(len(self._full_data_) / len_original_data
                                  * 100))
            self.__vb_helper(verbose, info_str_07)

        info_str_08 = "Total number of data points: {0}\n".format(
            len(self._full_data_))
        self.__vb_helper(verbose, info_str_08)

        # column names in fully concatenated data
        self.col_names = np.array(self._full_data_.columns)

        # input and output data
        for usr_input in [input_cols, output_cols]:
            for a_name in usr_input:
                if a_name not in self.col_names:
                    info_str_09 = ("\tWarning: No columns found with name "
                                   "'{0}'. Check your data.".format(a_name))
                    self.__vb_helper(verbose, info_str_09)
                    # raise Warning( info_str_09 ) TODO
        input_cols = [i for i in input_cols if i in self.col_names]
        output_cols = [i for i in output_cols if i in self.col_names]
        self._input_ = self._full_data_[input_cols]
        self._output_ = self._full_data_[output_cols]

        # max and min for input values
        self._max_input_vals = []
        self._min_input_vals = []
        for data_column in self._input_.values.T:
            self._max_input_vals.append(max(data_column))
            self._min_input_vals.append(min(data_column))

        # --------- classification variables ---------
        try:
            self._unique_class_keys_ = np.unique(
                self._full_data_[class_col_name])
            self._class_dtype_ = type(self._unique_class_keys_[0])
        except KeyError as class_col_name:
            info_str_10 = "Class column {0} not in {1}".format(
                class_col_name, np.array(self._full_data_.keys()).astype(str)
            )
            self.__vb_helper(False, info_str_10)
            raise KeyError(info_str_10)

        self.num_classes = len(self._unique_class_keys_)
        self.class_ids = np.arange(0, self.num_classes, 1, dtype=int)

        info_str_11 = (
            "Input columns: {0}".format(len(input_cols))
            + "\n"
            + "Output columns: {0}".format(len(output_cols))
            + "\n"
            + "Unique classes found in '{0}': {1}".format(
                class_col_name, self.num_classes
            )
        )
        self.__vb_helper(verbose, info_str_11)

        # mapping dict - forward & backward
        self._class_id_mapping_ = OrderedDict()
        for i in self.class_ids:
            self._class_id_mapping_[i] = self._unique_class_keys_[i]
            self._class_id_mapping_[self._unique_class_keys_[i]] = i

        # classification column replaced with class_id
        self._class_col_ = self._full_data_[class_col_name].values.astype(
            self._class_dtype_
        )
        self._class_col_to_ids_ = []
        for cl in self._class_col_:
            self._class_col_to_ids_.append(self._class_id_mapping_[cl])

        if my_colors is not None:
            self._class_colors_ = my_colors
            self.__vb_helper(verbose, "Using custom class colors.")
        else:
            self._class_colors_ = ["#EC6666", "#90A245", "#F5C258", "#1668E8",
                                   "#473335", "#98C0CB", "C0", "C1", "C2",
                                   "C3", "C4", "C5", "C6", "C7", "C8"]
            self.__vb_helper(verbose, "Using default class colors.")

        # --------- regression variables ---------

        # make input and output dicts for each class
        self._regr_inputs_ = OrderedDict()
        self._regr_outputs_ = OrderedDict()
        for cls_name in self._unique_class_keys_:
            rows_where_class = np.where(self._output_[class_col_name]
                                        == cls_name)
            if len(rows_where_class[0]) == 0:
                info_str_extra_1 = ("\tWarning: no output with class name {0}."
                                    .format(cls_name))
                self.__vb_helper(verbose, info_str_extra_1)
            self._regr_inputs_[cls_name] = self._input_.iloc[rows_where_class]
            self._regr_outputs_[cls_name] = self._output_.iloc[
                rows_where_class]

        info_str_12 = ("\nFinding values to regress:\n"
                       + "Num output(s) \t Class Name")
        self.__vb_helper(verbose, info_str_12)

        # find valid regression data and link it to a class
        self.regr_names = self._output_.columns
        self._num_outputs_to_regr_ = []
        self._regr_dfs_per_class_ = OrderedDict()
        # for each class - find columns which can be converted to floats
        for i, tuple_output in enumerate(self._regr_outputs_.items()):
            output_per_class = tuple_output[1]  # 0 returns key; 1 returns data

            # this makes sure the classification column is not part of
            # regression data in the case classes are integers or somthing
            where_no_classification = np.where(
                output_per_class.columns != class_col_name
            )[0]

            cols_with_float = []
            bad_cols = []
            # go through first row and try to convert each element to float
            # - also check if nan
            for col_num, val in zip(
                where_no_classification,
                output_per_class.iloc[0, where_no_classification],
            ):
                try:
                    converted_val = float(
                        val
                    )  # if fails -> straight to except (skips next line)
                    if pd.isna(converted_val):
                        bad_cols.append(col_num)
                    else:
                        cols_with_float.append(col_num)
                except Exception:
                    bad_cols.append(col_num)

            info_str_13 = "%7i \t '%s'" % (
                len(cols_with_float),
                self._unique_class_keys_[i],
            )
            self.__vb_helper(verbose, info_str_13)

            self._num_outputs_to_regr_.append(len(cols_with_float))
            regression_vals_df = output_per_class.iloc[
                :, cols_with_float
            ]  # has the regression DataFrame for a given class
            # if regressable elements - link class with the df of valid floats,
            # else - None
            if len(cols_with_float) > 0:
                self._regr_dfs_per_class_[
                    self._unique_class_keys_[i]
                ] = regression_vals_df
            else:
                self._regr_dfs_per_class_[self._unique_class_keys_[i]] = np.nan

        # take Nearest Neighbors differences
        self._undefined_p_change_val_ = undefined_p_change_val
        if n_neighbors is not None:
            info_str_14 = ("\nCalculate Average Distances & "
                           "Average Percent Change")
            self.__vb_helper(verbose, info_str_14)

            self._n_neighbors_ = [int(i) for i in n_neighbors]
            self._avg_dist_dfs_per_class_ = (
                OrderedDict()
            )  # stores the average distances

            # We need distances in input space, then compare %diff in output
            # space
            for key, input_df_per_class in self._regr_inputs_.items():
                self._avg_dist_dfs_per_class_[key] = OrderedDict()
                regr_data = self._regr_dfs_per_class_[
                    key
                ]  # either a DataFrame or np.nan

                if isinstance(regr_data, pd.DataFrame):
                    self.__vb_helper(verbose, "class: '{0}'".format(key))

                    # find nearest neighbors in input space
                    where_nearest_neighbors = OrderedDict()
                    for n_ in self._n_neighbors_:
                        try:
                            avg_dist, indicies = calc_avg_dist(
                                input_df_per_class.values, n_,
                                neighbor=neighbor)
                        except ValueError as err:
                            self.__vb_helper(verbose,
                                             "Skipping n = {0}\nError: {1}".
                                             format(n_, err))
                            avg_dist = None
                            continue

                        where_nearest_neighbors[n_] = indicies

                    for _key, _val in regr_data.items():  # regression data
                        data = np.array(_val.values, dtype=float)
                        where_zero = np.where(data == 0)[0]
                        # where_not_zero = np.where(data != 0)[0]
                        if len(where_zero) > 0:
                            info_str_15 = "\t -- {0} zeros in '{1}'".format(
                                len(where_zero), _key
                            )
                            self.__vb_helper(verbose, info_str_15)
                            if len(where_zero) == len(data):
                                self.__vb_helper(
                                    verbose, "Skipping percent change for "
                                    "{0}... All data is zero.".format(_key))
                                continue
                        # Take percent difference
                        avg_p_change = calc_avg_p_change(
                            data, where_nearest_neighbors,
                            undefined_p_change_val=self.
                            _undefined_p_change_val_)
                        if avg_p_change is None:
                            self.__vb_helper(verbose, "None in avg_p_change!? "
                                             "This shouldn't happen...")
                        else:
                            # update into regr_dfs_per_class
                            # - all data available for regression
                            my_kwargs = OrderedDict()
                            for i in range(np.shape(avg_p_change)[0]):
                                new_col_str = "APC{0}_{1}".format(
                                    self._n_neighbors_[i], _key
                                )
                                self.__vb_helper(verbose, "\t" + new_col_str)
                                my_kwargs[new_col_str] = avg_p_change[i]
                            self._regr_dfs_per_class_[key] = (
                                self._regr_dfs_per_class_[key].assign(
                                    **my_kwargs))

                        self._avg_dist_dfs_per_class_[key][_key] = avg_dist

                else:
                    self.__vb_helper(
                        verbose, "No regression data in '{0}'.".format(key)
                    )

        info_str_17 = "::: TableData created in {0:.2f} seconds :::".format(
            time.time() - start_time
        )
        self.__vb_helper(verbose, info_str_17)

    def find_n_neighbors(self, input_data, n_neighbors, neighbor=None,
                         return_avg_dists=False, **kwargs):
        """Find the N nearest neighbors of a given set of arbitrary points.

        Given a set of arbitrary input points, find the N nearest neighbors
        to each point, not including themselves. Can also return average
        distance from a point to its N nearest neighbors.

        Parameters
        ----------
        input_data : ndarray, pandas DataFrame
            Data points where their nearest neighbors will be found
        n_neighbors : list
            List of integers with number of neighbors to calculate
        neighbor : instance of NearestNeightbors class
            For passing your own object.
        return_avg_dists : bool, optional
            If True, return the a dictionary with average distances between
            the N nearest neighbors listed in 'n_neighbors'.
        **kwargs
            class_key : str
                If a DataFrame is given, specifies class column name.

        Returns
        -------
        where_nearest_neighbors : dict
            Dictionary containing the n nearest neighbors for every point
            in the input data.
        avg_distances : dict
            Returned if return_avg_dists is True.
            The average distances between the nearest neighbors.

        """
        class_key = kwargs.pop("class_key", None)
        if class_key is not None:
            input_data = self._regr_dfs_per_class_[class_key]
        if isinstance(input_data, pd.DataFrame):
            input_data = input_data.values
        elif isinstance(input_data, list):
            input_data = np.array(input_data)
        elif isinstance(input_data, np.ndarray):
            pass
        else:
            raise ValueError(
                "input_data must be pandas DataFrame or numpy array")

        where_nearest_neighbors = OrderedDict()
        avg_distances = OrderedDict()
        for n_ in n_neighbors:
            avg_dist, indicies = calc_avg_dist(input_data, n_,
                                               neighbor=neighbor)
            where_nearest_neighbors[n_] = indicies
            avg_distances[n_] = avg_dist
        if return_avg_dists:
            return where_nearest_neighbors, avg_distances
        else:
            return where_nearest_neighbors

    def get_data(self, what_data="full", return_df=False):
        """Get all data contained in TableData object.

        The data is returned after omission of columns and rows containing
        specified values (if given) and taking a subset (if given) of the
        original data set read in as a csv or given directly as a pandas
        DataFrame. (Before data processing for classification and regression.)

        Parameters
        ----------
        what_data: str, list, optional
            Default is 'full' with other options 'input', or 'output'.
            'full' - original data table (after omission and subsets)
            'input' - only data identified as inputs from 'full' data set
            'output' - only data identified as outputs from 'full' data set
        return_df: bool, optional
            If True, return a pandas DataFrame object.
            If False (default), return a numpy array.

        Returns
        -------
        data: tuple, ndarray or DataFrame
            Data before classification and regression data sorting is done.
            Will return tuple of len(what_data) if a list is passed.

        """
        data_dict = {
            "input": self._input_,
            "output": self._output_,
            "full": self._full_data_,
        }
        if isinstance(what_data, str):
            my_data = data_dict[what_data.lower()]
        elif isinstance(what_data, list):
            my_data = tuple([data_dict[i.lower()] for i in what_data])
        else:
            raise ValueError("'{0}' not supported. Try {1}.".
                             format(what_data, data_dict.keys()))
        if return_df:
            return my_data
        else:
            return my_data.values

    def get_binary_mapping_per_class(self,):
        """Get binary mapping (0 or 1) of the class data.

        Get binary mapping (0 or 1) of the class data for each unique
        classification. For each classification, a value of 1 is given
        if the class data matches that classification. If they do not
        match, then a value of 0 is given.

        Example:
        classifications -> A, B, C
        class data      -> [ A, B, B, A, B, C ]
        binary mapping  -> [[1, 0, 0, 1, 0, 0]   (for class A)
                            [0, 1, 1, 0, 1, 0]   (for class B)
                            [0, 0, 0, 0, 0, 1]]  (for class C)

        Returns
        -------
        binary_class_data: ndarray
            N by M array where N is the number of classes and M is the
            number of classifications in the data set.
            Order is determined by '_unique_class_keys_'.

        """
        cls = self._class_col_.copy()
        # create a different array for each class - one against all
        binary_class_data = []
        for i in self.class_ids:
            where_class_is = np.where(
                cls == self._unique_class_keys_[i], 1, 0
            )  # 1 for True, 0 for False
            binary_class_data.append(np.concatenate(where_class_is, axis=None))
        return np.array(binary_class_data)

    def _return_data_(self, name, verbose=False):
        """Return a regular or hidden class variable."""
        hidden_name = "_" + name + "_"
        if name in self.__dict__:
            return self.__dict__[name]
        elif hidden_name in self.__dict__:
            return self.__dict__[hidden_name]
        else:
            if verbose:
                print(
                    "{0}, {1} not in class variables: \n{2}".format(
                        name, hidden_name, self.__dict__
                    )
                )
            return None

    def get_class_data(self, what_data="full"):
        """Get data related to classification.

        Parameters
        ----------
        what_data, str, list, optional
            'class_col' (array)
                - Original classification data.
            'unique_class_keys' (array)
                - Unique classes found in the classification data.
            'class_col_to_ids' (array)
                - Original classification data replaced with their
                  respective class IDs (integers).
            'class_id_mapping' (dict)
                - Mapping between a classification from the original
                  data set and its class ID.
            'binary_data' (ndarray)
                - Iterating over the unique classes, classification data
                  is turned into 1 or 0 if it matches the given class.
                  See method 'get_binary_mapping_per_class()'.
            'full' (tuple)
                - All options listed above. (Default)

        Returns
        -------
        class_data: tuple, ndarray, dict
            An object containing the specified classification data.
            Will return tuple of len(what_data) if a list is passed. Default: 5
        """
        binary_class_data = self.get_binary_mapping_per_class()
        data_dict = {
            # simply the class column
            "class_col": self._class_col_,
            # unique values in class column
            "unique_class_keys": self._unique_class_keys_,
            # class column but class strings turned into integers (class IDs)
            "class_col_to_ids": self._class_col_to_ids_,
            # maps between a class ID and the original class string
            "class_id_mapping": self._class_id_mapping_,
            # for all classes, 1 where that class, 0 else
            "binary_data": binary_class_data,
            "full": tuple(
                [
                    self._class_col_,
                    self._unique_class_keys_,
                    self._class_col_to_ids_,
                    self._class_id_mapping_,
                    binary_class_data,
                ]
            ),
        }
        if isinstance(what_data, str):
            class_data = data_dict[what_data.lower()]
        elif isinstance(what_data, list):
            class_data = tuple([data_dict[i.lower()] for i in what_data])
        else:
            raise ValueError("'{0}' not supported. Try {1}.".
                             format(what_data, data_dict.keys()))
        return class_data

    def get_regr_data(self, what_data="full"):
        """Get data related to regression all sorted by class in dictionaries.

        Parameters
        ----------
        what_data: str, list, optional
            'input' - For each class, the input data with no cleaning.
            'raw_output' - For each class, the output data with no cleaning.
            'output' - For each class, the cleaned output data.
            'full' - All options listed above in that respective order.

        Returns
        -------
        data: tuple, ndarray or DataFrame
            An object containing the specified regression data.
            Will return tuple of len(what_data) if a list is passed. Default: 3
        """
        data_dict = {
            "input": self._regr_inputs_,
            "raw_output": self._regr_outputs_,
            "output": self._regr_dfs_per_class_,
            "full": tuple([self._regr_inputs_,
                           self._regr_outputs_,
                           self._regr_dfs_per_class_]),
        }
        if isinstance(what_data, str):
            regr_data = data_dict[what_data.lower()]
        elif isinstance(what_data, list):
            regr_data = tuple([data_dict[i.lower()] for i in what_data])
        else:
            raise ValueError("'{0}' not supported. Try {1}.".
                             format(what_data, data_dict.keys()))
        return regr_data

    def info(self):
        """Print info for the instance of TableData object.

        For output descriptions see the method 'get_info()'.

        """
        print("File List: \n{0}".format(np.array(self._files_)))
        print("df Index Keys: \n{0}\n".format(np.array(self._df_index_keys_)))
        print("---- VERBOSE OUTPUT START ----")
        print(*self._for_info_, sep="\n")
        print("---- VERBOSE OUTPUT  END  ----")

    def get_info(self):
        """Return what info is printed in the 'info()' method.

        Returns
        -------
        files: list
            File paths where data was loaded from.
        df_index_keys: list
            Index keys added to the DataFrame object once
            multiple files are joined together such that one can
            access data by file after they were joined.
        for_info: list
            Running list of print statements that include but
            are not limited to what is shown if 'verbose=True'.
        """
        return tuple([self._files_, self._df_index_keys_, self._for_info_])

    def plot_3D_class_data(
        self,
        axes=None,
        fig_size=(4, 5),
        mark_size=12,
        which_val=0,
        save_fig=False,
        plt_str="0",
        color_list=None,
    ):
        """Plot the 3D classification data in a 2D plot.

        3 input axis with classification output.

        Parameters
        ----------
        axes: list, optional
            By default it will order the axes as [x,y,z] in the original order
            the input axis were read in. To change the ordering, pass a list
            with the column names.
            Example:
            The default orderd is col_1, col_2, col_3.
            To change the horizontal axis from col_1 to col_2 you would use:
                'axes = ["col_2", "col_1", "col_3"]'
        fig_size: tuple, optional, default = (4,5)
            Size of the figure. (Matplotlib figure kwarg 'fig_size')
        mark_size: float, optional, default = 12
            Size of the scatter plot markers. (Matplotlib scatter kwarg 's')
        which_val: int, default = 0
            Integer choosing what unique value to 'slice' on in the
            3D data such that it can be plotted on 2D.
            (If you had x,y,z data you need to choose a z value)
        save_fig: bool, default = False
            Save the figure in the local directory.
        plt_str: str, default = '0'
            If you are saving multiple figures you can pass a string
            which will be added to the end of the default:
            "data_plot_{plt_str}.pdf"
        color_list: list, default = None

        Returns
        -------
        matplotlib figure

        """
        full_data = self._full_data_
        input_data = self._input_

        if len(input_data) == 2:
            print("2D input data")
        if len(input_data) == 3:
            print("3D input data")

        if axes is not None:
            first_axis, second_axis, third_axis = axes
        else:
            first_axis, second_axis, third_axis = input_data.keys()
        print("Axes: {0}, {1}, {2}".
              format(first_axis, second_axis, third_axis))

        if color_list:
            colors = color_list
            print(
                "To set default colors for this object,\
                    re-istantiate using the option 'my_colors'."
            )
        else:
            colors = self._class_colors_

        color_dict = OrderedDict()
        for j, color_str in enumerate(colors):
            color_dict[j] = color_str

        # MAKE LEGEND
        legend_elements = []
        for i, name in enumerate(self._unique_class_keys_):
            legend_elements.append(
                Line2D(
                    [],
                    [],
                    marker="s",
                    color=color_dict[i],
                    label=name,
                    linestyle="None",
                    markersize=8,
                )
            )
        # --------------

        # Specify value of other '3rd' axis not in 2D plot
        slice_value = np.unique(full_data[third_axis])[which_val]

        # boolean data frame
        where_to_slice = full_data[third_axis] == slice_value

        # Find all indicies (rows) that have the slice value
        what_index = np.where(where_to_slice)[0]
        IDs = np.array(self._class_col_to_ids_)
        data_to_plot = full_data[where_to_slice]

        # class IDs (ints 0->n) to color code
        class_to_colors = [color_dict[val] for val in IDs[what_index]]

        fig = plt.figure(figsize=fig_size, dpi=120)
        plt.title(third_axis + "= %f" % (slice_value))
        plt.scatter(
            data_to_plot[first_axis],
            data_to_plot[second_axis],
            c=class_to_colors,
            cmap=None,
            s=mark_size,
            marker="s",
        )

        plt.xlabel(first_axis)
        plt.ylabel(second_axis)
        plt.legend(handles=legend_elements, bbox_to_anchor=(1.03, 1.02))
        if save_fig:
            plt.savefig("data_plot_{0}.pdf".format(plt_str),
                        bbox_inches="tight")
        return fig

    def make_class_data_plot(self, fig, ax, axes_keys, my_slice_vals=None,
                             my_class_colors=None, return_legend_handles=False,
                             verbose=False, **kwargs):
        """Plot classification data on a given axis and figure.

        Parameters
        ----------
        fig : Matplotlib Figure object
            Figure on which the plot will be made.
        ax : Matplotlib Axis object
            Axis on which a scatter plot will be drawn.
        axes_keys : list
            List containing two names of input data columns to use as
            horizontal and verital axes.
        my_slice_vals : optional, list, dict
            List giving values on which to slice in the axes not being plotted.
            Default (None) uses the first unique value found in each axis.
            If instead of individual values, a range is desired (e.g. 10 +/- 1)
            then a dict can be given with integer keys mapping to a tuple with
            the lower and upper range. ( e.g. {0:(9,11)} )
        my_class_colors : optional, list
            List of colors used to represent different classes. Default (None)
            uses the default class colors.
        return_legend_handles : optional, bool
            Returns a list of handles that connect classes to colors.
        verbose : optional, bool
            Print useful information during runtime.
        **kwargs : optional
            Kwargs for matplotlib.pyplot.scatter() .

        Returns
        -------
        fig : Matplotlib Figure object
            Figure object.
        ax : Matplotlib Axis object
            Updated axis after plotting.
        handles : optional, list
            List of legend handles connecting classes and their colors.
            Returned if 'return_legend_handles' is True. Default is False.

        """
        # Converting class ids to colors for plotting
        if my_class_colors is None:
            class_colors = self._class_colors_
        else:
            if not isinstance(my_class_colors, (list, np.ndarray)):
                raise ValueError("'my_class_colors' takes a list of strings")
            if verbose:
                print("Using my_class_colors...")
            class_colors = my_class_colors

        color_dict = OrderedDict()
        for j, color_str in enumerate(class_colors):
            color_dict[j] = color_str
        class_to_colors = np.array([color_dict[val]
                                    for val in self._class_col_to_ids_])

        # Legend Handles (if desired)
        # For single plot, recomended kwarg: bbox_to_anchor = (1.03, 1.02)
        legend_elements = []
        for i, name in enumerate(self._unique_class_keys_):
            legend_elements.append(
                Line2D(
                    [],
                    [],
                    marker="s",
                    color=color_dict[i],
                    label=name,
                    linestyle="None",
                    markersize=8,
                )
            )

        # Seperate the plotting and slicing axes
        all_input_axes = self._input_.keys()
        num_dim = len(all_input_axes)
        # plotting_axes = [i for i in all_input_axes if i in axes_keys]
        slicing_axes = [i for i in all_input_axes if i not in axes_keys]
        # Check that user provided valid axes
        for ax_key in axes_keys:
            if ax_key not in all_input_axes:
                raise KeyError(ax_key)

        # Slicing algorithm
        # By default slice along the first unique element in all non-plotting
        # axes. We want to plot data in a constant hyperplane. That is all
        # points with the same values in all extra dimensions.
        if num_dim >= 3:
            running_bool_list = []
            for j, slice_axis in enumerate(slicing_axes):
                unique_vals = np.unique(self._input_[slice_axis])
                if verbose:
                    print("Slice Axis: '{0}'".format(slice_axis))
                    print("\tUnique values: {0}".format(len(unique_vals)))
                if my_slice_vals is None:
                    if verbose:
                        print("\tSlice value: {0}, '{1}'".
                              format(slice_axis, unique_vals[0]))
                    where_this_val = np.array(
                        self._input_[slice_axis] == unique_vals[0]
                    )
                else:
                    if len(my_slice_vals) != len(slicing_axes):
                        raise ValueError(
                            "Need {0} slice value(s) and {1} given.".format(
                                len(slicing_axes), len(my_slice_vals)
                            )
                        )
                    if verbose:
                        print("\tSlice val: {0}, '{1}'".
                              format(slice_axis, my_slice_vals[j]))
                    if isinstance(my_slice_vals, (list, np.ndarray)):
                        where_this_val = np.array(
                            self._input_[slice_axis] == my_slice_vals[j]
                        )
                    elif isinstance(my_slice_vals, (dict, OrderedDict)):
                        where_this_val = np.array(
                            (self._input_[slice_axis] >= my_slice_vals[j][0])
                            & (self._input_[slice_axis] <= my_slice_vals[j][1])
                        )
                    else:
                        raise ValueError(
                            "`my_slice_vals` must be list or dict, {0} given.".
                            format(type(my_slice_vals)))
                running_bool_list.append(where_this_val)  # for one axis
            all_results = np.array(
                running_bool_list
            ).T  # this is now (N_points x N_slice_axes)
            # For a point (row), True if all other axes contain the slice
            # values
            rows_where_all_true = np.array([i.all() for i in all_results])

            x_data = self._input_[axes_keys[0]].loc[rows_where_all_true]
            y_data = self._input_[axes_keys[1]].loc[rows_where_all_true]
            class_to_colors = class_to_colors[rows_where_all_true]
            if verbose:
                print("\tN total points: {0}".
                      format(len(self._input_[axes_keys[0]])))
                print("\tN points in hyperplane: {0}\n".
                      format(np.sum(rows_where_all_true)))
        else:
            x_data = self._input_[axes_keys[0]]
            y_data = self._input_[axes_keys[1]]

        # The actual plotting
        ax.scatter(x_data, y_data, c=class_to_colors, **kwargs)
        ax.set_xlabel(axes_keys[0])
        ax.set_ylabel(axes_keys[1])

        if return_legend_handles:
            return fig, ax, legend_elements
        else:
            return fig, ax
