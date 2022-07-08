"""Module for scaling data before and after interpolation."""


__authors__ = [
    "Juanga Serra Perez <jgserra@northwestern.edu>",
]


import numpy as np


class DataScaler:
    """Data Normalization class.

    This class provides normalizing tools for float 1D arrays. Features can be
    standarized or scaled to a range depending on the method chosen when
    calling the fit ot fit_and_transform functions.
    """

    def __init__(self):
        """Initialize the data scaler.

        No parameters are expected. After instantiation use methods `fit` or
        `fit_and_transform` first to fit a scaling object to a given vector of
        values.

        Example
        -------
        >>> sc = DataScaler()

        """
        self.params = None
        self.method = None
        # These two parameters will only be used when method is (log)_min_max
        self.lower = None
        self.upper = None

    def fit(self, x, method='none', lower=-1.0, upper=1.0):
        """Fit a transform of 1D numpy array x.

        Computes the parameters that define the transform.

        Parameters
        ----------
        x : numpy.ndarray
            expects a 1D array and finds norm values of columns

        method : str
            Scaling method. Possible values: 'min_max', 'max_abs', 'standarize'
            and their log versions 'log_min_max', 'neg_log_min_max',
            'log_max_abs', 'log_standarize', 'neg_log_standarize'.

        lower : float
            lower range value of x_t after (log_)min_max scaling
        upper : float
            upper range value of x_t after (log_)min_max scaling

        """
        # make sure that we have a 1D NumPy arrays
        assert isinstance(x, np.ndarray) and len(x.shape) == 1
        # this is not rechecked in transform/inv_transform,
        # so the same transform could potentially be used for ndarrays

        # store method as class attribute
        self.method = method

        # compute scaling parameters associated with each method
        if method == 'min_max':
            assert upper > lower, "upper must be greater than lower"
            self.lower, self.upper = lower, upper
            self.params = [x.min(axis=0), x.max(axis=0)]
        elif method == 'log_min_max':
            assert upper > lower, "upper must be greater than lower"
            self.lower, self.upper = lower, upper
            self.params = [np.log10(x.min(axis=0)), np.log10(x.max(axis=0))]
        elif method == 'neg_log_min_max':
            assert upper > lower, "upper must be greater than lower"
            self.lower, self.upper = lower, upper
            self.params = [np.log10((-x).min(axis=0)),
                           np.log10((-x).max(axis=0))]
        elif method == 'max_abs':
            self.params = [np.abs(x).max(axis=0)]
        elif method == 'log_max_abs':
            self.params = [np.abs(np.log10(x)).max(axis=0)]
        elif method == 'standarize':
            self.params = [x.mean(axis=0), x.std(axis=0)]
        elif method == 'log_standarize':
            # log will be computed in transform again
            self.params = [np.log10(x).mean(axis=0), np.log10(x).std(axis=0)]
        elif method == 'neg_log_standarize':    # log(-x)
            self.params = [np.log10(-x).mean(axis=0), np.log10(-x).std(axis=0)]
        elif method == 'log':
            self.params = []
        elif method == 'none':  # no transformation
            self.params = []
        else:
            raise ValueError(f"Unknown method `{method}` for data scaler.")

    def transform(self, x):
        """Transform x using the already obtained normalization values.

        `self.fit()`` must be called first.
        lower/upper will only be taken into account for (log_)min_max
        normalization. In this case, the transformed x will have
        min(x_transf) = lower, max(x_transf) = upper

        Parameters
        ----------
        x : numpy.ndarray
            values to normalize

        Returns
        -------
        numpy.ndarray
            transformed version of x

        """
        # Make sure that "self.fit" was called first
        if self.method is None:
            raise AssertionError(
                "You have to fit a scaling object to a feature vector first")

        if self.method == 'min_max':
            x_t = ((x - self.params[0]) / (self.params[1] - self.params[0])
                   * (self.upper - self.lower) + self.lower)
        elif self.method == 'log_min_max':
            x_t = ((np.log10(x) - self.params[0])
                   / (self.params[1] - self.params[0])
                   * (self.upper - self.lower) + self.lower)
        elif self.method == 'neg_log_min_max':
            x_t = ((np.log10(-x) - self.params[0])
                   / (self.params[1] - self.params[0])
                   * (self.upper - self.lower) + self.lower)
        elif self.method == 'max_abs':
            x_t = x / self.params[0]
        elif self.method == 'log_max_abs':
            x_t = np.log10(x) / self.params[0]
        elif self.method == 'standarize':
            x_t = (x - self.params[0]) / self.params[1]
        elif self.method == 'log_standarize':
            # log will be computed in transform again
            x_t = (np.log10(x) - self.params[0]) / self.params[1]
        elif self.method == 'neg_log_standarize':
            x_t = (np.log10(-x) - self.params[0]) / self.params[1]
        elif self.method == 'log':
            x_t = np.log10(x)
        else:  # no transformation
            x_t = x

        return x_t

    def fit_and_transform(self, x, method='none', lower=-1, upper=1):
        """Fit and transform the array x according to the chosen scaling.

        lower/upper will only be taken into account for (log_)min_max
        normalization. In this case, the transformed x will have
        min(x_transf) = lower, max(x_transf) = upper

        Parameters
        ----------
        x : numpy.ndarray
            expects a 1D array and finds norm values of columns
        method : str
            scaling method. Possible values: 'min_max', 'max_abs', 'standarize'
            and the log versions 'log_min_max', 'log_max_abs', 'log_standarize'

        lower : float
            lower range value of x_t after (log_)min_max scaling
        upper : float
            upper range value of x_t after (log_)min_max scaling

        Returns
        -------
        numpy.ndarray
            transformed version of x

        """
        self.fit(x, method, lower, upper)
        return self.transform(x)

    def inv_transform(self, x_t):
        """Revert the scaling using the stored transform parameters.

        Parameters
        ----------
        x_t : numpy.ndarray
            expects a 1D array to unnormalize given the fitted transform.

        Returns
        -------
        numpy.ndarray
            denormalized x using the stored parameters.

        """
        if self.method is None:
            raise AssertionError("Transformation not defined yet. "
                                 "Fit a scaling object first.")

        if self.method == 'min_max':
            x = ((x_t - self.lower)
                 / (self.upper - self.lower)
                 * (self.params[1] - self.params[0]) + self.params[0])
        elif self.method == 'log_min_max':
            x = 10 ** ((x_t - self.lower) / (self.upper - self.lower)
                       * (self.params[1] - self.params[0]) + self.params[0])
        elif self.method == 'neg_log_min_max':
            x = -10 ** ((x_t - self.lower) / (self.upper - self.lower)
                        * (self.params[1] - self.params[0]) + self.params[0])
        elif self.method == 'max_abs':
            x = x_t * self.params[0]
        elif self.method == 'log_max_abs':
            x = 10 ** (x_t * self.params[0])
        elif self.method == 'standarize':
            x = x_t * self.params[1] + self.params[0]
        elif self.method == 'log_standarize':
            x = 10 ** (x_t * self.params[1] + self.params[0])
        elif self.method == 'neg_log_standarize':
            x = -10 ** (x_t * self.params[1] + self.params[0])
        elif self.method == 'log':
            x = 10 ** x_t
        else:  # no transformation
            x = x_t

        return x
