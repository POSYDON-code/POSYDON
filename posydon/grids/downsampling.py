"""Module for downsampling accross a path in multi-dimensional space.

Let a d-dimensional parametric curve in of the form:

    x_1 = x_1(t)
    x_2 = x_2(t)
    ...
    x_d = x_d(t)

sampled at N time steps, t_1, t_2, ..., t_N. It is possible to select a subset
of the time steps - essentially defining a downsampled version of the curve -
so that the original curve can be reconstructed with minimum error.

The provided class `TrackDownsampler` performs the downsampling. It requires
two input arrays. The first is the "independent" variable, t, which does not
have to be temporal (e.g., spatial), but must be strictly increasing. The 2nd
array is the track, in which the different columns correspond to the different
dimensions:

    --- parameters --->

  |     x_1(t_1)    x_2(t_1)    ...     x_d(t_1)
time    x_1(t_2)    x_2(t_2)    ...     x_d(t_2)
 or         .           .       ...         .
 any        .           .       ...         .
  |         .           .       ...         .
  V     x_1(t_N)    x_2(t_N)    ...     x_d(t_N)


 Usage example
 -------------
 TD = TrackDownsample(t, X)
 t_new, X_new = TD.downsample(max_err=0.001)


"""


__authors__ = [
    "Konstantinos Kovlakas <Konstantinos.Kovlakas@unige.ch>",
]


import numpy as np
import warnings
import sys
from posydon.utils.posydonwarning import InappropriateValueWarning


sys.setrecursionlimit(100000)


def rescale_from_0_to_1(array, return_minmax=False):
    """Rescale a matrix so that all dimensions span from 0 to 1."""
    arr = np.array(array)
    if len(array.shape) != 2:
        raise ValueError("Scaler works for matrices only.")

    minima, maxima = np.min(arr, axis=0), np.max(arr, axis=0)

    for col, minimum, maximum in zip(range(arr.shape[1]), minima, maxima):
        if minimum == maximum:
            arr[:, col] = 0.0
        else:
            arr[:, col] = (arr[:, col] - minimum) / (maximum - minimum)

    if return_minmax:
        return arr, minima, maxima
    return arr


class TrackDownsampler:
    """Class performing downsampling of multi-dimensional paths."""

    def __init__(self, independent, dependent, verbose=False):
        """Get, reduce and rescale data."""
        self.verbose = verbose
        self.say("Initializing downsampler.")
        self.keep = None    # boolean array inidicating rows to keep

        if not np.iterable(independent):
            raise ValueError("`independent` is not iterable.")
        if np.iterable(independent[0]):
            raise ValueError("`independent` should be one-dimensional.")
        self.t = np.array(independent)

        if not np.iterable(dependent):
            raise ValueError("`dependent` is not iterable.")

        self.X = np.array(dependent)

        if not np.iterable(dependent[0]):
            self.X = self.X.reshape((-1, 1))

        if np.iterable(self.X[0, 0]):
            raise ValueError("Number of dimensions in `dependent` > 2.")

        # ensure `independent` is strictly increasing
        if np.any(np.diff(self.t) <= 0):
            raise ValueError("`independent` should be strictly increasing.")

        self.rescale()

    def say(self, message):
        """Speak to the standard output, only if `.verbose` is True."""
        if self.verbose:
            print(message)

    def rescale(self):
        """Rescale the parameter space from 0 to 1."""
        self.say("Rescaling data.")
        self.X, self.minima, self.maxima = rescale_from_0_to_1(
            self.X, return_minmax=True)
        self.say("    done.")

    def extract_downsample(self, scale_back=True):
        """Extract the downsampled array."""
        t, X = self.t[self.keep], self.X[self.keep, :]
        if scale_back:
            t, X = t, (X * (self.maxima - self.minima)) + self.minima
        return t, X

    def find_downsample(self, max_err=None, max_interval=None):
        """Find the rows of the data that constitute the "down-sample".

        Note: if `max_interval` is negative, it's value is the relative ratio
              between maximum allowed dm and total star mass.
        """
        t, X = self.t, self.X
        N = len(t)

        if max_err is None or N < 3:
            return np.ones_like(t, dtype=bool)
        if max_err < 0.0 or not np.isfinite(max_err):
            raise ValueError("`max_err` must be a non-negative finite number.")
        if max_err >= 1.0:
            warnings.warn("`max_err` >= 1.0 is like disabling downsampling.",
                          InappropriateValueWarning)

        keep = np.zeros_like(t, dtype=bool)  # initially keep no row

        def DS(i, j, tolerance, max_interval=None):
            if j <= i + 1:
                return
            X_i, X_j = X[i, :], X[j, :]    # points at the ends

            t_i, dt = t[i], t[j] - t[i]    # t-coordinate at the start and dt
            t_intermediate = t[i+1:j]      # t values for intermediate points
            X_intermediate = X[i+1:j]      # intermediate points

            # interpolate all coordinates simulataneously
            X_intepolated = X_i + np.outer(t_intermediate-t_i, X_j-X_i) / dt
            del X_j, X_i, t_intermediate, t_i, dt

            # compute the error (as in 2-norm) for all intermediate points
            X_err = np.linalg.norm(X_intermediate - X_intepolated, axis=1)
            del X_intermediate, X_intepolated
            # and find the maximum error and the corresponding index
            argmax = np.argmax(X_err)
            e_max = X_err[argmax]
            k = i + 1 + argmax
            del X_err, argmax

            # if below tolerance search no more,
            # otherwise mark the point and continue at its left and right
            keep_it = False
            if max_interval is not None:
                interval = abs(t[i] - t[j])
                if max_interval < 0:
                    keep_it = interval > abs(max_interval * (t[0] - t[-1]))
                else:
                    keep_it = interval > max_interval
            if e_max > tolerance:
                keep_it = True

            if keep_it:
                keep[k] = True
                DS(i, k, tolerance, max_interval=max_interval)
                DS(k, j, tolerance, max_interval=max_interval)

        # keep the first and last point
        keep[0], keep[N-1] = True, True
        # ...and the intermediate ones that can't be accurately interpolated
        DS(0, N-1, tolerance=max_err, max_interval=max_interval)

        self.keep = keep

    def downsample(self, max_err=None, scale_back=True, max_interval=None):
        """Perform the downsampling and return the result."""
        self.find_downsample(max_err=max_err, max_interval=max_interval)
        return self.extract_downsample(scale_back=scale_back)
