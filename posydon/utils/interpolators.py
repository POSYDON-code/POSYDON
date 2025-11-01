"""Own interpolator classes."""


__authors__ = [
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>",
    "Seth Gossage <seth.gossage@northwestern.edu>"
]


import copy

import numpy as np
from scipy.interpolate import PchipInterpolator


class interp1d:
    """Interpolation class for one dimensional interpolation."""
    def __init__(self, x, y, kind='linear', **kwargs):
        """Initialize the interpolator.

        Parameters
        ----------
        x : ndarray-like
            X-coordinates of the data to interpolate.
        y : ndarray-like
            Y-coordinates of the data to interpolate. Needs to have the same
            length as `x`.
        kind : str (default: 'linear')
            Currently only 'linear' is allowed.
                - 'linear' : uses np.interp
        **kwargs : dict (optional)
            Dictionary containing extra parameters.
                - 'left' : float
                    Y-value to return for evaluations below the data `x`.
                - 'right' : float
                    Y-value to return for evaluations above the data `x`.
                - 'fill_value' : 'extrapolate' or tuple of size 1 or 2
                    Superseded by 'left' or 'right'
                    - 'extrapolate' : use the first and last elements in the
                         data to determine the extrapolation below and above,
                         respectively; ignores fixed values for below and above
                    - size 1 : use y-value for below and above
                    - size 2 : use first y-value for below and second for above

        """
        if kind not in ['linear']:
            raise NotImplementedError(f"kind = {kind} is not supported")
        self.kind = kind
        self.x = np.array(x)
        self.y = np.array(y)
        # check that x is increasing
        if not np.all(np.diff(self.x) > 0):
            # if instead being strictly decreasing, flip data
            if np.all(np.diff(self.x) < 0):
                self.x = np.flip(self.x)
                self.y = np.flip(self.y)
            else:
                # make them strictly increasing if neither
                indices = np.argsort(self.x)
                self.x = self.x[indices]
                self.y = self.y[indices]

        self.below = None
        self.above = None
        self.extrapolate = False
        if 'fill_value' in kwargs:
            if kwargs['fill_value'] == 'extrapolate':
                self.extrapolate = True
            elif (isinstance(kwargs['fill_value'], tuple) and
                (len(kwargs['fill_value']) == 1)):
                self.below = kwargs['fill_value'][0]
                self.above = kwargs['fill_value'][0]
            elif  (isinstance(kwargs['fill_value'], tuple) and
                   (len(kwargs['fill_value']) == 2)):
                self.below = kwargs['fill_value'][0]
                self.above = kwargs['fill_value'][1]
            else:
                raise NotImplementedError("fill_value has to be 'extrapolate' "
                                          "or a tuple with 1 or 2 elements")
        if 'left' in kwargs:
            self.below = kwargs['left']
        if 'right' in kwargs:
            self.above = kwargs['right']
        # check that the interpolator can be called
        assert self.__call__(x[0]) == y[0]

    def __call__(self, x_new):
        """Use the interpolator.

        Parameters
        ----------
        x_new : float or ndarray
            X-coordinate(s) to get interpolated y-value(s) for.

        Returns
        -------
        ndarray
            Interpolated y-values.

        """

        x_n = np.array(x_new)
        if self.kind == 'linear':
            y_interp = np.array(np.interp(x=x_n, xp=self.x, fp=self.y,
                                          left=self.below, right=self.above))
            if (self.extrapolate and (len(self.x)>1)):
                below_mask = x_n < self.x[0]
                above_mask = x_n > self.x[-1]
                if np.any(below_mask):
                    slope_below = ((self.y[1]-self.y[0])
                                   / (self.x[1]-self.x[0]))
                    y_interp[below_mask] = self.y[0] + slope_below * (
                                            x_n[below_mask]-self.x[0])
                if np.any(above_mask):
                    slope_above = ((self.y[-1]-self.y[-2])
                                   / (self.x[-1]-self.x[-2]))
                    y_interp[above_mask] = self.y[-1] + slope_above * (
                                            x_n[above_mask]-self.x[-1])
            return y_interp
        else:
            raise NotImplementedError(f"kind = {self.kind} is not supported")


class StellarInterpolator:
    """Interpolation class for stellar evolution tracks."""

    def __init__(self, t, y_data, y_keys, positives=None, derivatives=None):
        """Initialize the interpolator.

        Parameters
        ----------
        t : 1D-ndarray
            Time coordinates of the data to interpolate.
        y_data : list of 1D-ndarray
            List of Y-coordinates of the data to interpolate. Each entry
            needs to have the same length as `t`.
        y_leys : list of str
            List of keys corresponding to each y_data entry.
        positives : list of bool (default: None)
            List indicating for each y_data entry whether to request positive
            y-values. If `None`, all entries are assumed to be non-positive.
        derivatives : list of bool (default: None)
            List indicating for each y_data entry whether to return the
            derivative of the interpolator. If `None`, all entries are assumed
            to return the interpolated y-values.
            """

        # get all positive flags
        if positives is None:
            positives = [False for _ in y_data]
        if derivatives is None:
            derivatives = [False for _ in y_data]

        # Create interpolators for each combination
        self._interpolators = {}
        self._keys = {}
        for is_positive in [True, False]:
            for is_derivative in [True, False]:
                # Filter data for this combination
                mask = [(p == is_positive and d == is_derivative)
                        for p, d in zip(positives, derivatives)]
                filtered_data = [y for i, y in enumerate(y_data) if mask[i]]
                filtered_keys = [y_keys[i] for i in range(len(y_keys)) if mask[i]]

                if not filtered_data:
                    continue

                # Create interpolator only if there is data for this combination
                key = (is_positive, is_derivative)
                if filtered_data:
                    base_interp = PchipInterpolator(t, np.array(filtered_data).T)
                    self._keys[key] = filtered_keys
                    if is_derivative:
                        self._interpolators[key] = base_interp.derivative()
                    else:
                        self._interpolators[key] = base_interp

        self.t_max = t.max()
        self.max_time = np.inf
        self.t0 = 0.0
        self.m0 = 0.0

        # potential offset (for example, we've matched a MESA sim to another
        # but need to offset the matched track's age to fast forward to current
        # age in the evolution)
        self.offset = 0.0
        self._positives = positives
        self._derivatives = derivatives


    @property
    def keys(self):
        """List of all keys in the interpolator."""
        return [key for key_list in self._keys.values() for key in key_list]


    def __call__(self, t):
        """Use the interpolator.

        Parameters
        ----------
        t : float or ndarray
            Time coordinate(s) to get interpolated y-value(s) for.

        Returns
        -------
        list of ndarray-like
            List of interpolated y-values for each y_data entry.

        """
        results = {}

        # offset t

        t = t - self.offset

        for (is_positive, is_derivative), interp in self._interpolators.items():
            # Evaluate interpolator
            values = interp(t)

            # Apply positive constraint if needed
            if is_positive:
                values = np.maximum(values, 0.0)

            # Map values to keys
            keys = self._keys[(is_positive, is_derivative)]

            if np.ndim(values) == 1:
                for i, key in enumerate(keys):
                    results[key] = values[i]
            else:
                for i, key in enumerate(keys):
                    results[key] = values[:, i]

        return results
