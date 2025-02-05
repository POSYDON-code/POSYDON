"""Own interpolator classes."""


__authors__ = [
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>",
]


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
        self.x = np.array(x)
        self.y = np.array(y)
        if kind not in ['linear']:
            raise NotImplementedError(f"kind = {kind} is not supported")
        # check that the data is sorted;
        # if not, check if it's strictly decreasing and flip arrays accordingly
        if not np.all(np.diff(self.x) > 0):
            if np.all(np.diff(self.x) < 0):
                self.x = self.x[::-1]
                self.y = self.y[::-1]
            else:
                raise ValueError("x values must be strictly increasing or strictly decreasing.")
        
        self.kind = kind
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


class PchipInterpolator2:
    """Interpolation class."""
    def __init__(self, *args, positive=False, **kwargs):
        """Initialize the interpolator.

        Parameters
        ----------
        positive : bool (default: False)
            If `True` request all y-values to be positive.
        *args, **kwargs : dict (optional)
            Arguments passed to original PchipInterpolator init.

        """
        self.interpolator = PchipInterpolator(*args, **kwargs)
        self.positive = positive

    def __call__(self, *args, **kwargs):
        """Use the interpolator.

        Parameters
        ----------
        *args, **kwargs : dict (optional)
            Arguments passed to original PchipInterpolator call.

        Returns
        -------
        ndarray-like
            Interpolated values.

        """
        result = self.interpolator(*args, **kwargs)
        if self.positive:
            result = np.maximum(result, 0.0)
        return result

