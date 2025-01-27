"""Own interpolator classes."""


__authors__ = [
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>",
]


import numpy as np
from scipy.interpolate import PchipInterpolator


class interp1d:
    def __init__(self, x, y, kind='linear'):
        self.x = x
        self.y = y
        if kind not in ['linear']:
            raise NotImplementedError(f"kind = {kind} is not supported")
        self.kind = kind
        assert self.__call__(x[0]) == y[0]

    def __call__(self, x_new):
        if self.kind == 'linear':
            return np.interp(x_new, self.x, self.y)
        else:
            raise NotImplementedError(f"kind = {self.kind} is not supported")


class PchipInterpolator2:
    """Interpolation class."""

    def __init__(self, *args, positive=False, **kwargs):
        """Initialize the interpolator."""
        self.interpolator = PchipInterpolator(*args, **kwargs)
        self.positive = positive

    def __call__(self, *args, **kwargs):
        """Use the interpolator."""
        result = self.interpolator(*args, **kwargs)
        if self.positive:
            result = np.maximum(result, 0.0)
        return result

