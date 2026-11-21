"""Unit tests of helpers in posydon/binary_evol/DT/step_detached.py"""

__authors__ = [
    "Seth Gossage <seth.gossage@northwestern.edu>"
]

import numpy as np

import posydon.binary_evol.DT.step_detached as totest


def test_non_negative_matches_np_max():
    values = [0.0, -0.0, 1e-300, -1e-300, 5e-324, -5e-324, 1.0, -1.0,
              np.inf, -np.inf, np.nan, 1e308, -1e308]
    for x in values:
        v = np.float64(x)
        ref = np.max([v, 0])
        got = np.float64(totest._non_negative(v))
        if np.isnan(ref):
            assert np.isnan(got)
        else:
            assert got == ref and np.signbit(got) == np.signbit(ref)
