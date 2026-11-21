"""Unit tests of the single-point classifier fast path in
posydon/interpolation/IF_interpolation.py
"""

__authors__ = [
    "Seth Gossage <seth.gossage@northwestern.edu>"
]

from unittest import mock

import numpy as np
from pytest import fixture
from sklearn.neighbors import KNeighborsClassifier

import posydon.interpolation.IF_interpolation as totest

LABELS = np.array(["stable_MT", "unstable_MT", "no_MT", "initial_MT",
                   "not_converged"], dtype=object)


def _fit(X, y, k, **kw):
    return KNeighborsClassifier(n_neighbors=k, weights="distance", **kw).fit(X, y)


class TestPredictSingleDistanceWeighted:

    @fixture
    def rng(self):
        return np.random.default_rng(0)

    def test_matches_sklearn(self, rng):
        n_checked = 0
        for trial in range(20):
            n = int(rng.integers(30, 800))
            k = int(rng.integers(1, 15))
            X = rng.uniform(size=(n, 3))
            if trial % 3 == 0:          # gridded inputs -> distance ties
                X = np.round(X * 5) / 5
            y = LABELS[rng.integers(0, int(rng.integers(2, 6)), size=n)]
            clf = _fit(X, y, k)
            queries = np.vstack([rng.uniform(-0.2, 1.2, size=(60, 3)),
                                 X[:20],                       # zero distance
                                 np.round(rng.uniform(size=(20, 3)) * 10) / 10])
            for q in queries:
                q = q[None, :]
                fast = totest._predict_single_distance_weighted(clf, q)
                ref = clf.predict(q)
                assert fast is not None
                assert fast.dtype == ref.dtype and fast.shape == ref.shape
                assert (fast == ref).all()
                n_checked += 1
        assert n_checked == 20 * 100

    def test_not_applicable(self, rng):
        X = rng.uniform(size=(50, 3))
        y = LABELS[rng.integers(0, 3, size=50)]
        clf = _fit(X, y, 4)
        q = X[:1] + 0.01
        f = totest._predict_single_distance_weighted
        assert f(clf, np.vstack([q, q])) is None                    # 2 rows
        assert f(clf, q[0]) is None                                 # 1-D
        assert f(clf, q.tolist()) is None                           # not array
        assert f(clf, np.array([[np.nan, 0.5, 0.5]])) is None       # non-finite
        assert f(clf, q[:, :2]) is None                             # wrong dim
        uniform = KNeighborsClassifier(n_neighbors=4).fit(X, y)
        assert f(uniform, q) is None                                # weights
        brute = _fit(X, y, 4, algorithm="brute")
        assert f(brute, q) is None                                  # no tree
        y2 = np.stack([y, y[::-1]], axis=1)
        assert f(_fit(X, y2, 4), q) is None                         # 2-D y
        too_many = _fit(X[:3], y[:3], 3)
        too_many.n_neighbors = 5
        assert f(too_many, q) is None                               # k > n

    def test_knnclassifier_uses_fast_path(self, rng):
        X = rng.uniform(size=(100, 3))
        y = LABELS[rng.integers(0, 4, size=100)]
        c = totest.KNNClassifier()
        c.train(X, y, K=5)
        q = rng.uniform(size=(1, 3))
        with mock.patch.object(c.classifier, "predict",
                               wraps=c.classifier.predict) as p:
            out = c.predict(q)
            assert p.call_count == 0
            assert (out == _fit(X, y, 5).predict(q)).all()
            c.predict(np.vstack([q, q]))                # falls back
            assert p.call_count == 1


class TestEvaluateMatClassifiesOnce:

    def test_c_key_classes_are_reused(self):
        interp = object.__new__(totest.BaseIFInterpolator)
        interp.n_in = 3
        interp.c_key = "interpolation_class"
        interp.class_method = "kNN"
        interp.interp_method = "linear"
        ycat = {"interpolation_class": np.array(["stable_MT"], dtype=object),
                "S1_state": np.array(["BH"], dtype=object)}
        interp.test_classifiers = mock.Mock(return_value=ycat)
        interp.test_interpolator = mock.Mock(return_value=np.ones((1, 2)))
        Xt = np.array([10.0, 5.0, 3.0])
        ynum, ycat_out = interp.evaluate_mat(Xt)
        interp.test_classifiers.assert_called_once()
        args, kwargs = interp.test_interpolator.call_args
        assert kwargs["classes"] is ycat["interpolation_class"]
        assert ycat_out is ycat
