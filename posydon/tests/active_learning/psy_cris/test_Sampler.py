"""Unit test for posydon.active_learning.psy_cris classes
"""
import math
import unittest

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from posydon.active_learning.psy_cris.classify import Classifier
from posydon.active_learning.psy_cris.data import TableData
from posydon.active_learning.psy_cris.regress import Regressor
from posydon.active_learning.psy_cris.sample import Sampler
from posydon.active_learning.psy_cris.synthetic_data.synth_data_3D import get_output_3D
from posydon.active_learning.psy_cris.utils import (
    get_random_grid_df,
    get_regular_grid_df,
)

SKIP_TEST_PLOTS = False
SHOW_PLOTS = False


class TestSampler(unittest.TestCase):
    """Test Sampler class on the 3d synthetic data set."""

    @classmethod
    def setUpClass(cls):
        np.random.seed(12345)
        cls.TEST_DATA_GRID = get_regular_grid_df(N=10 ** 3, dim=3)
        cls.TEST_DATA_RAND = get_random_grid_df(N=10 ** 3, dim=3)
        cls.UNIQUE_CLASSES = [1, 2, 3, 4, 6, 8]

        cls.TEST_INPUT_POINTS = np.array([[0, 0, 0], [-0.5, 0.5, 0.5]])
        cls.TEST_OUTPUT_TRUTH = get_output_3D(*cls.TEST_INPUT_POINTS.T)

    def setUp(self):
        my_kwargs = {"n_neighbors": [2, 3, 5]}
        self.table_grid = self.create_TableData(self.TEST_DATA_GRID, **my_kwargs)
        self.regr_grid = Regressor(self.table_grid)
        self.cls_grid = Classifier(self.table_grid)

    def create_TableData(self, data_frame, **kwargs):
        files = None
        input_cols = ["input_1", "input_2", "input_3"]
        output_cols = ["class", "output_1"]
        class_col_name = "class"
        table_obj = TableData(
            files,
            input_cols,
            output_cols,
            class_col_name,
            my_DataFrame=data_frame,
            verbose=False,
            **kwargs
        )
        return table_obj

    def train_everything_grid(self, cls_names, regr_names):
        if cls_names is not None:
            self.cls_grid.train_everything(cls_names)
        if regr_names is not None:
            self.regr_grid.train_everything(regr_names)

    def test_init_0(self):
        test_cases = [
            (None, None),
            (self.cls_grid, self.regr_grid),
            (None, self.regr_grid),
        ]
        for i, tup_input in enumerate([(None, None), ()]):
            with self.subTest("Sampler init", iter=i):
                samp = Sampler(*tup_input)

    def test_mcmc(self):
        self.train_everything_grid(["rbf"], None)
        samp = Sampler(classifier=self.cls_grid, regressor=None)

        steps, acc, rej = samp.run_MCMC(
            15, 0.25, [0, 0, 0], samp.TD_classification, "rbf", T=1, **{"TD_BETA": 2}
        )
        self.assertTrue(len(steps) == (acc + 1), msg="steps taken should match acc.")

    def test_ptmcmc(self):
        self.train_everything_grid(["rbf"], ["rbf"])
        samp = Sampler(classifier=self.cls_grid, regressor=self.regr_grid)

        chain_step_hist, T_list = samp.run_PTMCMC(
            5,
            15,
            samp.TD_classification_regression,
            ("rbf", "rbf"),
            init_pos=[0, 0, 0],
            alpha=0.25,
            verbose=False,
            trace_plots=False,
            TD_BETA=1,
        )
        # try with default values
        chain_step_hist, T_list = samp.run_PTMCMC(
            10, 15, samp.TD_classification_regression, ("rbf", "rbf"),
            verbose=False, trace_plots=False)


    def test_simple_density_logic(self):
        self.cls_grid.train("rbf")
        samp = Sampler(classifier=self.cls_grid, regressor=None)
        steps, acc, rej = samp.run_MCMC(
            200, 0.25, [0, 0, 0], samp.TD_classification, "rbf", T=1
        )
        acc_pts, rej_pts = samp.do_simple_density_logic(steps, 10, 0.05)
        return samp, steps

    def test_get_proposed_points(self):
        N = 10
        samp, step_hist = self.test_simple_density_logic()
        prop_points, kappa = samp.get_proposed_points(step_hist, N, 0.046)
        self.assertTrue(len(prop_points) == N)

    @unittest.skipIf(SKIP_TEST_PLOTS, "Plotting C, C+R target distributions")
    def test_TD_plots(self):
        self.train_everything_grid(["rbf"], ["rbf"])
        samp = Sampler(classifier=self.cls_grid, regressor=self.regr_grid)

        N = 70 if SHOW_PLOTS else 5
        zed = 0
        x, y = np.meshgrid(np.linspace(-1, 1, N), np.linspace(-1, 1, N))
        z = np.ones(x.shape) * zed
        data_points = np.concatenate(
            (x.flatten()[:, None], y.flatten()[:, None], z.flatten()[:, None]), axis=1
        )

        max_probs, pos, cls_keys = samp.get_TD_classification_data("rbf", data_points)

        kwargs = {"TD_BETA": 2, "TD_TAU": 0.5}
        cls_regr_td_vals = [
            float(samp.TD_classification_regression(["rbf", "rbf"], dat, **kwargs))
            for dat in data_points
        ]

        fig, subs = plt.subplots(1, 2, figsize=(13, 5))
        subs[0].set_title("TD_classification at z = {}".format(zed))
        cls_plot = subs[0].pcolormesh(
            x, y, (1 - max_probs).reshape(N, N), shading="auto"
        )

        subs[1].set_title("TD_classification_regression at z = {}".format(zed))
        cls_regr_plot = subs[1].pcolormesh(
            x, y, np.array(cls_regr_td_vals).reshape(N, N), shading="auto"
        )

        fig.colorbar(cls_plot, ax=subs[0])
        fig.colorbar(cls_regr_plot, ax=subs[1])

        for i in range(2):
            subs[i].set_xlabel("input_1")
            subs[1].set_ylabel("input_2")
            subs[i].axis("equal")

        if SHOW_PLOTS:
            plt.show()
        plt.close()


if __name__ == "__main__":
    unittest.main()
