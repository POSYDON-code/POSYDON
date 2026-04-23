"""Unit test for posydon.active_learning.psy_cris classes
"""
import math
import unittest

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from posydon.active_learning.psy_cris.data import TableData
from posydon.active_learning.psy_cris.regress import Regressor
from posydon.active_learning.psy_cris.synthetic_data.synth_data_3D import get_output_3D
from posydon.active_learning.psy_cris.utils import (
    get_random_grid_df,
    get_regular_grid_df,
)

# True for faster runtime ~ 1s vs 3s
SKIP_GP_TESTS = False

SKIP_TEST_PLOTS = False
SHOW_PLOTS = False


class TestRegressor(unittest.TestCase):
    """Test Regressor class on the 3d synthetic data set."""

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
        self.table_rand = self.create_TableData(self.TEST_DATA_RAND, **my_kwargs)
        self.regr_grid = Regressor(self.table_grid)
        self.regr_rand = Regressor(self.table_rand)

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

    def train_grid_regressors(self, names, *args, **kwargs):
        for i, regr_name in enumerate(names):
            self.regr_grid.train(regr_name, *args, **kwargs)

    def test_train_regressors_1(self):
        # di can be [16,19,24,32,41,48,64] to avoid multi class error for gp
        classes_to_train = [ [1,2,3,4,6,8], [1,6,8], [1,6,8] ]
        col_keys = ["output_1"]
        for i, regr in enumerate(["rbf", "linear", "gp"]):
            with self.subTest("Train grid:", regressor=regr):
                self.regr_grid.train(regr, classes_to_train[i], col_keys,
                                     verbose=False)

    def test_train_regressors_2(self):
        # skipping gp here
        classes_to_train = [ [1,2,3,4,6,8], [1,6,8], [1,6,8] ]
        col_keys = ["output_1"]
        for i, regr in enumerate(["rbf", "linear"]):
            with self.subTest("Train random:", regressor=regr):
                self.regr_rand.train(regr, classes_to_train[i],
                                     col_keys, verbose=False)

    def test_predictions(self):
        """Checking for consistency only, not true values."""
        self.train_grid_regressors(["rbf", "linear"], [6], ["output_1"])

        regr_out = self.regr_grid.get_predictions(
            ["rbf", "linear"], [6], ["output_1"], self.TEST_INPUT_POINTS
        )
        # RBF check
        for i, corr_ans in enumerate([-0.78295781, -0.7401497]):
            with self.subTest("RBF regr", correct_ans=corr_ans):
                pred = regr_out["RBF"][6]["output_1"][i]
                self.assertAlmostEqual(pred, corr_ans, places=5)

        # LinearNDInterpolator check
        for i, corr_ans in enumerate([0.13898644, float("Nan")]):
            with self.subTest("LinearNDInterpolator regr", correct_ans=corr_ans):
                pred = regr_out["LinearNDInterpolator"][6]["output_1"][i]
                if i == 1:
                    self.assertTrue(math.isnan(pred), msg="Prediction should be Nan. {}".format(pred))
                else:
                    self.assertAlmostEqual(pred, corr_ans, places=5)

    @unittest.skipIf(SKIP_GP_TESTS, "GP train / predict - longer runtime.")
    def test_predictions_gp(self):
        self.train_grid_regressors(["gp"], [6], ["output_1"], di=None)

        regr_out = self.regr_grid.get_predictions(
            ["gp"], [6], ["output_1"], self.TEST_INPUT_POINTS
        )
        for i, corr_ans in enumerate([0,0]):
            with self.subTest("GaussianProcessRegressor", correct_ans=corr_ans):
                pred = regr_out["GaussianProcessRegressor"][6]["output_1"][i]
                self.assertAlmostEqual(pred, corr_ans, places=5)


    def test_pred_train_err(self):
        # Trying to predict without training
        names = ["grid", "random"]
        for i, regressor in enumerate([self.regr_grid, self.regr_rand]):
            with self.subTest(classifier_name=names[i]):
                with self.assertRaisesRegex(
                    Exception, "No trained interpolators exist"
                ):
                    regressor.get_predictions( ["linear"], [6], ["output_1"], [[0, 0, 0]])

    def test_cross_val(self):
        corr_ans = [-16.528141760898365, -61.41327730988214, -10.621903342287425]
        for index, cls in enumerate([1,6,8]):
            with self.subTest("Cross Validation Regression", class_key=cls):
                perc_diffs, actual_diffs = self.regr_grid.cross_validate("rbf", cls, "output_1", 0.5 )
                self.assertAlmostEqual( np.mean(perc_diffs), corr_ans[index], places=5)

                plt.hist(perc_diffs, bins=40, density=True, range=(-300,300),
                         histtype="step", label="class "+str(cls))
        plt.xlabel("Regression Percent Difference")
        plt.title("Test Regression CV")
        plt.legend()
        if SHOW_PLOTS:
            plt.show()
        plt.close()


    @unittest.skipIf(SKIP_GP_TESTS, "GP cross_val - longer runtime.")
    def test_cross_val_gp(self):
        corr_ans = [-36.08266156603506, -100.0, -70.79557229587994]
        for index, cls in enumerate([1,6,8]):
            with self.subTest("Cross Validation Regression GP", class_key=cls, ans=corr_ans[index]):
                perc_diffs, actual_diffs = self.regr_grid.cross_validate("gp", cls, "output_1", 0.5 )
                self.assertAlmostEqual( np.mean(perc_diffs), corr_ans[index], places=5)

                plt.hist(perc_diffs, bins=40, density=True, range=(-300,300),
                         histtype="step", label="class "+str(cls))
        plt.xlabel("Regression Percent Difference")
        plt.title("Test GP Regression CV")
        plt.legend()
        if SHOW_PLOTS:
            plt.show()
        plt.close()

    @unittest.skipIf(SKIP_TEST_PLOTS, "All regression data plot.")
    def test_max_cls_plot(self):
        class_key = 1
        fig = self.regr_grid.plot_regr_data(class_key)
        if SHOW_PLOTS:
            plt.show()
        else:
            print("To show plots set SHOW_PLOTS to True.")
        plt.close()


if __name__ == "__main__":
    unittest.main()
