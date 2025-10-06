"""Unit test for posydon.active_learning.psy_cris classes
"""
import unittest

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from posydon.active_learning.psy_cris.classify import Classifier
from posydon.active_learning.psy_cris.data import TableData
from posydon.active_learning.psy_cris.synthetic_data.synth_data_3D import get_output_3D
from posydon.active_learning.psy_cris.utils import (
    get_random_grid_df,
    get_regular_grid_df,
)

# True for faster runtime ~ 3s vs 15s
SKIP_GP_TESTS = True

SKIP_TEST_PLOTS = False
SHOW_PLOTS = False


class TestClassifier(unittest.TestCase):
    """Test Classifier class on the 3d synthetic data set."""

    @classmethod
    def setUpClass(cls):
        np.random.seed(12345)
        cls.TEST_DATA_GRID = get_regular_grid_df(N=10 ** 3, dim=3)
        cls.TEST_DATA_RAND = get_random_grid_df(N=10 ** 3, dim=3)
        cls.UNIQUE_CLASSES = [1, 2, 3, 4, 6, 8]

        cls.TEST_INPUT_POINTS = np.array([[0, 0, 0], [-0.5, 0.5, 0.5]])
        cls.TEST_OUTPUT_TRUTH = get_output_3D(*cls.TEST_INPUT_POINTS.T)

    def setUp(self):
        my_kwargs = {"n_neighbors": [2, 3]}
        self.table_grid = self.create_TableData(self.TEST_DATA_GRID, **my_kwargs)
        self.table_rand = self.create_TableData(self.TEST_DATA_RAND, **my_kwargs)
        self.cls_obj_grid = Classifier(self.table_grid)
        self.cls_obj_rand = Classifier(self.table_rand)

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

    def test_train_classifiers_1(self):
        # di can be [16,19,24,32,41,48,64] to avoid multi class error for gp
        data_range = np.arange(0, 1000)[::64]
        for cls in ["rbf", "linear", "gp"]:
            with self.subTest("Train grid:", classifier=cls):
                self.cls_obj_grid.train(cls, di=data_range, verbose=False)

    def test_train_classifiers_2(self):
        # skipping gp here
        data_range = np.arange(0, 1000)[::10]
        for cls in ["rbf", "linear"]:
            with self.subTest("Train random:", classifier=cls):
                self.cls_obj_rand.train(cls, di=data_range, verbose=False)

    def train_grid_classifiers(self, cls_names, **kwargs):
        for name in cls_names:
            self.cls_obj_grid.train(name, **kwargs)

    def test_predictions(self):
        self.train_grid_classifiers(["linear", "rbf"])
        correct_probabilities = [0.98519722, 1.00000, 0.5883452]
        for i, cls_name in enumerate(["rbf", "linear"]):
            with self.subTest("Get class predictions:", classifier=cls_name):
                tup_out = self.cls_obj_grid.get_class_predictions(
                    cls_name, self.TEST_INPUT_POINTS, return_ids=False
                )
                class_pred, probs, where_not_nan = tup_out

                self.assertTrue(
                    (3 in class_pred) and (8 in class_pred),
                    msg="All predictions should contain class 3 and 8",
                )
                self.assertAlmostEqual(probs[0], correct_probabilities[i], places=3)
                self.assertTrue(len(where_not_nan) == 2, msg="Should not get any nans.")

    @unittest.skipIf(SKIP_GP_TESTS, "GP train / predict - long runtime.")
    def test_predictions_gp(self):
        self.train_grid_classifiers(["gp"], di=np.arange(0, 1000)[::6])
        tup_out = self.cls_obj_grid.get_class_predictions(
            "gp", self.TEST_INPUT_POINTS, return_ids=False
        )
        class_pred, probs, where_not_nan = tup_out

        self.assertTrue(
            (3 in class_pred) and (8 in class_pred),
            msg="All predictions should contain class 3 and 8",
        )
        self.assertAlmostEqual(probs[0], 0.5883452, places=3)
        self.assertTrue(len(where_not_nan) == 2, msg="Should not get any nans.")

    def test_pred_train_err(self):
        # Trying to predict without training
        names = ["grid", "random"]
        for i, classifier in enumerate([self.cls_obj_grid, self.cls_obj_rand]):
            with self.subTest(classifier_name=names[i]):
                with self.assertRaisesRegex(
                    Exception, "No trained interpolators exist"
                ):
                    classifier.get_class_predictions("linear", [[0, 0, 0]])

    def test_pred_linear_err(self):
        self.cls_obj_rand.train("linear", di=np.arange(0, 1000, 50))
        tup_out = self.cls_obj_rand.get_class_predictions(
            "lin", [[-1, -1, -1], [1, 1, 1]], return_ids=False
        )
        self.assertTrue(len(tup_out[2]) == 0, msg="Should return no valid values.")

    # def test_cross_val(self):
    #     correct_ans = [67.36842105263158, 66.66666666666666]
    #     acc, times = self.cls_obj_grid.cross_validate(
    #         ["rbf", "linear"], 0.05, verbose=False
    #     )
    #     for i, percent_acc in enumerate(acc):
    #         with self.subTest("Cross Val", i=i, percent_acc=percent_acc):
    #             self.assertAlmostEqual(acc[i], correct_ans[i], places=3)

    @unittest.skipIf(SKIP_GP_TESTS, "GP cross_val - long runtime.")
    def test_cross_val_gp(self):
        correct_ans = [73.76470588235294]
        acc, times = self.cls_obj_grid.cross_validate(["gp"], 0.15, verbose=False)
        for i, percent_acc in enumerate(acc):
            with self.subTest("Cross Val", i=i, percent_acc=percent_acc):
                self.assertAlmostEqual(acc[i], correct_ans[i], places=3)

    @unittest.skipIf(SKIP_TEST_PLOTS, "Skipping maximum class P plot.")
    def test_max_cls_plot(self):
        N = int(2e4) if SHOW_PLOTS else 100
        self.train_grid_classifiers(["rbf"])
        fig, axes = self.cls_obj_grid.make_max_cls_plot(
            "rbf", ("input_1", "input_2"), N=N, s=3, alpha=0.6, cmap="bone"
        )
        if SHOW_PLOTS:
            fig.show()
        else:
            print("To show plots set SHOW_PLOTS to True.")
            plt.close(fig)


if __name__ == "__main__":
    unittest.main()
