"""Unit test for posydon.active_learning.psy_cris classes
"""
import unittest
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from posydon.active_learning.psy_cris.data import TableData

from posydon.active_learning.psy_cris.utils import (
    get_regular_grid_df,
    get_random_grid_df,
)
from posydon.active_learning.psy_cris.synthetic_data.synth_data_3D import get_output_3D

SKIP_TEST_PLOTS = False
SHOW_PLOTS = False


class TestTableData(unittest.TestCase):
    """Test TableData class on the 3d synthetic data set."""

    @classmethod
    def setUpClass(cls):
        np.random.seed(12345)
        cls.TEST_DATA_GRID = get_regular_grid_df(N=10 ** 3, dim=3)
        cls.TEST_DATA_RAND = get_random_grid_df(N=10 ** 3, dim=3)
        cls.UNIQUE_CLASSES = [1, 2, 3, 4, 6, 8]

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

    def test_init_0(self):
        td = self.create_TableData(self.TEST_DATA_GRID)
        self.assertTrue(isinstance(td, TableData))

    def test_init_1_grid(self):
        my_kwargs = {"n_neighbors": [2, 3]}
        table_grid = self.create_TableData(self.TEST_DATA_GRID, **my_kwargs)
        # N classes
        self.assertTrue(
            table_grid.num_classes == len(self.UNIQUE_CLASSES),
            msg="Should find 6 classes. Found {}".format(table_grid.num_classes),
        )
        # Unique classes + APC cols
        regr_data = table_grid.get_regr_data(what_data="output")
        for cls_key in regr_data.keys():
            with self.subTest("Checking data by class.", cls_key=cls_key):
                self.assertIn(cls_key, self.UNIQUE_CLASSES)
                self.assertIn("APC2_output_1", regr_data[cls_key].columns)
                if cls_key != 2:
                    self.assertIn("APC3_output_1", regr_data[cls_key].columns)

    def test_init_2_rand(self):
        my_kwargs = {"n_neighbors": [2, 3]}
        table_rand = self.create_TableData(self.TEST_DATA_RAND, **my_kwargs)

        # N classes
        self.assertTrue(
            table_rand.num_classes == len(self.UNIQUE_CLASSES),
            msg="Should find 6 classes. Found {}".format(table_rand.num_classes),
        )
        # Unique classes + APC cols
        regr_data = table_rand.get_regr_data(what_data="output")
        for cls_key in regr_data.keys():
            with self.subTest("Checking data by class.", cls_key=cls_key):
                self.assertIn(cls_key, self.UNIQUE_CLASSES)
                self.assertIn("APC2_output_1", regr_data[cls_key].columns)
                self.assertIn("APC3_output_1", regr_data[cls_key].columns)

    def test_init_3_clean_data(self):
        my_kwargs = {"n_neighbors": [2, 3], "omit_vals": [-1]}
        table_grid = self.create_TableData(self.TEST_DATA_GRID, **my_kwargs)

        self.assertTrue(
            len(table_grid.get_data()) == 729,
            msg="Should remove 271 rows from grid data set with value -1.",
        )

    def test_init_4_clean_data(self):
        my_kwargs = {"n_neighbors": [2, 3], "omit_vals": [-1]}
        table_rand = self.create_TableData(self.TEST_DATA_RAND, **my_kwargs)

        self.assertTrue(
            len(table_rand.get_data()) == 1000,
            msg="Should remove 0 rows from random data set with value -1.",
        )

    def test_classification_data(self):
        my_kwargs = {"n_neighbors": [2, 3], "omit_vals": [-1]}
        table_grid = self.create_TableData(self.TEST_DATA_GRID, **my_kwargs)

        binary_classification_data = table_grid.get_binary_mapping_per_class()
        self.assertTrue(
            binary_classification_data.shape == (len(self.UNIQUE_CLASSES), 729)
        )
        self.assertTrue(
            all(
                [
                    sum((row == 1) + (row == 0)) == 729
                    for row in binary_classification_data
                ]
            )
        )

    def test_regression_data_1(self):
        my_kwargs = {
            "n_neighbors": [2, 3],
        }
        table_grid = self.create_TableData(self.TEST_DATA_GRID, **my_kwargs)
        output_dat = table_grid.get_regr_data(what_data="output")
        for cls in [5, 7]:
            with self.subTest(cls=cls):
                # Raise KeyError for classes that shouldn't exist
                with self.assertRaisesRegex(KeyError, str(cls)):
                    output_dat[cls]

    def test_regression_data_2(self):
        my_kwargs = {"n_neighbors": [2, 3, 5, 10, 50]}
        table_grid = self.create_TableData(self.TEST_DATA_GRID, **my_kwargs)
        table_rand = self.create_TableData(self.TEST_DATA_RAND, **my_kwargs)

        for key, val in table_grid._regr_dfs_per_class_.items():
            col_names = list(val.columns)
            with self.subTest(key=key, data="grid"):
                self.assertTrue(any(["APC2" in item for item in col_names]))
                if key in [1, 3, 4, 8]:
                    self.assertTrue(any(["APC50" in item for item in col_names]))

        for key, val in table_rand._regr_dfs_per_class_.items():
            col_names = list(val.columns)
            with self.subTest(key=key, data="random"):
                self.assertTrue(any(["APC2" in item for item in col_names]))
                self.assertTrue(any(["APC3" in item for item in col_names]))
                if key in [1, 3, 6, 8]:
                    self.assertTrue(any(["APC50" in item for item in col_names]))

    @unittest.skipIf( SKIP_TEST_PLOTS, "Test by plotting nearest neighbors skipped")
    def test_nearest_neighbhors(self):
        my_kwargs = {
            "n_neighbors": [2, 3],
        }
        table_grid = self.create_TableData(self.TEST_DATA_GRID, **my_kwargs)
        n_neigh = 3

        dat = np.random.uniform(low=(-1, -1), high=(1, 1), size=(15, 2))
        output = table_grid.find_n_neighbors(dat, [n_neigh])

        plt.figure(figsize=(4, 4), dpi=100)
        plt.title("NearestNeighbors test")
        plt.plot(
            dat.T[0][0], dat.T[1][0], "+", color="r", markersize=10, label="reference"
        )
        plt.scatter(dat.T[0], dat.T[1], label="data")
        for i in range(n_neigh):
            plt.scatter(
                dat.T[0][output[n_neigh][0][i]],
                dat.T[1][output[n_neigh][0][i]],
                marker="+",
                color="lime",
                s=29,
                label="nearest",
            )
        plt.axis("equal")
        if SHOW_PLOTS:
            plt.show()
        else:
            print("To show plots set SHOW_PLOTS to True.")
        plt.close()

    @unittest.skipIf( SKIP_TEST_PLOTS, "Test for general plotting skipped")
    def test_plotting_1(self):
        my_kwargs = {
            "n_neighbors": [2, 3],
        }
        table_grid = self.create_TableData(self.TEST_DATA_GRID, **my_kwargs)
        table_rand = self.create_TableData(self.TEST_DATA_RAND, **my_kwargs)

        zed_val = 1
        N = 40
        X, Y = np.meshgrid(np.linspace(-1, 1, N), np.linspace(-1, 1, N))
        Z = np.ones(X.shape) * zed_val
        f_out = get_output_3D(X, Y, Z)
        print("ZED VAL: {}".format(zed_val))

        fig, subs = plt.subplots(1, 3, figsize=(14, 4), dpi=100)
        subs[0].set_title("TableData - even grid")
        fig, subs[0], handles = table_grid.make_class_data_plot(
            fig,
            subs[0],
            ["input_1", "input_2"],
            my_slice_vals={0: (0.9, 1.1)},
            return_legend_handles=True,
        )
        subs[0].legend(
            handles, table_grid._unique_class_keys_, bbox_to_anchor=(-0.25, 0.5)
        )

        subs[1].set_title("TableData - random points")
        fig, subs[1], handles = table_rand.make_class_data_plot(
            fig,
            subs[1],
            ["input_1", "input_2"],
            my_slice_vals={0: (0.9, 1.1)},
            return_legend_handles=True,
        )

        subs[2].set_title("Analytic Classification")
        subs[2].pcolormesh(X, Y, f_out["class"].values.reshape(N, N), shading="auto")
        for i in range(3):
            subs[i].axis("equal")
        if SHOW_PLOTS:
            plt.show()
        else:
            print("To show plots set SHOW_PLOTS to True.")
        plt.close()

    @unittest.skipIf( SKIP_TEST_PLOTS, "Test for plotting 3d skipped")
    def test_plotting_2(self):
        my_kwargs = {
            "n_neighbors": [2, 3],
        }
        table_grid = self.create_TableData(self.TEST_DATA_GRID, **my_kwargs)
        table_grid.plot_3D_class_data()
        plt.title("plot_3D_class_data")
        if SHOW_PLOTS:
            plt.show()
        else:
            print("To show plots set SHOW_PLOTS to True.")
        plt.close()


if __name__ == "__main__":
    unittest.main()
