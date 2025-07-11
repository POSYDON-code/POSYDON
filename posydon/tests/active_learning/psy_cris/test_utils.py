"""Unit test for posydon.active_learning.psy_cris classes
"""
import unittest
import numpy as np
import pandas as pd
import os
from posydon.config import PATH_TO_POSYDON

from posydon.active_learning.psy_cris.utils import (
    parse_inifile,
    get_new_query_points,
    check_dist,
    get_regular_grid_df,
    get_random_grid_df,
)


class TestUtils(unittest.TestCase):
    """Test methods in utils."""

    @classmethod
    def setUpClass(cls):
        np.random.seed(12345)
        psy_cris_dir = os.path.join(
            PATH_TO_POSYDON, "posydon/active_learning/psy_cris")
        cls.INI_FILE_PATH = os.path.join(psy_cris_dir,
                                         "run_params/psycris_default.ini")

    def test_parse_inifile(self):
        self.assertTrue(os.path.isfile(self.INI_FILE_PATH), msg="Can't find file.")
        my_kwargs = parse_inifile(self.INI_FILE_PATH)
        self.assertTrue(isinstance(my_kwargs, dict))
        return my_kwargs

    def test_get_new_query_points(self):
        my_kwargs = self.test_parse_inifile()
        holder = my_kwargs["TableData_kwargs"]
        holder["my_DataFrame"] = get_regular_grid_df(N=10 ** 3, dim=3)
        my_kwargs["TableData_kwargs"] = holder

        holder_1 = my_kwargs["Sampler_kwargs"]
        holder_1["N_tot"] = 50
        holder_1["T_max"] = 5
        holder_1["verbose"] = False
        my_kwargs["Sampler_kwargs"] = holder_1

        query_pts, preds = get_new_query_points(3, **my_kwargs)
        self.assertTrue(len(query_pts) == 3)

    def test_check_dist(self):
        original_pts = np.random.uniform(
            low=(-1, -1, -1), high=(1, 1, 1), size=(500, 3)
        )
        proposed_pts = get_regular_grid_df(N=10 ** 3, dim=3).values[:, 0:3]
        result = check_dist(original_pts, proposed_pts, threshold=1e-2)
        self.assertTrue(
            sum(result) == len(proposed_pts),
            msg="All points should not be within 1e-2 of eachother.",
        )

    def test_get_regular_grid_df(self):
        for config in [dict(N=1000, dim=3), dict(N=50, dim=2), dict(jitter=True)]:
            with self.subTest("regular_grid_df", config=config):
                df = get_regular_grid_df(**config)
                self.assertTrue(isinstance(df, pd.DataFrame))

    def test_get_random_grid_df(self):
        for config in [dict(N=1000, dim=3), dict(N=50, dim=2)]:
            with self.subTest("random_grid_df", config=config):
                df = get_random_grid_df(**config)
                self.assertTrue(isinstance(df, pd.DataFrame))


if __name__ == "__main__":
    unittest.main()
