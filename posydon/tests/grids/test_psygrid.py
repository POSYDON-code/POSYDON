"""Unit test module for PSyGrid class."""

# import unittest
# import os
# import shutil
# import zipfile
# from posydon.config import PATH_TO_POSYDON
# from posydon.grids.psygrid import PSyGrid
#
# 
# class TestPSyGrid(unittest.TestCase):
#     """Class for unit-testing the PSyGrid object."""
#
#     @classmethod
#     def setUpClass(cls):
#         """Prepare the data for the individual tests."""
#         # find the input/output paths
#         path_to_workdir = os.path.join(
#             PATH_TO_POSYDON, "posydon/tests/data/POSYDON-UNIT-TESTS/grids/")
#         path_to_zipfile = os.path.join(
#             path_to_workdir, "sample_HMSHMS_grid.zip")
#         cls.path_to_extract = os.path.join(path_to_workdir, "tmp")
#         cls.path_to_psygrid = os.path.join(path_to_workdir, "tmp.h5")
#         cls.path_to_psyslim = os.path.join(path_to_workdir, "tmp_slim.h5")
#
#         # unzip the test data
#         with zipfile.ZipFile(path_to_zipfile) as zipf:
#             zipf.extractall(path=cls.path_to_extract)
#
#         # try to create the PSyGrid objects, and store whether they failed
#         try:
#             psygrid = PSyGrid()
#             psygrid.create(cls.path_to_extract, cls.path_to_psygrid)
#             del psygrid
#             cls.failure_msg = None
#         except Exception as e:
#             cls.failure_msg = str(e)
#         try:
#             psygrid = PSyGrid()
#             psygrid.create(cls.path_to_extract, cls.path_to_psyslim, slim=True)
#             del psygrid
#             cls.failure_msg_slim = None
#         except Exception as e:
#             cls.failure_msg_slim = str(e)
#
#     @classmethod
#     def tearDownClass(cls):
#         """Remove the unzipped an newly created files."""
#         shutil.rmtree(cls.path_to_extract)
#         os.remove(cls.path_to_psygrid)
#         os.remove(cls.path_to_psyslim)
#
#     def setUp(self):
#         """Open the grids before each test."""
#         # if failed to the create the grid objects, do not continue
#         if self.failure_msg is not None:
#             self.fail("Cannot create a PSyGrid object: "
#                       + self.failure_msg)
#
#         if self.failure_msg_slim is not None:
#             self.fail("Cannot create a slim PSyGrid object: "
#                       + self.failure_msg_slim)
#
#         # load the grids before any testing
#         try:
#             self.psygrid = PSyGrid(self.path_to_psygrid)
#         except Exception as e:
#             self.fail("Cannot load the PSyGrid object: " + str(e))
#
#         try:
#             self.psygrid_slim = PSyGrid(self.path_to_psyslim)
#         except Exception as e:
#             self.fail("Cannot load the slim PSyGrid object: " + str(e))
#
#         # keep them together for use in loops
#         self.grids = [self.psygrid, self.psygrid_slim]
#
#     def tearDown(self):
#         """Close and delete the grid objects after each test."""
#         del self.psygrid
#         del self.psygrid_slim
#
#     def test_filename(self):
#         """Check grid paths, and that PSyGrid objects know them."""
#         self.assertEqual(self.psygrid.filepath, self.path_to_psygrid)
#         self.assertEqual(self.psygrid_slim.filepath, self.path_to_psyslim)
#
#     def test_print_and_str(self):
#         """Test whether print works with grids (or `str` method)."""
#         for grid in self.grids:
#             s = str(grid)
#             self.assertTrue(len(s) > 0)
#
#     def test_config(self):
#         """Check that the ConfigFile instance works and reports HDF5 format."""
#         for grid in self.grids:
#             self.assertEqual(grid.config.format, "hdf5")
#
#     def test_number_of_runs(self):
#         """The normal grid must have run data, while the slim one must not."""
#         self.assertTrue(self.psygrid.n_runs > 0)
#         self.assertTrue(self.psygrid_slim.n_runs == 0)
#
#     def test_in_keyword(self):
#         """Check `in` keyword and extreme values of run indices."""
#         N = self.psygrid.n_runs
#         self.assertTrue(0 in self.psygrid)
#         self.assertTrue(N not in self.psygrid)
#         self.assertTrue(0 not in self.psygrid_slim)
#
#     def test_len_list_set(self):
#         """Test the `len`, `list` and `set` methods."""
#         N = self.psygrid.n_runs
#         self.assertEqual(N, len(self.psygrid))
#         self.assertEqual(N, len(list(self.psygrid)))
#         self.assertEqual(N, len(set(self.psygrid)))
#
#     def test_data_integrity(self):
#         """Test whether the run data agree with the initial/finial values."""
#         self.assertEqual(self.psygrid[0].binary_history["age"][0],
#                          self.psygrid.initial_values[0]["age"])
#         self.assertEqual(self.psygrid.initial_values["age"][0],
#                          self.psygrid.initial_values[0]["age"])
#         self.assertEqual(self.psygrid.final_values["age"][0],
#                          self.psygrid[0].final_values["age"])
#         self.assertEqual(self.psygrid_slim.initial_values["age"][0],
#                          self.psygrid_slim.initial_values[0]["age"])
#
#
# if __name__ == "__main__":
#     unittest.main()
