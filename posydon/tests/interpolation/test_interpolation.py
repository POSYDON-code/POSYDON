# from unittest import TestCase
#
# import numpy as np
# import posydon.grids.psygrid as psg
# import posydon.interpolation.interpolation as psi
#
# try:
#     import gpflow
# except ImportError:
#     print("Import Error for TensorFlow and/or GPFlow, most, if not all "
#           "features of the psyInterp class will not work, please check your installation "
#           "of gpflow or tensorflow or install the correct gpflow by running pip install .[ml]")
#
#
# class Interpolation_test(TestCase):
#     def setUp(self):
#         # FIX PATH YOU CANNOT GIVE A LOCAL PATH
#         self.grid = psg.PSyGrid()
#         self.grid.load("/home/juanga/Desktop/data/grid_BH_He_star.h5")
#         self.input_keys = self.grid.initial_values.dtype.names
#         self.output_keys = self.grid.final_values.dtype.names[2:4]
#         self.input_norms = ['log_min_max', 'log_min_max', 'log_min_max']
#         self.output_norms = ['log_standarize', 'log_standarize']
#
#     def test_init(self):
#         m = psi.psyInterp(grid=self.grid,
#                           in_keys=self.input_keys,
#                           out_keys=self.output_keys,
#                           in_scaling=self.input_norms,
#                           out_scaling=self.output_norms)
#         self.assertEqual(len(m.in_keys), len(self.input_keys))
#         self.assertEqual(len(m.out_keys), len(self.output_keys))
#         self.assertEqual(m.XYT.shape[0], m.N)
#         self.assertEqual(m.XYT.shape[1], m.n_in+m.n_out)
#
# class SGPInterp_test(TestCase):
#     pass
