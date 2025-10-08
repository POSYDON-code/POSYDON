from unittest import TestCase

import numpy as np

from posydon.interpolation.data_scaling import DataScaler


class DataScaler_test(TestCase):
    def setUp(self):
        self.sc = DataScaler()
        self.x = np.array([1,2,3,4])
        self.y = -self.x.copy()

    def test_fit(self):
        # not a 1D array
        with self.assertRaises(AssertionError):
            self.sc.fit([12,2,3]) #list
        with self.assertRaises(AssertionError):
            self.sc.fit(np.ones((5,1)))  # list
        # default value 'none'
        with self.subTest(i=0):
            self.sc.fit(self.x)
            self.assertIsInstance(self.sc.params, list)
            self.assertEqual(self.sc.method,'none')
            self.assertEqual(len(self.sc.params),0)
        # min_max
        with self.subTest(i=1):
            self.sc.fit(self.x, method='min_max')
            self.assertEqual(self.sc.method, 'min_max')
            self.assertEqual(len(self.sc.params),2)
            self.assertEqual(self.sc.params[0], 1)
            self.assertEqual(self.sc.params[1], 4)
            self.assertEqual(self.sc.lower, -1)
            self.assertEqual(self.sc.upper, 1)
        # min_max modifying lower/upper
        with self.subTest(i=2):
            with self.assertRaises(AssertionError):
                self.sc.fit(self.x, method='min_max', lower=2)
            self.sc.fit(self.x, method='min_max', lower=-2, upper=0.5)
            self.assertEqual(self.sc.params[0], 1)
            self.assertEqual(self.sc.params[1], 4)
            self.assertEqual(self.sc.lower, -2)
            self.assertEqual(self.sc.upper, 0.5)
        # max_abs
        with self.subTest(i=3):
            self.sc.fit(self.x, method='max_abs')
            self.assertEqual(self.sc.method, 'max_abs')
            self.assertEqual(len(self.sc.params), 1)
            self.assertEqual(self.sc.params[0], 4)
            self.sc.fit(self.y, method='max_abs') # check with negative numbers
            self.assertEqual(self.sc.params[0], 4)
        # standarize
        with self.subTest(i=4):
            self.sc.fit(self.x, method='standarize')
            self.assertEqual(self.sc.method, 'standarize')
            self.assertEqual(len(self.sc.params), 2)
            self.assertEqual(self.sc.params[0], np.mean(self.x))
            self.assertEqual(self.sc.params[1], np.std(self.x))
        # log_min_max
        with self.subTest(i=5):
            self.sc.fit(self.x, method='log_min_max')
            self.assertEqual(self.sc.method, 'log_min_max')
            self.assertEqual(len(self.sc.params), 2)
            self.assertEqual(self.sc.params[0], 0)
            self.assertEqual(self.sc.params[1], np.log10(4))
            self.assertEqual(self.sc.lower, -1)
            self.assertEqual(self.sc.upper, 1)
        # log_min_max modifying lower/upper
        with self.subTest(i=6):
            with self.assertRaises(AssertionError):
                self.sc.fit(self.x, method='log_min_max', lower=2)
            self.sc.fit(self.x, method='log_min_max', lower=-2, upper=0.5)
            self.assertEqual(self.sc.params[0], 0)
            self.assertEqual(self.sc.params[1], np.log10(4))
            self.assertEqual(self.sc.lower, -2)
            self.assertEqual(self.sc.upper, 0.5)
        # log_max_abs
        with self.subTest(i=7):
            self.sc.fit(self.x, method='log_max_abs')
            self.assertEqual(self.sc.method, 'log_max_abs')
            self.assertEqual(len(self.sc.params), 1)
            self.assertEqual(self.sc.params[0], np.log10(4))
            self.sc.fit(self.y, method='log_max_abs')  # check with negative numbers
            self.assertTrue(np.isnan(self.sc.params[0]))
        # log_standarize
        with self.subTest(i=8):
            self.sc.fit(self.x, method='log_standarize')
            self.assertEqual(self.sc.method, 'log_standarize')
            self.assertEqual(len(self.sc.params), 2)
            self.assertEqual(self.sc.params[0], np.mean(np.log10(self.x)))
            self.assertEqual(self.sc.params[1], np.std(np.log10(self.x)))
        # wrong method string
        with self.assertRaises(ValueError):
            self.sc.fit(self.x, method='wrong')

    def test_transform(self):
        # check .fit has been run first
        with self.assertRaises(AssertionError):
            sc = DataScaler()
            sc.transform(self.x)
        # default value 'none'
        with self.subTest(i=0):
            self.sc.fit(self.x)
            xt = self.sc.transform(self.x)
            self.assertIsInstance(xt, np.ndarray)
            self.assertEqual(len(xt.shape),1)
            self.assertEqual(np.sum(np.abs(xt-self.x)),0)
        # min_max
        with self.subTest(i=1):
            self.sc.fit(self.x, method='min_max')
            xt = self.sc.transform(self.x)
            self.assertAlmostEqual(xt.min(),self.sc.lower)
            self.assertAlmostEqual(xt.max(), self.sc.upper)
        # min_max modifying lower/upper
        with self.subTest(i=2):
            self.sc.fit(self.x, method='min_max', lower=-2, upper=0.5)
            xt = self.sc.transform(self.x)
            self.assertAlmostEqual(xt.min(), self.sc.lower)
            self.assertAlmostEqual(xt.max(), self.sc.upper)
        # max_abs
        with self.subTest(i=3):
            self.sc.fit(self.x, method='max_abs')
            xt = self.sc.transform(self.x)
            self.assertEqual(np.abs(xt).max(), 1)
            self.assertGreaterEqual(np.abs(xt).min(),-1)
            self.sc.fit(self.y, method='max_abs') # check with negative numbers
            xt = self.sc.transform(self.x)
            self.assertEqual(np.abs(xt).max(), 1)
            self.assertGreaterEqual(np.abs(xt).min(), -1)
        # standarize
        with self.subTest(i=4):
            self.sc.fit(self.x, method='standarize')
            xt = self.sc.transform(self.x)
            self.assertAlmostEqual(xt.mean(),0)
            self.assertAlmostEqual(xt.std(), 1)
        # log_min_max
        with self.subTest(i=5):
            self.sc.fit(self.x, method='log_min_max')
            xt = self.sc.transform(self.x)
            self.assertAlmostEqual(xt.min(), self.sc.lower)
            self.assertAlmostEqual(xt.max(), self.sc.upper)
        # log_min_max modifying lower/upper
        with self.subTest(i=6):
            self.sc.fit(self.x, method='log_min_max', lower=-2, upper=0.5)
            xt = self.sc.transform(self.x)
            self.assertAlmostEqual(xt.min(), self.sc.lower)
            self.assertAlmostEqual(xt.max(), self.sc.upper)
        # log_max_abs
        with self.subTest(i=7):
            self.sc.fit(self.x, method='log_max_abs')
            xt = self.sc.transform(self.x)
            self.assertEqual(np.abs(xt).max(), 1)
            self.assertGreaterEqual(np.abs(xt).min(), -1)
        # log_standarize
        with self.subTest(i=8):
            self.sc.fit(self.x, method='log_standarize')
            xt = self.sc.transform(self.x)
            self.assertAlmostEqual(xt.mean(), 0)
            self.assertAlmostEqual(xt.std(), 1)

    def test_fit_and_transform(self):

        # default value 'none'
        with self.subTest(i=0):
            xt = self.sc.fit_and_transform(self.x)
            self.assertIsInstance(xt, np.ndarray)
            self.assertEqual(len(xt.shape),1)
            self.assertEqual(np.sum(np.abs(xt-self.x)),0)
        # min_max
        with self.subTest(i=1):
            xt = self.sc.fit_and_transform(self.x, method='min_max')
            self.assertAlmostEqual(xt.min(),self.sc.lower)
            self.assertAlmostEqual(xt.max(), self.sc.upper)
        # min_max modifying lower/upper
        with self.subTest(i=2):
            xt = self.sc.fit_and_transform(self.x, method='min_max', lower=-2, upper=0.5)
            self.assertAlmostEqual(xt.min(), self.sc.lower)
            self.assertAlmostEqual(xt.max(), self.sc.upper)
        # max_abs
        with self.subTest(i=3):
            xt = self.sc.fit_and_transform(self.x, method='max_abs')
            self.assertEqual(np.abs(xt).max(), 1)
            self.assertGreaterEqual(np.abs(xt).min(),-1)
            xt = self.sc.fit_and_transform(self.y, method='max_abs') # check with negative numbers
            self.assertEqual(np.abs(xt).max(), 1)
            self.assertGreaterEqual(np.abs(xt).min(), -1)
        # standarize
        with self.subTest(i=4):
            xt = self.sc.fit_and_transform(self.x, method='standarize')
            self.assertAlmostEqual(xt.mean(),0)
            self.assertAlmostEqual(xt.std(), 1)
        # log_min_max
        with self.subTest(i=5):
            xt = self.sc.fit_and_transform(self.x, method='log_min_max')
            self.assertAlmostEqual(xt.min(), self.sc.lower)
            self.assertAlmostEqual(xt.max(), self.sc.upper)
        # log_min_max modifying lower/upper
        with self.subTest(i=6):
            xt = self.sc.fit_and_transform(self.x, method='log_min_max', lower=-2, upper=0.5)
            self.assertAlmostEqual(xt.min(), self.sc.lower)
            self.assertAlmostEqual(xt.max(), self.sc.upper)
        # log_max_abs
        with self.subTest(i=7):
            xt = self.sc.fit_and_transform(self.x, method='log_max_abs')
            self.assertEqual(np.abs(xt).max(), 1)
            self.assertGreaterEqual(np.abs(xt).min(), -1)
        # log_standarize
        with self.subTest(i=8):
            xt = self.sc.fit_and_transform(self.x, method='log_standarize')
            self.assertAlmostEqual(xt.mean(), 0)
            self.assertAlmostEqual(xt.std(), 1)

    def test_inv_transform(self):
        # default value 'none'
        with self.subTest(i=0):
            xt = self.sc.fit_and_transform(self.x)
            self.assertAlmostEqual(np.sum(np.abs(self.sc.inv_transform(xt) - self.x)), 0)
        # min_max
        with self.subTest(i=1):
            xt = self.sc.fit_and_transform(self.x, method='min_max')
            self.assertAlmostEqual(np.sum(np.abs(self.sc.inv_transform(xt) - self.x)), 0)
        # min_max modifying lower/upper
        with self.subTest(i=2):
            xt = self.sc.fit_and_transform(self.x, method='min_max', lower=-2, upper=0.5)
            self.assertAlmostEqual(np.sum(np.abs(self.sc.inv_transform(xt) - self.x)), 0)
        # max_abs
        with self.subTest(i=3):
            xt = self.sc.fit_and_transform(self.x, method='max_abs')
            self.assertAlmostEqual(np.sum(np.abs(self.sc.inv_transform(xt) - self.x)), 0)
            xt = self.sc.fit_and_transform(self.y, method='max_abs')  # check with negative numbers
            self.assertAlmostEqual(np.sum(np.abs(self.sc.inv_transform(xt) - self.y)), 0)
        # standarize
        with self.subTest(i=4):
            xt = self.sc.fit_and_transform(self.x, method='standarize')
            self.assertAlmostEqual(np.sum(np.abs(self.sc.inv_transform(xt) - self.x)), 0)
        # log_min_max
        with self.subTest(i=5):
            xt = self.sc.fit_and_transform(self.x, method='log_min_max')
            self.assertAlmostEqual(np.sum(np.abs(self.sc.inv_transform(xt) - self.x)), 0)
        # log_min_max modifying lower/upper
        with self.subTest(i=6):
            xt = self.sc.fit_and_transform(self.x, method='log_min_max', lower=-2, upper=0.5)
            self.assertAlmostEqual(np.sum(np.abs(self.sc.inv_transform(xt) - self.x)), 0)
        # log_max_abs
        with self.subTest(i=7):
            xt = self.sc.fit_and_transform(self.x, method='log_max_abs')
            self.assertAlmostEqual(np.sum(np.abs(self.sc.inv_transform(xt) - self.x)), 0)
        # log_standarize
        with self.subTest(i=8):
            xt = self.sc.fit_and_transform(self.x, method='log_standarize')
            self.assertAlmostEqual(np.sum(np.abs(self.sc.inv_transform(xt) - self.x)), 0)
