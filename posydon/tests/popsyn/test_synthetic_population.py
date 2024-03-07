import unittest
import os
import numpy as np
import pandas as pd


from posydon.popsyn.synthetic_population import (
    PopulationRunner,
    History,
    Oneline,
    PopulationIO,
    parameter_array,
    Population
)


# Test the PopulationRunner class
class TestPopulationRunner(unittest.TestCase):
    # Test the initialisation of the PopulationRunner class
    def test_init(self):
        # Test the initialisation of the PopulationRunner class
        poprun = PopulationRunner('../../popsyn/population_params_default.ini', verbose=True)

        # Check if the verbose attribute is set correctly
        self.assertEqual(poprun.verbose, True, 'Verbose attribute is not set correctly')

        # Check if the solar_metallicities attribute is a list
        self.assertIsInstance(poprun.solar_metallicities, list, 'solar_metallicities attribute is not a list')

        # Check if the binary_populations attribute is a list
        self.assertIsInstance(poprun.binary_populations, list, 'binary_populations attribute is not a list')

    # Test the evolve method
    def test_changed_binarypop(self):
        poprun = PopulationRunner('../../popsyn/population_params_default.ini')
        # test 
        self.assertEqual(poprun.binary_populations[0].metallicity, 0.0001)
        # Check if the temp_directory attribute is set correctly
        self.assertEqual(poprun.binary_populations[0].kwargs['temp_directory'],
                         '1e-04_Zsun_batches',
                         'temp_directory attribute is not set correctly')
        
        # Check if the binary_populations attribute is not empty after calling the evolve method
        self.assertTrue(poprun.binary_populations, 'binary_populations attribute is empty after calling the evolve method')


# Test the History class
class TestHistory(unittest.TestCase):
    
    @classmethod
    def setUpClass(self):
        # Set up a test HDF5 file using pandas HDFStore
        self.filename = 'test_population.h5'
        with pd.HDFStore(self.filename, 'w') as store:
            # Create a history dataframe
            history_data = pd.DataFrame({'time': [1, 2, 3], 'event': ['ZAMS','oRLO1', 'CEE']})
            store.append('history',history_data, data_columns=True)

    def setUp(self):
        self.history = History(self.filename, verbose=False, chunksize=10000)

    def test_init(self):
        history = History(self.filename, verbose=False, chunksize=10000)
        self.assertEqual(history.filename, self.filename, 'Filename is not set correctly')
        self.assertFalse(history.verbose, 'Verbose attribute is not set correctly')
        self.assertEqual(history.chunksize, 10000, 'Chunksize attribute is not set correctly')
        
        expected_lengths = pd.DataFrame(index=[0, 1, 2],data={'index': [1, 1, 1]})
        expected_lengths.index.name = 'index'
        pd.testing.assert_frame_equal(history.lengths, expected_lengths, 'Lengths attribute is not equal to the expected dataframe')

        #self.assertEquals(history.lengths, expected_lengths, 'Lengths attribute is not equal to the expected dataframe')
        
        self.assertEqual(history.number_of_systems, 3, 'Number of systems attribute is not None')
        self.assertEqual(history.columns.to_list(), ['time', 'event'], 'Columns attribute is not None')
        
        self.assertIsInstance(history.indices, np.ndarray, 'Indices attribute is not an ndarray')
        np.testing.assert_array_equal(history.indices, np.array([0, 1, 2]), 'Indices attribute is not equal to the expected list')
        
        with self.assertRaises(KeyError):
            History('invalid_filename.h5', verbose=False, chunksize=10000)


    def test_getitem_single_index(self):
        df = self.history[0]
        self.assertIsInstance(df, pd.DataFrame, 'Returned object is not a DataFrame')
        self.assertEqual(len(df), 1, 'Returned DataFrame does not have the correct length')
        

    def test_getitem_multiple_indices(self):
        df = self.history[[0, 1, 2]]
        self.assertIsInstance(df, pd.DataFrame, 'Returned object is not a DataFrame')
        self.assertEqual(len(df), 3, 'Returned DataFrame does not have the correct length')

    def test_getitem_index_array(self):
        indices = np.array([0, 1, 2])
        df = self.history[indices]
        self.assertIsInstance(df, pd.DataFrame, 'Returned object is not a DataFrame')
        self.assertEqual(len(df), 3, 'Returned DataFrame does not have the correct length')

    def test_getitem_single_column(self):
        column = 'time'
        df = self.history[column]
        self.assertIsInstance(df, pd.DataFrame, 'Returned object is not a DataFrame')
        self.assertEqual(len(df.columns), 1, 'Returned DataFrame does not have the correct number of columns')

    def test_getitem_boolean_mask_numpy(self):
        mask = (self.history['time'] > 1).to_numpy()
        df = self.history[mask]
        self.assertIsInstance(df, pd.DataFrame, 'Returned object is not a DataFrame')

    def test_getitem_boolean_mask_pandas(self):
        mask = self.history['time'] > 1
        df = self.history[mask]
        self.assertIsInstance(df, pd.DataFrame, 'Returned object is not a DataFrame')

    def test_getitem_multiple_columns(self):
        columns = ['time', 'event']
        df = self.history[columns]
        self.assertIsInstance(df, pd.DataFrame, 'Returned object is not a DataFrame')
        self.assertEqual(len(df.columns), 2, 'Returned DataFrame does not have the correct number of columns')

    def test_getitem_invalid_key(self):
        with self.assertRaises(ValueError):
            self.history['invalid_key']

    def test_len(self):
        length = len(self.history)
        self.assertIsInstance(length, int, 'Returned object is not an integer')
        self.assertEqual(length, 3, 'Returned length is not correct')
        
    def test_head(self):
        n = 2
        df = self.history.head(n)
        self.assertIsInstance(df, pd.DataFrame, 'Returned object is not a DataFrame')
        self.assertEqual(len(df), n, 'Returned DataFrame does not have the correct length')

    def test_tail(self):
        n = 2
        df = self.history.tail(n)
        self.assertIsInstance(df, pd.DataFrame, 'Returned object is not a DataFrame')
        self.assertEqual(len(df), n, 'Returned DataFrame does not have the correct length')

    def test_repr(self):
        representation = self.history.__repr__()
        self.assertIsInstance(representation, str, 'Returned object is not a string')

    def test_repr_html(self):
        html_representation = self.history._repr_html_()
        self.assertIsInstance(html_representation, str, 'Returned object is not a string')


    def test_select(self):
        df = self.history.select(where="time > 1", start=0, stop=10, columns=['event', 'time'])
        self.assertIsInstance(df, pd.DataFrame, 'Returned object is not a DataFrame')
        self.assertEqual(len(df.columns), 2, 'Returned DataFrame does not have the correct number of columns')
        self.assertEqual(len(df), 2, 'Returned DataFrame does not have the correct length')
        
        
     
    @classmethod   
    def tearDownClass(self):
        os.remove(self.filename)

# Test the Oneline class
class TestOneline(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        # Set up a test HDF5 file using pandas HDFStore
        cls.filename = 'test_oneline.h5'
        with pd.HDFStore(cls.filename, 'w') as store:
            # Create a oneline dataframe
            oneline_data = pd.DataFrame({'time': [1, 2, 3], 'S1_mass_i': ['30','30', '70']})
            store.append('oneline', oneline_data, data_columns=True)

    def setUp(self):
        self.oneline = Oneline(self.filename, verbose=False, chunksize=10000)

    def test_init(self):
        oneline = Oneline(self.filename, verbose=True, chunksize=5000)

        self.assertEqual(oneline.filename, self.filename, 'Filename is not set correctly')
        self.assertTrue(oneline.verbose, 'Verbose attribute is not set correctly')
        self.assertEqual(oneline.chunksize, 5000, 'Chunksize attribute is not set correctly')
        self.assertEqual(oneline.number_of_systems, 3, 'Number of systems attribute is not set correctly')
        self.assertEqual(oneline.columns.to_list(), ['time', 'S1_mass_i'], 'Columns attribute is not set correctly')
        self.assertEqual(oneline.number_of_systems, 3, 'Number of systems attribute is not set correctly')
        
        self.assertIsInstance(oneline.indices, np.ndarray, 'Indices attribute is not an ndarray')
        np.testing.assert_array_equal(oneline.indices, np.array([0, 1, 2]), 'Indices attribute is not equal to the expected list')

        with self.assertRaises(KeyError):
            Oneline('invalid_filename.h5', verbose=False, chunksize=10000)

    def test_getitem_single_index(self):
        df = self.oneline[0]
        self.assertIsInstance(df, pd.DataFrame, 'Returned object is not a DataFrame')
        self.assertEqual(len(df), 1, 'Returned DataFrame does not have the correct length')
        
    def test_getitem_multiple_indices(self):
        df = self.oneline[[0, 1, 2]]
        self.assertIsInstance(df, pd.DataFrame, 'Returned object is not a DataFrame')
        self.assertEqual(len(df), 3, 'Returned DataFrame does not have the correct length')
        
    def test_getitem_index_array(self):
        indices = np.array([0, 1, 2])
        df = self.oneline[indices]
        self.assertIsInstance(df, pd.DataFrame, 'Returned object is not a DataFrame')
        self.assertEqual(len(df), 3, 'Returned DataFrame does not have the correct length')
        pd.testing.assert_frame_equal(df, 
                                      pd.DataFrame({'time': [1, 2, 3],
                                                    'S1_mass_i': ['30','30', '70']}), 
                                      'Returned DataFrame is not equal to the expected DataFrame')

        
    def test_getitem_single_column(self):
        column = 'time'
        df = self.oneline[column]
        self.assertIsInstance(df, pd.DataFrame, 'Returned object is not a DataFrame')
        self.assertEqual(len(df.columns), 1, 'Returned DataFrame does not have the correct number of columns')
    
    def test_getitem_boolean_mask_numpy(self):
        mask = (self.oneline['time'] > 1).to_numpy().flatten()
        df = self.oneline[mask]
        self.assertIsInstance(df, pd.DataFrame, 'Returned object is not a DataFrame')
        
    def test_getitem_boolean_mask_pandas(self):
        mask = self.oneline['time'] > 1
        print(mask)
        df = self.oneline[mask]
        self.assertIsInstance(df, pd.DataFrame, 'Returned object is not a DataFrame')
        
    def test_getitem_multiple_columns(self):
        columns = ['time', 'S1_mass_i']
        df = self.oneline[columns]
        self.assertIsInstance(df, pd.DataFrame, 'Returned object is not a DataFrame')
        self.assertEqual(len(df.columns), 2, 'Returned DataFrame does not have the correct number of columns')
        
    def test_getitem_invalid_key(self):
        with self.assertRaises(ValueError):
            self.oneline['invalid_key']
            
    def test_len(self):
        length = len(self.oneline)
        self.assertIsInstance(length, int, 'Returned object is not an integer')
        self.assertEqual(length, 3, 'Returned length is not correct')
        
    def test_head(self):
        n = 2
        df = self.oneline.head(n)
        self.assertIsInstance(df, pd.DataFrame, 'Returned object is not a DataFrame')
        self.assertEqual(len(df), n, 'Returned DataFrame does not have the correct length')
        
    def test_tail(self):
        n = 2
        df = self.oneline.tail(n)
        self.assertIsInstance(df, pd.DataFrame, 'Returned object is not a DataFrame')
        self.assertEqual(len(df), n, 'Returned DataFrame does not have the correct length')
        
    def test_repr(self):
        representation = self.oneline.__repr__()
        self.assertIsInstance(representation, str, 'Returned object is not a string')
    
    def test_repr_html(self):
        html_representation = self.oneline._repr_html_()
        self.assertIsInstance(html_representation, str, 'Returned object is not a string')
    
    def test_select(self):
        df = self.oneline.select(where="time > 1", start=0, stop=10, columns=['S1_mass_i', 'time'])
        self.assertIsInstance(df, pd.DataFrame, 'Returned object is not a DataFrame')
        self.assertEqual(len(df.columns), 2, 'Returned DataFrame does not have the correct number of columns')
        self.assertEqual(len(df), 2, 'Returned DataFrame does not have the correct length')
    
        
    @classmethod
    def tearDownClass(self):
        os.remove(self.filename)


# Test the PopulationIO class
class TestPopulationIO(unittest.TestCase):

    def setUp(self):
        self.filename = "test_population.hdf5"

    def tearDown(self):
        if os.path.exists(self.filename):
            os.remove(self.filename)
            
    def test_init(self):
        pop_io = PopulationIO()
        self.assertEqual(pop_io.verbose, False, "Verbose attribute is not set correctly")

    def test_save_and_load_mass_per_met(self):
        population_io = PopulationIO()
        population_io.mass_per_met = pd.DataFrame({"metallicity": [0.02, 0.04], "mass": [1.0, 2.0]})
        population_io._save_mass_per_met(self.filename)
        
        loaded_io = PopulationIO()
        loaded_io._load_mass_per_met(self.filename)
        pd.testing.assert_frame_equal(population_io.mass_per_met, loaded_io.mass_per_met)

    def test_save_and_load_ini_params(self):
        population_io = PopulationIO()
        population_io.ini_params = {i:10 for i in parameter_array}
        population_io._save_ini_params(self.filename)
        
        loaded_io = PopulationIO()
        loaded_io._load_ini_params(self.filename)

        self.assertEqual(population_io.ini_params, 
                         loaded_io.ini_params, "Loaded ini_params are not equal to the saved ini_params")



if __name__ == '__main__':
    unittest.main()
    