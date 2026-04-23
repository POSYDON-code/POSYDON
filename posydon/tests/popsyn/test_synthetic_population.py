import os
import tempfile

import numpy as np
import pandas as pd
import pytest

from posydon.config import PATH_TO_POSYDON
from posydon.popsyn.synthetic_population import (
    History,
    Oneline,
    Population,
    PopulationIO,
    PopulationRunner,
    parameter_array,
)
from posydon.utils.constants import Zsun


# Test the PopulationRunner class
class TestPopulationRunner:
    # Test the initialisation of the PopulationRunner class
    def test_init(self):
        # Test the initialisation of the PopulationRunner class
        poprun = PopulationRunner(PATH_TO_POSYDON+'posydon/popsyn/population_params_default.ini', verbose=True)

        # Check if the verbose attribute is set correctly
        assert poprun.verbose == True, 'Verbose attribute is not set correctly'

        # Check if the solar_metallicities attribute is a list
        assert isinstance(poprun.solar_metallicities, list), 'solar_metallicities attribute is not a list'

        # Check if the binary_populations attribute is a list
        assert isinstance(poprun.binary_populations, list), 'binary_populations attribute is not a list'

    def test_init_invalid_ini_file(self):
        with pytest.raises(ValueError):
            PopulationRunner('invalid_file')

    def test_single_metallicity(self):
        # copy the default ini file to a new file
        new_ini_file = 'test_population_params.ini'
        with open(PATH_TO_POSYDON+'posydon/popsyn/population_params_default.ini', 'r') as file:
            data = file.read()
        start = data.find('metallicity')
        replace_str = 'metallicity = 0.0001'
        data_new = data[:start] + replace_str + data[start+22:]
        with open(new_ini_file, 'w') as file:
            file.write(data_new)

        poprun = PopulationRunner(new_ini_file)
        assert poprun.binary_populations[0].metallicity == 0.0001

    def test_evolve(self, mocker):
        mocker.patch('posydon.popsyn.binarypopulation.BinaryPopulation.evolve', return_value=None)
        mocker.patch('posydon.popsyn.binarypopulation.BinaryPopulation.combine_saved_files', return_value=None)

        poprun = PopulationRunner(PATH_TO_POSYDON+'/posydon/popsyn/population_params_default.ini')
        # set population to 1 binary
        for pop in poprun.binary_populations:
            pop.number_of_systems = 1

        # create a temporary directory with 1e-04_Zsun_batches
        os.makedirs('1e-04_Zsun_batches', exist_ok=True)

        poprun.evolve()
        assert poprun.binary_populations, 'binary_populations attribute is empty after calling the evolve method'


    def test_evolve_file_exists(self, mocker):
        mocker.patch('posydon.popsyn.binarypopulation.BinaryPopulation.evolve', return_value=None)
        mocker.patch('posydon.popsyn.binarypopulation.BinaryPopulation.combine_saved_files', return_value=None)

        # Create a temporary file with the 1e-04_ZSun_population.h5 name
        open('1e-04_Zsun_population.h5', 'w').close()

        poprun = PopulationRunner(PATH_TO_POSYDON+'/posydon/popsyn/population_params_default.ini')
        # set population to 1 binary
        for pop in poprun.binary_populations:
            pop.number_of_systems = 1

        with pytest.raises(FileExistsError):
            poprun.evolve()


    # Test the evolve method
    def test_changed_binarypop(self):
        poprun = PopulationRunner(PATH_TO_POSYDON+'/posydon/popsyn/population_params_default.ini')
        # test
        assert poprun.binary_populations[0].metallicity == 0.0001
        # Check if the temp_directory attribute is set correctly
        assert poprun.binary_populations[0].kwargs['temp_directory'] == '1e-04_Zsun_batches', 'temp_directory attribute is not set correctly'

        # Check if the binary_populations attribute is not empty after calling the evolve method
        assert poprun.binary_populations, 'binary_populations attribute is empty after calling the evolve method'

    @classmethod
    def teardown_class(cls):
        if os.path.exists('1e-04_Zsun_batches'):
            os.rmdir('1e-04_Zsun_batches')
        if os.path.exists('test_population_params.ini'):
            os.remove('test_population_params.ini')
        if os.path.exists('1e-04_Zsun_population.h5'):
            os.remove('1e-04_Zsun_population.h5')



# Test the History class
class TestHistory:

    @classmethod
    def setup_class(cls):
        # Set up a test HDF5 file using pandas HDFStore
        cls.filename = 'test_population.h5'
        with pd.HDFStore(cls.filename, 'w') as store:
            # Create a history dataframe
            history_data = pd.DataFrame({'time': [1, 2, 3], 'event': ['ZAMS','oRLO1', 'CEE']})
            store.append('history',history_data, data_columns=True)

        cls.filename2 = 'test_population2.h5'
        with pd.HDFStore(cls.filename2, 'w') as store:
            # Create a history dataframe
            history_data = pd.DataFrame({'time': [1, 2, 3], 'event': ['ZAMS','oRLO1', 'CEE']})
            store.append('history',history_data, data_columns=True)

    @classmethod
    def teardown_class(cls):
        os.remove(cls.filename)

    def setup_method(self):
        self.history = History(self.filename, verbose=False, chunksize=10000)

    def test_init(self):
        history = History(self.filename2, verbose=True, chunksize=10000)
        assert history.filename == self.filename2, 'Filename is not set correctly'
        assert history.verbose == True, 'Verbose attribute is not set correctly'
        assert history.chunksize == 10000, 'Chunksize attribute is not set correctly'

        expected_lengths = pd.DataFrame(index=[0, 1, 2],data={'index': [1, 1, 1]})
        expected_lengths.index.name = 'index'
        pd.testing.assert_frame_equal(history.lengths, expected_lengths, 'Lengths attribute is not equal to the expected dataframe')

        assert history.number_of_systems == 3, 'Number of systems attribute is not None'
        assert history.columns.to_list() == ['time', 'event'], 'Columns attribute is not None'

        assert isinstance(history.indices, np.ndarray), 'Indices attribute is not an ndarray'
        np.testing.assert_array_equal(history.indices, np.array([0, 1, 2]), 'Indices attribute is not equal to the expected list')

        with pytest.raises(FileNotFoundError):
            History('invalid_filename.h5', verbose=False, chunksize=10000)

    def test_init_verbose_true(self):
        history = History(self.filename, chunksize=10000)
        assert history.verbose == False, 'Verbose attribute is not set correctly'


    def test_getitem_single_index(self):
        df = self.history[0]
        assert isinstance(df, pd.DataFrame), 'Returned object is not a DataFrame'
        assert len(df) == 1, 'Returned DataFrame does not have the correct length'


    def test_getitem_multiple_indices(self):
        df = self.history[[0, 1, 2]]
        assert isinstance(df, pd.DataFrame), 'Returned object is not a DataFrame'
        assert len(df) == 3, 'Returned DataFrame does not have the correct length'

    def test_getitem_index_array(self):
        indices = np.array([0, 1, 2])
        df = self.history[indices]
        assert isinstance(df, pd.DataFrame), 'Returned object is not a DataFrame'
        assert len(df) == 3, 'Returned DataFrame does not have the correct length'

    def test_getitem_single_column(self):
        column = 'time'
        df = self.history[column]
        assert isinstance(df, pd.DataFrame), 'Returned object is not a DataFrame'
        assert len(df.columns) == 1, 'Returned DataFrame does not have the correct number of columns'

    def test_getitem_invalid_column(self):
        column = 'invalid_column'
        with pytest.raises(ValueError):
            self.history[column]

    def test_getitem_invalid_keys(self):
        columns = ['time', 'invalid_column']
        with pytest.raises(ValueError):
            self.history[columns]

    def test_getitem_boolean_mask_numpy(self):
        mask = (self.history['time'] > 1).to_numpy()
        df = self.history[mask]
        assert isinstance(df, pd.DataFrame), 'Returned object is not a DataFrame'

    def test_getitem_boolean_mask_pandas(self):
        mask = self.history['time'] > 1
        df = self.history[mask]
        assert isinstance(df, pd.DataFrame), 'Returned object is not a DataFrame'

    def test_getitem_multiple_columns(self):
        columns = ['time', 'event']
        df = self.history[columns]
        assert isinstance(df, pd.DataFrame), 'Returned object is not a DataFrame'
        assert len(df.columns) == 2, 'Returned DataFrame does not have the correct number of columns'

    def test_getitem_invalid_key(self):
        with pytest.raises(ValueError):
            self.history[{1: 2}]

    def test_len(self):
        length = len(self.history)
        assert isinstance(length, int), 'Returned object is not an integer'
        assert length == 3, 'Returned length is not correct'

    def test_head(self):
        n = 2
        df = self.history.head(n)
        assert isinstance(df, pd.DataFrame), 'Returned object is not a DataFrame'
        assert len(df) == n, 'Returned DataFrame does not have the correct length'

    def test_tail(self):
        n = 2
        df = self.history.tail(n)
        assert isinstance(df, pd.DataFrame), 'Returned object is not a DataFrame'
        assert len(df) == n, 'Returned DataFrame does not have the correct length'

    def test_repr(self):
        representation = self.history.__repr__()
        assert isinstance(representation, str), 'Returned object is not a string'

    def test_repr_html(self):
        html_representation = self.history._repr_html_()
        assert isinstance(html_representation, str), 'Returned object is not a string'


    def test_select(self):
        df = self.history.select(where="time > 1", start=0, stop=10, columns=['event', 'time'])
        assert isinstance(df, pd.DataFrame), 'Returned object is not a DataFrame'
        assert len(df.columns) == 2, 'Returned DataFrame does not have the correct number of columns'
        assert len(df) == 2, 'Returned DataFrame does not have the correct length'



# Test the Oneline class
class TestOneline:

    @classmethod
    def setup_class(cls):
        # Set up a test HDF5 file using pandas HDFStore
        cls.filename = 'test_oneline.h5'
        with pd.HDFStore(cls.filename, 'w') as store:
            # Create a oneline dataframe
            oneline_data = pd.DataFrame({'time': [1, 2, 3], 'S1_mass_i': ['30','30', '70']})
            store.append('oneline', oneline_data, data_columns=True)

    def setup_method(self):
        self.oneline = Oneline(self.filename, verbose=False, chunksize=10000)

    def test_init(self):
        oneline = Oneline(self.filename, verbose=True, chunksize=5000)

        assert oneline.filename == self.filename, 'Filename is not set correctly'
        assert oneline.verbose == True, 'Verbose attribute is not set correctly'
        assert oneline.chunksize == 5000, 'Chunksize attribute is not set correctly'
        assert oneline.number_of_systems == 3, 'Number of systems attribute is not set correctly'
        assert oneline.columns.to_list() == ['time', 'S1_mass_i'], 'Columns attribute is not set correctly'
        assert oneline.number_of_systems == 3, 'Number of systems attribute is not set correctly'

        assert isinstance(oneline.indices, np.ndarray), 'Indices attribute is not an ndarray'
        np.testing.assert_array_equal(oneline.indices, np.array([0, 1, 2]), 'Indices attribute is not equal to the expected list')

        with pytest.raises(FileNotFoundError):
            Oneline('invalid_filename.h5', verbose=False, chunksize=10000)

    def test_getitem_single_index(self):
        df = self.oneline[0]
        assert isinstance(df, pd.DataFrame), 'Returned object is not a DataFrame'
        assert len(df) == 1, 'Returned DataFrame does not have the correct length'

    def test_getitem_multiple_indices(self):
        df = self.oneline[[0, 1, 2]]
        assert isinstance(df, pd.DataFrame), 'Returned object is not a DataFrame'
        assert len(df) == 3, 'Returned DataFrame does not have the correct length'

    def test_getitem_index_array(self):
        indices = np.array([0, 1, 2])
        df = self.oneline[indices]
        assert isinstance(df, pd.DataFrame), 'Returned object is not a DataFrame'
        assert len(df) == 3, 'Returned DataFrame does not have the correct length'
        pd.testing.assert_frame_equal(df,
                                      pd.DataFrame({'time': [1, 2, 3],
                                                    'S1_mass_i': ['30','30', '70']}),
                                      'Returned DataFrame is not equal to the expected DataFrame')

    def test_getitem_slice(self):
        df = self.oneline[0:2]
        assert isinstance(df, pd.DataFrame), 'Returned object is not a DataFrame'
        assert len(df) == 2, 'Returned DataFrame does not have the correct length'
        pd.testing.assert_frame_equal(df,
                                      pd.DataFrame({'time': [1, 2],
                                                    'S1_mass_i': ['30','30']}),)
    def test_getitem_endslice(self):
        df = self.oneline[:2]
        assert isinstance(df, pd.DataFrame), 'Returned object is not a DataFrame'
        assert len(df) == 2, 'Returned DataFrame does not have the correct length'
        pd.testing.assert_frame_equal(df,
                                      pd.DataFrame({'time': [1, 2],
                                                    'S1_mass_i': ['30','30']}),)
    def test_getitem_beginslice(self):
        df = self.oneline[1:]
        assert isinstance(df, pd.DataFrame), 'Returned object is not a DataFrame'
        assert len(df) == 2, 'Returned DataFrame does not have the correct length'
        pd.testing.assert_frame_equal(df,
                                      pd.DataFrame(index=[1, 2],
                                                   data={'time': [2, 3],
                                                         'S1_mass_i': ['30', '70']}),)
    def test_getitem_float_indices(self):
        with pytest.raises(ValueError):
            self.oneline[[0.5, 1.2]]



    def test_getitem_single_column(self):
        column = 'time'
        df = self.oneline[column]
        assert isinstance(df, pd.DataFrame), 'Returned object is not a DataFrame'
        assert len(df.columns) == 1, 'Returned DataFrame does not have the correct number of columns'

    def test_getitem_boolean_mask_numpy(self):
        mask = (self.oneline['time'] > 1).to_numpy().flatten()
        df = self.oneline[mask]
        assert isinstance(df, pd.DataFrame), 'Returned object is not a DataFrame'

    def test_getitem_boolean_mask_pandas(self):
        mask = self.oneline['time'] > 1
        print(mask)
        df = self.oneline[mask]
        assert isinstance(df, pd.DataFrame), 'Returned object is not a DataFrame'

    def test_getitem_multiple_columns(self):
        columns = ['time', 'S1_mass_i']
        df = self.oneline[columns]
        assert isinstance(df, pd.DataFrame), 'Returned object is not a DataFrame'
        assert len(df.columns) == 2, 'Returned DataFrame does not have the correct number of columns'

    def test_getitems_multiple_columns_invalid(self):
        columns = ['time', 'invalid_column']
        with pytest.raises(ValueError):
            self.oneline[columns]

    def test_getitem_invalid_key_type(self):
        with pytest.raises(ValueError):
            self.oneline[{1: 2}]

    def test_getitem_invalid_key(self):
        with pytest.raises(ValueError):
            self.oneline['invalid_key']

    def test_len(self):
        length = len(self.oneline)
        assert isinstance(length, int), 'Returned object is not an integer'
        assert length == 3, 'Returned length is not correct'

    def test_head(self):
        n = 2
        df = self.oneline.head(n)
        assert isinstance(df, pd.DataFrame), 'Returned object is not a DataFrame'
        assert len(df) == n, 'Returned DataFrame does not have the correct length'

    def test_tail(self):
        n = 2
        df = self.oneline.tail(n)
        assert isinstance(df, pd.DataFrame), 'Returned object is not a DataFrame'
        assert len(df) == n, 'Returned DataFrame does not have the correct length'

    def test_repr(self):
        representation = self.oneline.__repr__()
        assert isinstance(representation, str), 'Returned object is not a string'

    def test_repr_html(self):
        html_representation = self.oneline._repr_html_()
        assert isinstance(html_representation, str), 'Returned object is not a string'

    def test_select(self):
        df = self.oneline.select(where="time > 1", start=0, stop=10, columns=['S1_mass_i', 'time'])
        assert isinstance(df, pd.DataFrame), 'Returned object is not a DataFrame'
        assert len(df.columns) == 2, 'Returned DataFrame does not have the correct number of columns'
        assert len(df) == 2, 'Returned DataFrame does not have the correct length'


    @classmethod
    def teardown_class(cls):
        os.remove(cls.filename)

# Test the PopulationIO class
class TestPopulationIO:

    def setup_method(self):
        self.filename = "test_population.h5"

    def teardown_method(self):
        if os.path.exists(self.filename):
            os.remove(self.filename)

    def test_init(self):
        pop_io = PopulationIO()
        assert pop_io.verbose == False, "Verbose attribute is not set correctly"


    def test_invalid_filename(self):
        pop_io = PopulationIO()
        with pytest.raises(ValueError):
            pop_io._load_metadata("invalid_filename")

    def test_save_and_load_mass_per_met(self):
        population_io = PopulationIO()
        population_io.verbose = True
        population_io.mass_per_metallicity = pd.DataFrame({"metallicity": [0.02, 0.04], "mass": [1.0, 2.0]})
        population_io._save_mass_per_metallicity(self.filename)

        loaded_io = PopulationIO()
        loaded_io.verbose = True
        loaded_io._load_mass_per_metallicity(self.filename)
        pd.testing.assert_frame_equal(population_io.mass_per_metallicity, loaded_io.mass_per_metallicity)

    def test_save_and_load_ini_params(self):
        population_io = PopulationIO()
        population_io.ini_params = {i:10 for i in parameter_array}
        population_io._save_ini_params(self.filename)

        loaded_io = PopulationIO()
        loaded_io._load_ini_params(self.filename)

        assert population_io.ini_params == loaded_io.ini_params, "Loaded ini_params are not equal to the saved ini_params"

    def test_save_and_load_metadata(self):
        pop_io = PopulationIO()
        pop_io.verbose = True
        pop_io.ini_params = {i:10 for i in parameter_array}
        pop_io._save_ini_params(self.filename)
        pop_io.mass_per_metallicity = pd.DataFrame({"metallicity": [0.02, 0.04], "mass": [1.0, 2.0]})
        pop_io._save_mass_per_metallicity(self.filename)

        load_io = PopulationIO()
        load_io._load_metadata(self.filename)

        assert pop_io.ini_params == load_io.ini_params, "Loaded ini_params are not equal to the saved ini_params"
        assert pop_io.mass_per_metallicity.equals(load_io.mass_per_metallicity), "Loaded mass_per_metallicity is not equal to the saved mass_per_metallicity"



class TestPopulation:
    def setup_method(self):
        pass

    def teardown_method(self):
        # Clean up any resources used by the test
        pass

    def setup_class(self):
        self.filename1 = "no_mass_per_met_population.h5"
        self.filename2 = "history_population.h5"
        self.filename3 = "oneline_population.h5"
        self.history_data = pd.DataFrame({'time': [1, 2, 3], 'event': ['ZAMS','oRLO1', 'CEE']})
        self.oneline_data = pd.DataFrame({'time': [1, 2, 3], 'S1_mass_i': [30, 30, 70], 'S2_mass_i': [30, 30, 70.]})
        self.formation_channels = pd.DataFrame({'channel': ['channel1', 'channel2', 'channel3'], 'channel_debug':['debug1', 'debug2', 'debug3']})

        # create a file with only history and oneline data
        with pd.HDFStore(self.filename1, 'w') as store:
            store.append('history',self.history_data, data_columns=True)
            store.append('oneline', self.oneline_data, data_columns=True)

        with pd.HDFStore(self.filename2, 'w') as store:
            store.append('history', self.history_data, data_columns=True)

        with pd.HDFStore(self.filename3, 'w') as store:
            store.append('oneline', self.oneline_data, data_columns=True)

    def teardown_class(self):
        if os.path.exists(self.filename1):
            os.remove(self.filename1)
        if os.path.exists(self.filename2):
            os.remove(self.filename2)
        if os.path.exists(self.filename3):
            os.remove(self.filename3)

    @pytest.fixture
    def mass_per_met_pop(self):
        self.filename = "mass_per_met_population.h5"
        with pd.HDFStore(self.filename, 'w') as store:
            store.append('history', self.history_data, data_columns=True)
            store.append('oneline', self.oneline_data, data_columns=True)

        pop = Population(self.filename, verbose=True, metallicity=0.02, ini_file=PATH_TO_POSYDON+'/posydon/popsyn/population_params_default.ini')
        yield
        if os.path.exists(self.filename):
            os.remove(self.filename)

    @pytest.fixture
    def no_mass_per_met_pop(self):
        self.filename = "no_mass_per_met_population.h5"
        with pd.HDFStore(self.filename, 'w') as store:
            store.append('history', self.history_data, data_columns=True)
            store.append('oneline', self.oneline_data, data_columns=True)
        yield
        if os.path.exists(self.filename):
            os.remove(self.filename)

    @pytest.fixture
    def mass_per_met_pop_channels(self):
        self.filename = "mass_per_met_population.h5"
        with pd.HDFStore(self.filename, 'w') as store:
            store.append('history', self.history_data, data_columns=True)
            store.append('oneline', self.oneline_data, data_columns=True)
            store.append('formation_channels', self.formation_channels, data_columns=True)
        pop = Population(self.filename, verbose=True, metallicity=0.02, ini_file=PATH_TO_POSYDON+'/posydon/popsyn/population_params_default.ini')
        yield
        if os.path.exists(self.filename):
            os.remove(self.filename)

    @pytest.fixture
    def clean_up_selection_file(self):
        self.outfile = "test_selection.h5"
        yield
        if os.path.exists(self.outfile):
            os.remove(self.outfile)

    def test_init_invalid_file(self):
        with pytest.raises(ValueError):
            pop = Population('invalid_filename')

    def test_init_no_history(self):
        with pytest.raises(ValueError):
            pop = Population(self.filename3)

    def test_init_no_oneline(self):
        with pytest.raises(ValueError):
            pop = Population(self.filename2)

    def test_init_no_mass_per_met(self):
        with pytest.raises(ValueError):
            pop = Population(self.filename1, verbose=True)


    def test_init_mass_per_met_calc(self, no_mass_per_met_pop: None):
        pop = Population(self.filename, verbose=True, metallicity=1., ini_file=PATH_TO_POSYDON+'/posydon/popsyn/population_params_default.ini')
        # check that the history and oneline data are read correctly
        pd.testing.assert_frame_equal(pop.history[:], self.history_data)
        pd.testing.assert_frame_equal(pop.oneline[:], self.oneline_data)
        assert pop.solar_metallicities == [1.]
        assert pop.metallicities == [1*Zsun]
        pd.testing.assert_frame_equal(pop.mass_per_metallicity, pd.DataFrame(index=[1.], data={'simulated_mass': [260.], 'underlying_mass': [1462.194834], 'number_of_systems': [3]}))

        pop = Population(self.filename, verbose=True, metallicity=1., ini_file=PATH_TO_POSYDON+'/posydon/popsyn/population_params_default.ini')
        # check that the history and oneline data are the same
        pd.testing.assert_frame_equal(pop.history[:], self.history_data)
        pd.testing.assert_frame_equal(pop.oneline[:], self.oneline_data)
        assert pop.solar_metallicities == [1.]
        assert pop.metallicities == [1*Zsun]
        pd.testing.assert_frame_equal(pop.mass_per_metallicity, pd.DataFrame(index=[1.], data={'simulated_mass': [260.], 'underlying_mass': [1462.194834], 'number_of_systems': [3]}))


    def test_init(self,mass_per_met_pop: None):
        pop = Population(self.filename)
        # check that the history and oneline data are read correctly
        pd.testing.assert_frame_equal(pop.history[:], self.history_data)
        pd.testing.assert_frame_equal(pop.oneline[:], self.oneline_data)
        assert pop.metallicities == [0.02*Zsun]
        assert pop.solar_metallicities == [0.02]
        tmp_df = pd.DataFrame(index=[0, 1, 2], data={'index': [1, 1, 1]})
        tmp_df.index.name = 'index'
        pd.testing.assert_frame_equal(pop.history_lengths, tmp_df)
        pd.testing.assert_frame_equal(pop.mass_per_metallicity, pd.DataFrame(index=[0.02], data={'simulated_mass': [260.], 'underlying_mass': [1462.194834], 'number_of_systems': [3]}))


    def test_read_formation_channels(self, mass_per_met_pop_channels: None):
        pop = Population(self.filename)
        # check that the formation channels are read correctly
        assert pop.formation_channels.equals(self.formation_channels)

    def test_export_selection(self, mass_per_met_pop: None, clean_up_selection_file: None):
        selection = [1, 2]
        chunksize = 1000
        pop = Population(self.filename)
        pop.export_selection(selection, self.outfile, chunksize)
        assert os.path.exists(self.outfile)
        assert pd.read_hdf(self.outfile, 'history').shape[0] == 2
        assert pd.read_hdf(self.outfile, 'oneline').shape[0] == 2

    def test_bad_name_export_selection(self, mass_per_met_pop: None):
        selection = [1, 2]
        chunksize = 1000
        pop = Population(self.filename)
        with pytest.raises(ValueError):
            pop.export_selection(selection, 'test_selection.csv', history_chunksize=chunksize)

    def test_append_selection(self, mass_per_met_pop: None, clean_up_selection_file: None):
        selection = [1, 2]
        chunksize = 1000
        pop = Population(self.filename)
        pop.export_selection(selection, self.outfile, overwrite=True, history_chunksize=chunksize)
        pop.export_selection(selection, self.outfile, overwrite=False, history_chunksize=chunksize)

        assert pd.read_hdf(self.outfile, 'history').shape[0] == 4
        assert pd.read_hdf(self.outfile, 'oneline').shape[0] == 4

    def test_no_formation_channels(self, mass_per_met_pop: None):
        pop = Population(self.filename, verbose=True)
        assert pop.formation_channels is None

    def test_len(self, mass_per_met_pop: None):
        pop = Population(self.filename)
        assert len(pop) == 3

    def test_columns(self, mass_per_met_pop: None):
        pop = Population(self.filename)
        columns = pop.columns

        assert columns['history'].tolist() == self.history_data.columns.tolist()
        assert columns['oneline'].tolist() ==  self.oneline_data.columns.tolist()




    # Test formation channel calculation, I need a specific test file for this,
    # since it requires specific columns to be present in the oneline and history dataframes

    # Test create_transient_population method requires a specific test file for this,
    # since it requires specific columns to be present in the oneline and history dataframes


class TestTransientPopulation:
    pass
    # to implement



class TestRates:
    pass
    # to implement

# Run the tests

if __name__ == '__main__':
    pytest.main()
