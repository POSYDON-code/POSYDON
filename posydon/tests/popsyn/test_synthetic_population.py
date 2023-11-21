import unittest
import os
import numpy as np
from posydon.popsyn.synthetic_population import PopulationRunner, ParsedPopAttrs, ParsedPopulation, SyntheticPopulation

from posydon.config import PATH_TO_POSYDON_DATA


# test the PopulationRunner class
@unittest.skip
class TestPopulationRunner(unittest.TestCase):
    
    # Test the initialisation of the PopulationRunner class
    def test_init(self):
        # Test the initialisation of the PopulationRunner class
        poprun = PopulationRunner('../../popsyn/population_params_default.ini')

        # Parameters are based on the default values in the ini file
        assert poprun.synthetic_pop_params['star_formation'] == 'burst', 'Default star formation is not correct'
        assert poprun.synthetic_pop_params['number_of_binaries'] == 10, 'Default number of binaries is not correct'
        assert poprun.synthetic_pop_params['binary_fraction_scheme'] == 'const', 'Default binary fraction scheme is not correct'
        assert np.allclose(poprun.solar_metallicities, [0.0001]), 'Default metallicities are not correct'
        assert poprun.binary_populations == None, 'Default binary populations are not correctly initialised'

    # Test the create_population method
    def test_create_population(self):
        # Test the create_population method
        poprun = PopulationRunner('../../popsyn/population_params_default.ini')
        poprun.create_binary_populations()
        assert poprun.binary_populations != None, 'Binary populations have not been created'
        assert len(poprun.binary_populations) == 1, 'The binary_populations list is not the correct length'
        
    def test_create_population_multi_met(self):
        poprun = PopulationRunner('../../popsyn/population_params_default.ini')
        # add additional metallicity manually 
        poprun.solar_metallicities = [0.0001, 0.001]
        poprun.create_binary_populations()
        assert poprun.binary_populations != None, 'Binary populations have not been created'
        assert len(poprun.binary_populations) == 2, 'The binary_populations list is not the correct length'
    
    
    # Test the evolve
    def test_evolve_pop(self):
        # This only test the population creation for a single process.
        poprun = PopulationRunner('../../popsyn/population_params_default.ini')
        poprun.evolve()
        # Check that the population has been created
        # check the batch folder
        for met in poprun.solar_metallicities:
            folder_name = poprun.create_met_prefix(met) + poprun.synthetic_pop_params['temp_directory']
            assert os.path.isdir(folder_name), 'Batch folder has not been created'                    


    # test the merge parallel runs
    def test_merge_parallel_runs(self):
        poprun = PopulationRunner('../../popsyn/population_params_default.ini')
        paths_to_batches = []
        for met in poprun.solar_metallicities:
            paths_to_batches.append(poprun.create_met_prefix(met) + poprun.synthetic_pop_params['temp_directory'])
        
        poprun.merge_parallel_runs(paths_to_batches)
        # Check that the population has been created
        for met in poprun.solar_metallicities:
            file_name = poprun.create_met_prefix(met) + 'population.h5'
            assert os.path.isfile(file_name), 'Population file has not been created'
        
    
    @classmethod
    def tearDownClass(cls):
        # remove files created by the tests
        poprun = PopulationRunner('../../popsyn/population_params_default.ini')
        for met in poprun.solar_metallicities:
            file_name = poprun.create_met_prefix(met) + 'population.h5'
            if os.path.isfile(file_name):
                os.remove(file_name)
            batch_folder = poprun.create_met_prefix(met) + poprun.synthetic_pop_params['temp_directory'] 
            # remove files + folders if a test failed
            if os.path.isdir(batch_folder):
                os.remove(batch_folder+'/evolution.combined')
                os.rmdir(batch_folder)
    

# Test the ParsedPopAttrs class
class TestParsedPopAttrs(unittest.TestCase):
    
    # Test the initialisation of the ParsedPopAttrs class
    # Are all variables present?
    def test_ini(self):
        pop = ParsedPopAttrs()
        assert pop.underlying_mass_per_met == None, 'Underlying mass per metallicity is not None'
        assert pop.simulated_mass_per_met == None, 'Simulated mass per metallicity is not None'
        assert pop.metallicities == None, 'Metallicities are not None'
        assert pop.solar_metallicities == None, 'Solar metallicities are not None'
        assert pop.verbose == False, 'Verbose is not False'
        assert pop.synthetic_pop_params == None, 'Synthetic population parameters are not None'
    



class TestParsedPop(unittest.TestCase):
    
    # Test the initialisation of the ParsedPopulation class
    def test_ini(self):
        parsed_pop = ParsedPopulation('testing.h5',
                                      path_to_ini='../../popsyn/population_params_default.ini',
                                      verbose=True)
        
        # check metallicities
        assert np.allclose(parsed_pop.solar_metallicities,[0.0001]), 'Default solar metallicities are not correct'
        assert np.isclose(parsed_pop.metallicities[0], 1.42e-6), 'Absolute metallicity is incorrect'
        # check ini_params
        assert parsed_pop.ini_params != None, 'Ini params are not saved'
        parameter_array = [ "number_of_binaries",
                   'binary_fraction_scheme',
                   'binary_fraction_const',
                   'star_formation',
                   'max_simulation_time',
                   'primary_mass_scheme',
                   'primary_mass_min',                              
                   'primary_mass_max',                                  
                   'secondary_mass_scheme',
                   'secondary_mass_min',
                   'secondary_mass_max',
                   'orbital_scheme',
                   'orbital_period_scheme',
                   'orbital_period_min',
                   'orbital_period_max',
                   'eccentricity_scheme']
        for param in parameter_array:
            assert parsed_pop.ini_params[param] != None, f'{param} is not saved'
        assert parsed_pop.verbose == True, 'Verbose is not True'
    
    # load in a population file and parse BBH
    def test_parse_and_save_bbh(self):
        parsed_pop = ParsedPopulation('testing.h5',
                                      path_to_ini='../../popsyn/population_params_default.ini',)
        # TODO: remove hard coded link!
        path_to_data = os.path.join(PATH_TO_POSYDON_DATA, "POSYDON_data/tutorials/population-synthesis/example/")
        path_to_data += '/1.00e+00_Zsun_population.h5'
        
        parsed_pop.solar_metallicities = [1.00]
        parsed_pop.metallicities = [0.0142]
        
        parsed_pop.parse(path_to_data=[path_to_data],
                         S1_state='BH',
                         S2_state='BH',
                         binary_state='contact',
                         invert_S1S2=False,
                         chunksize=5000000)

        assert parsed_pop.total_systems[0] == 1000000, 'Parsed population is not the correct length'
        assert parsed_pop.selected_systems[0] == 233, 'Parsed population is not the correct length'
        # mass parameters
        assert np.all(parsed_pop.underlying_mass_per_met != None), 'Underlying mass per metallicity is None'
        assert np.all(parsed_pop.simulated_mass_per_met != None), 'Simulated mass per metallicity is None'
        assert len(parsed_pop.underlying_mass_per_met) == 1, 'Underlying mass per metallicity is not the correct length'
        # path to data
        assert parsed_pop.path_to_data[0] == path_to_data, f'{path_to_data} is not equal to {parsed_pop.path_to_data[0]}'
        # parse kwargs data
        assert parsed_pop.parse_kwargs != None, 'Parse kwargs are not saved'
        assert parsed_pop.parse_kwargs['S1_state'] == 'BH', 'S1_state is not saved'
        assert parsed_pop.parse_kwargs['S2_state'] == 'BH', 'S2_state is not saved'
        assert parsed_pop.parse_kwargs['binary_state'] == 'contact', 'binary_state is not saved'
        assert parsed_pop.parse_kwargs['binary_event'] == None, 'binary_event is not saved'
        assert parsed_pop.parse_kwargs['step_name'] == None, 'step_name is not saved'
        assert parsed_pop.parse_kwargs['invert_S1S2'] == False, 'invert_S1S2 is not saved'
        
        pp = ParsedPopulation('testing.h5')
        
        assert pp.parse_kwargs != None, 'Parse kwargs are not saved'
        assert len(pp.underlying_mass_per_met) == 1, 'Underlying mass per metallicity is None'
        assert len(pp.simulated_mass_per_met) == 1, 'Simulated mass per metallicity is None'
        assert pp.path_to_data[0] == path_to_data, 'Path to data is not saved'
        assert pp.total_systems[0] == 1000000, 'total_systems is not saved'
        assert pp.selected_systems[0] == 233, 'selected_systems is not saved'

    
        
    # def test_get_formation_channel(self):
    #     parsed_pop = ParsedPopulation('../../popsyn/population_params_default.ini')
    #     parsed_pop.load('test.h5')
    #     pre_shape = parsed_pop.df_oneline.shape
    #     parsed_pop.get_formation_channels(mt_history=True)
        
    #     assert pre_shape[0] == parsed_pop.df_oneline.shape[0], 'df_oneline has changed nr of rows'
    #     assert pre_shape[1] == parsed_pop.df_oneline.shape[1]-2, 'df_oneline does not have the additional columns'
    #     assert any(parsed_pop.df_oneline['channel']), 'channels have not been added'
    #     assert any(parsed_pop.df_oneline['channel_debug']), 'debug channels have not been added'
    #     assert parsed_pop.df_oneline['channel'].iloc[0] == 'ZAMS_oRLO1-contact_CC1_oRLO2_CC2_END', 'channel is not correct'
        
    # def test_write_read_formation_channels(self):
    #     parsed_pop = ParsedPopulation('../../popsyn/population_params_default.ini')
    #     parsed_pop.load('test.h5')
    #     parsed_pop.get_formation_channels(mt_history=True)
    #     parsed_pop.save('mt_test.h5')
    #     self.assertRaises(ValueError, parsed_pop.load, 'mt_test.h5')
    #     parsed_pop = ParsedPopulation('../../popsyn/population_params_default.ini')
    #     parsed_pop.load('mt_test.h5')
    #     assert any(parsed_pop.df_oneline['channel']), 'channels have not been added'
    #     assert any(parsed_pop.df_oneline['channel_debug']), 'debug channels have not been added'
    #     assert parsed_pop.df_oneline['channel'].iloc[0] == 'ZAMS_oRLO1-contact_CC1_oRLO2_CC2_END', 'channel is not correct'
        
    # def test__create_DCO_population(self):
        
    #     parsed_pop = ParsedPopulation('../../popsyn/population_params_default.ini')
    #     parsed_pop.load('mt_test.h5')
    #     DCO_population = parsed_pop.create_DCO_population()
        
    #     assert type(DCO_population) == SyntheticPopulation
        
        
        
    
if __name__ == '__main__':
    unittest.main()
    