import unittest
import os
from posydon.popsyn.synthetic_population import PopulationRunner, ParsedPopAttrs, ParsedPopulation



# test the PopulationRunner class

class TestPopulationRunner(unittest.TestCase):
    
    # Test the initialisation of the PopulationRunner class
    def test_init(self):
        # Test the initialisation of the PopulationRunner class
        poprun = PopulationRunner('../../popsyn/population_params_default.ini')

        # Parameters are based on the default values in the ini file
        assert poprun.synthetic_pop_params['star_formation'] == 'burst', 'Default star formation is not correct'
        assert poprun.synthetic_pop_params['number_of_binaries'] == 10, 'Default number of binaries is not correct'
        assert poprun.synthetic_pop_params['binary_fraction_scheme'] == 'const', 'Default binary fraction scheme is not correct'
        assert poprun.metallicities == [0.0001], 'Default metallicities are not correct'
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
        poprun.metallicities = [0.0001, 0.001]
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
        for met in poprun.metallicities:
            folder_name = poprun.create_met_prefix(met) + poprun.synthetic_pop_params['temp_directory']
            assert os.path.isdir(folder_name), 'Batch folder has not been created'                    


    # test the merge parallel runs
    def test_merge_parallel_runs(self):
        poprun = PopulationRunner('../../popsyn/population_params_default.ini')
        paths_to_batches = []
        for met in poprun.metallicities:
            paths_to_batches.append(poprun.create_met_prefix(met) + poprun.synthetic_pop_params['temp_directory'])
        
        poprun.merge_parallel_runs(paths_to_batches)
        # Check that the population has been created
        for met in poprun.metallicities:
            file_name = poprun.create_met_prefix(met) + 'population.h5'
            assert os.path.isfile(file_name), 'Population file has not been created'
        
    
    @classmethod
    def tearDownClass(cls):
        # remove files created by the tests
        poprun = PopulationRunner('../../popsyn/population_params_default.ini')
        for met in poprun.metallicities:
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
        pass
    
    
    
if __name__ == '__main__':
    unittest.main()
    