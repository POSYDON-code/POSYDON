"""Evolve multiple BinaryPopulations together.

e.g. with multiple solar_metallicities
"""

__authors__ = [
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Kyle Akira Rocha <kylerocha2024@u.northwestern.edu>",
    "Monica Gallegos-Garcia <monicagallegosgarcia2024@u.northwestern.edu>",
    "Max Briel < max.briel@gmail.com",
]

import warnings
import numpy as np
import pandas as pd
from tqdm import tqdm
import os
import copy
from matplotlib import pyplot as plt
import multiprocessing as mp

from posydon.utils.constants import Zsun
from posydon.popsyn.io import binarypop_kwargs_from_ini
from posydon.popsyn.binarypopulation import BinaryPopulation
from posydon.utils.common_functions import convert_metallicity_to_string
from posydon.popsyn.normalized_pop_mass import initial_total_underlying_mass
from posydon.popsyn.rate_calculation import Rates
import posydon.visualization.plot_pop as plot_pop
from posydon.popsyn.GRB import get_GRB_properties, GRB_PROPERTIES

# TODO: temp import, remove after TF2 classification is implemented in pop synth
from posydon.interpolation.IF_interpolation import IFInterpolator
from posydon.binary_evol.binarystar import BinaryStar
from posydon.binary_evol.singlestar import SingleStar
from posydon.utils.common_functions import convert_metallicity_to_string

from posydon.popsyn.binarypopulation import HISTORY_MIN_ITEMSIZE, ONELINE_MIN_ITEMSIZE

class PopulationRunner:
    
    def __init__(self, path_to_ini, verbose=False):
        '''A Synthetic Population to run, process, and parse BinaryPopulations.
        
        Parameters
        -----------
        path_to_ini : str
            path to the ini file to parse, which determines the populations.
        verbose : bool (optional; True)
            Whether to print progress statements.
        '''
        self.synthetic_pop_params = None
        self.solar_metallicities = None
        self.binary_populations = None
        self.verbose = verbose
        
        if '.ini' not in path_to_ini:
            raise ValueError('You did not provide a valid path_to_ini!')
        else:
            self.synthetic_pop_params = binarypop_kwargs_from_ini(path_to_ini)
            self.solar_metallicities = self.synthetic_pop_params['metallicity']
            if not isinstance( self.solar_metallicities, list):
                self.solar_metallicities = [self.solar_metallicities]
            self.binary_populations = None

    def create_binary_populations(self):
        """Create a list of BinaryPopulation objects."""
        self.binary_populations = []
        for met in self.solar_metallicities[::-1]:
            ini_kw = copy.deepcopy(self.synthetic_pop_params)
            ini_kw['metallicity'] = met
            ini_kw['temp_directory'] = self.create_met_prefix(met) + self.synthetic_pop_params['temp_directory']
            self.binary_populations.append(BinaryPopulation(**ini_kw))  
            
    def get_ini_kw(self):
        return self.synthetic_pop_params.copy()

    def evolve(self):
        """Evolve population(s) at given Z(s)."""
        if self.binary_populations is None:
            self.create_binary_populations()
        while self.binary_populations:
            pop =  self.binary_populations.pop()
            if self.verbose:
                print( f'Z={pop.kwargs["metallicity"]:.2e} Z_sun' )
            pop.evolve()
            del pop

    def merge_parallel_runs(self, path_to_batches):
        """
        Merge the folder or list of folders into a single file per metallicity.

        Parameters
        ----------
        path_to_batches : str or list of str
            Path to the folder(s) containing the batch folders.
        """
        # Max Briel: Should this function read in the ini file from the batch?
        
        if isinstance(path_to_batches, str):
            path_to_batches = [path_to_batches]
        # check if path_to_batches is the same length as the number of solar_metallicities
        if len(path_to_batches) != len(self.solar_metallicities):
            raise ValueError('The number of metallicity and batch directories do not match!')

        for met, path_to_batch in zip(self.solar_metallicities, path_to_batches):
            met_prefix = self.create_met_prefix(met)
            ini_kw = self.synthetic_pop_params.copy()
            ini_kw['metallicity'] = met
            tmp_files = [os.path.join(path_to_batch, f)     \
                         for f in os.listdir(path_to_batch) \
                            if os.path.isfile(os.path.join(path_to_batch, f))]
            
            BinaryPopulation(**ini_kw).combine_saved_files(met_prefix+ 'population.h5', tmp_files)
            if self.verbose:
                print(f'Population at Z={met:.2e} Z_sun successfully merged!')
            # Store the population ini parameters inside the file.
            # This is useful to keep track of the population parameters
            # when merging multiple populations.
            with pd.HDFStore(met_prefix+ 'population.h5', mode='a') as store:
                store.put('ini_parameters', pd.Series(ini_kw))
            
            if len(os.listdir(path_to_batch)) == 0:
                os.rmdir(path_to_batch)
            elif self.verbose:
                print(f'{path_to_batch} is not empty, it was not removed!')

    @staticmethod
    def create_met_prefix(met):
        """Append a prefix to the name of directories for batch saving."""
        return convert_metallicity_to_string(met) + '_Zsun_'
    
###############################################################################

class ParsedPopAttrs:
    
    def __init__(self):
        self.verbose = False
        self.parse_kwargs = None
        self.ini_params = None
        self.solar_metallicities = None
        self.metallicities = None
        self.underlying_mass_per_met = None
        self.simulated_mass_per_met = None

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

class ParsedPopulationIO():

    def __init__(self, path_to_parsed_pop_file=None, verbose=False):
        self.verbose = verbose
        self.parsed_pop_file = path_to_parsed_pop_file    
    
    def load_params_from_ini(self, parsed_pop_instance,  path_to_ini):
        '''load parameters from ini file'''
        synthetic_pop_params = binarypop_kwargs_from_ini(path_to_ini)
        parsed_pop_instance.solar_metallicities = np.array(synthetic_pop_params['metallicity'])
        parsed_pop_instance.metallicities = parsed_pop_instance.solar_metallicities * Zsun
        parsed_pop_instance.ini_params = {}
        for c in parameter_array:
            parsed_pop_instance.ini_params[c] = synthetic_pop_params[c]

    def save_ini_params(self, parsed_pop_instance):
        '''
        Store the ini parameters to the ParsedPopulation file'''
        with pd.HDFStore(self.parsed_pop_file, mode='a') as store:
            
            # write ini parameters to file
            tmp_df = pd.DataFrame()
            for c in parameter_array:
                tmp_df[c] = [parsed_pop_instance.ini_params[c]]
            store.put('ini_parameters', tmp_df)
            
            # write metallicity to file
            tmp_df = pd.DataFrame(index=parsed_pop_instance.metallicities)
            store.put('mass_per_met', tmp_df)
            
    def load_ini_params(self, parsed_pop_instance):
        # load ini parameters
        with pd.HDFStore(self.parsed_pop_file, mode='r', ) as store:
            tmp_df = store['ini_parameters']
            parsed_pop_instance.ini_params = {}
            for c in parameter_array:
                parsed_pop_instance.ini_params[c] = tmp_df[c][0]
            
            parsed_pop_instance.metallicities = store['mass_per_met'].index.to_numpy()
        
        # extract solar metallicities from metallicity
        parsed_pop_instance.solar_metallicities = parsed_pop_instance.metallicities / Zsun    

    def save_parse_kwargs(self, parsed_pop_instance):
        with pd.HDFStore(self.parsed_pop_file, mode='a') as store:
            tmp_df = pd.DataFrame()
            for c in parsed_pop_instance.parse_kwargs.keys():
                tmp_df[c] = [parsed_pop_instance.parse_kwargs[c]]
            
            store.put('parse_kwargs', tmp_df.astype(str))
            
    def load_parse_kwargs(self, parsed_pop_instance):
        """load the parse arguments used to create the population
        """    
        with pd.HDFStore(self.parsed_pop_file, mode='r') as store:
            parsed_pop_instance.parse_kwargs = store['parse_kwargs'].loc[0].to_dict()

    def save_per_metallicity_info(self, parsed_pop_instance):
        with pd.HDFStore(self.parsed_pop_file, mode='a') as store:
            tmp_df = pd.DataFrame(index=parsed_pop_instance.metallicities)
            tmp_df['underlying_mass'] = parsed_pop_instance.underlying_mass_per_met
            tmp_df['simulated_mass'] = parsed_pop_instance.simulated_mass_per_met
            tmp_df['selected_systems'] = parsed_pop_instance.selected_systems
            tmp_df['total_systems'] = parsed_pop_instance.total_systems
            store.put('mass_per_met', tmp_df)
    
    def save_path_to_data(self, parsed_pop_instance):
        with pd.HDFStore(self.parsed_pop_file, mode='a') as store:
            tmp_df = pd.Series(parsed_pop_instance.path_to_data)
            store.put('path_to_data', tmp_df)
        
    def load_per_metallicity_info(self, parsed_pop_instance):
        """Load the underlying mass, simulated mass, selected systems, 
        and total_systems"""
        with pd.HDFStore(self.parsed_pop_file, mode='r') as store:
            tmp_df = store['mass_per_met']
            parsed_pop_instance.underlying_mass_per_met = tmp_df['underlying_mass'].to_numpy()
            parsed_pop_instance.simulated_mass_per_met = tmp_df['simulated_mass'].to_numpy()
            parsed_pop_instance.selected_systems = tmp_df['selected_systems'].to_numpy()
            parsed_pop_instance.total_systems = tmp_df['total_systems'].to_numpy()
            
    def load_path_to_data(self, parsed_pop_instance):
        """Load the path_to_data array
        
        Contains which population files were used to create the
        Parsed Population file"""
        with pd.HDFStore(self.parsed_pop_file, model='r') as store:
            parsed_pop_instance.path_to_data = store['path_to_data'].to_numpy()

    def read_formation_channels(self, parsed_pop_instance):
        '''Read the formation channels from the ParsedPopulation file'''
        with pd.HDFStore(self.parsed_pop_file, mode='r') as store:
            tmp_df = store['formation_channels']
        return tmp_df
    
    def write_formation_channels(self, parsed_pop_instance, df):
        '''Write the formation channels to the ParsedPopulation file'''
        mask = df['channel_debug'].str.len() > 100
        if any(mask):
            if self.verbose:
                print('Warning: Channels are longer than 100 characters:')
                print('This should not be possible and there is likely a bug in the flow_chart')
                print('We will truncate the channels to 100 characters!')
            df.loc[mask, 'channel_debug'] = df.loc[mask, 'channel_debug'].str[:100]
            df.loc[mask, 'channel'] = df.loc[mask, 'channel'].str[:100]

        with pd.HDFStore(self.parsed_pop_file, mode='a') as store:
            store.append('formation_channels',
                      df,
                      format='table',
                      data_columns=True,
                      min_itemsize = {'channel': 100, 'channel_debug': 100})
            
    def load_parsed_pop_path(self, synthetic_pop_instance):
        '''Load the path to the parsed population file'''
        with pd.HDFStore(self.parsed_pop_file, model='r') as store:
            synthetic_pop_instance.parsed_pop_path = store['parsed_pop_path'].to_numpy()[0]
            
    def save_parsed_pop_path(self, parsed_pop_path):
        '''Save the path to the parsed population file'''
        with pd.HDFStore(self.parsed_pop_file, mode='a') as store:
            store.put('parsed_pop_path', pd.Series(parsed_pop_path))
    
class ParsedPopulation():
    
    def __init__(self, path_to_parsed_pop_file, path_to_ini=None, verbose=False, chunksize=500000):
        '''A parsed stellar population from raw BinaryPopulations.
        
        Calculates the total number of binaries in the population, and the
        total mass of the simulated and underlying stellar populations.
        
        The user has to make sure that the ini file and population files
        are consistent!
        
        Parameters
        -----------
        
        path_to_data : str or list of str
            Path to the data to parse.
            
        '''
        # initialise the parsed population class variables
        
        ParsedPopAttrs.__init__(self)
    
        self.verbose = verbose
        self.chunksize = chunksize
        # ini file given, create parsed_pop_file
        if path_to_ini is not None:
            # Read path_to_ini
            if '.ini' not in path_to_ini:
                raise ValueError('You did not provide a valid path_to_ini!')
            else:
                ParsedPopulationIO().load_params_from_ini(self,path_to_ini)
            
            # makes sure the parsed population isn't overwritten or double added
            if os.path.isfile(path_to_parsed_pop_file):
                raise FileExistsError('The parsed population file already exists!')                
            
            # create parsed_pop_file
            self.parsed_pop_file = path_to_parsed_pop_file
            self.io = ParsedPopulationIO(self.parsed_pop_file, self.verbose)
            self.io.save_ini_params(self)
            
        # if only parsed_pop_file is given, load the data from the file
        else:
            self.parsed_pop_file = path_to_parsed_pop_file
            self.io = ParsedPopulationIO(self.parsed_pop_file, self.verbose)
            self.io.load_ini_params(self)
            self.io.load_parse_kwargs(self)
            self.io.load_per_metallicity_info(self)
            self.io.load_path_to_data(self)
        
    def parse(self, path_to_data, S1_state=None, S2_state=None, binary_state=None,
               binary_event=None, step_name=None, invert_S1S2=False):
        """Given the path to the data, parse the datafiles for binaries which
        satisfy the given conditions.
        
        It also stores the underlying stellar mass and
        the initial simulated stellar mass for each metallicity.

        Parameters
        ----------
        S1_state : str
            Star1 stellar state.
        S2_state : str
            Star2 stellar state.
        binary_state : str
            Binary state.
        binary_event : str
            Binary event.
        step_name : str
            Name of posydon step.
        invert_S1S2 : bool
            If `True` isolated also sort S1_state=S2_state and S2_state=S1_state
            systems.
        chunksize : int
            Read the POSYDON binary population in chuncks to prevent OFM error.

        """
        # TODO: A check should be applied in the gien S1_state, S2_state, etc.
        # are valid.
        
        # if the user provided a single string instead of a list of strings
        if type(path_to_data) is str and ('.h5' in path_to_data):
            path_to_data = [path_to_data]
        
        # The ini file and path_to_data metallities should match.
        if len(path_to_data) != len(self.solar_metallicities):
            raise ValueError('The number of metallicities and data files do not match!')        
        
        # catch the case where the user did not provide a path to data 
        if (isinstance(path_to_data, list)):
            for path in path_to_data:
                if os.path.splitext(path)[-1] != '.h5':
                    raise ValueError('You did not provide a valid path_to_data!')
        else:
            raise ValueError('You did not provide a valid path_to_data!')
        
        #df_sel = pd.DataFrame()
        #df_sel_oneline = pd.DataFrame()
        parsed_population_file = pd.HDFStore(self.parsed_pop_file,
                                             mode='a',
                                             complevel=1,
                                             complib='lzo')
        
        count = 0
        tmp = 0
        index_shift = 0
        
        self.parse_kwargs = {'S1_state': S1_state,
                            'S2_state': S2_state,
                            'binary_state': binary_state,
                            'binary_event': binary_event,
                            'step_name': step_name,
                            'invert_S1S2': invert_S1S2}

        self.path_to_data = path_to_data
        self.underlying_mass_per_met = np.zeros(len(self.metallicities))
        self.simulated_mass_per_met = np.zeros(len(self.metallicities))
        self.selected_systems = np.zeros(len(self.metallicities))
        self.total_systems = np.zeros(len(self.metallicities))
        
        if self.verbose:
            print('Binary count with (S1_state, S2_state, binary_state, binary_event, step_name) equal')
            print(f'to ({S1_state}, {S2_state}, {binary_state}, {binary_event}, {step_name})')
            if invert_S1S2:
                print(f'and ({S2_state}, {S1_state}, {binary_state}, {binary_event}, {step_name})')
                
        for k, file in enumerate(path_to_data):
            indices_sel_met = []
            last_binary_df = None
            total = 0
            # read metallicity from path
            met = float(file.split('/')[-1].split('_Zsun')[0])*Zsun
            met_index = np.where(np.isclose(self.metallicities, met))[0][0]
            simulated_mass_for_met = 0.
            
            # get the history columns + min_itemsize
            history_cols = pd.read_hdf(file, key='history', start=0, stop=0).columns
            history_min_itemsize = {key: val for key, val in
                                    HISTORY_MIN_ITEMSIZE.items()
                                    if key in history_cols}
            
            oneline_cols = pd.read_hdf(file, key='oneline', start=0, end=0).columns
            oneline_min_itemsize = {key: val for key, val in 
                                    ONELINE_MIN_ITEMSIZE.items()
                                    if key in oneline_cols}
            
            for i, df in enumerate(pd.read_hdf(file,  key='history', chunksize=self.chunksize)):
                
                df = pd.concat([last_binary_df, df])
                    
                last_binary_df = df.loc[[df.index[-1]]]
                df.drop(df.index[-1], inplace=True)
                total += len(df.index.drop_duplicates())
                max_index = df.index.max()
                
                logic = self.apply_logic(df, 
                                         S1_state     = S1_state,
                                         S2_state     = S2_state,
                                         binary_state = binary_state,
                                         binary_event = binary_event,
                                         step_name    = step_name,
                                         invert_S1S2  = invert_S1S2)
                
                # select systems
                # remove duplicate indicies, e.g. if selecting 'contact' state it appears twice
                # if no specific event is selected (the second time is from the copied END event)
                sel = df.loc[logic].index.drop_duplicates()
                
                # count systems
                count += len(np.unique(sel))

                # read the simulated ZAMS mass
                sel_ZAMS = df['event'] == 'ZAMS'
                mass = sum(df['S1_mass'][sel_ZAMS]+df['S2_mass'][sel_ZAMS])
                simulated_mass_for_met += mass

                # sort systems
                if any(sel):
                    df_tmp = pd.DataFrame()
                    df_tmp = df.loc[sel]
                    # store metallicity
                    df_tmp['metallicity'] = met
                    indices_sel_met.extend(sel.values)

                    # shift index
                    df_tmp.index += index_shift
                    parsed_population_file.append('history',
                                                  df_tmp,
                                                  data_columns=True,
                                                  min_itemsize=history_min_itemsize)
                    
                    del df_tmp

            # check last binary if it should be included
            if last_binary_df is not None:
                max_index = np.max([max_index, last_binary_df.index.max()])
                total += 1
                logic = self.apply_logic(last_binary_df, 
                                    S1_state     = S1_state,
                                    S2_state     = S2_state,
                                    binary_state = binary_state,
                                    binary_event = binary_event,
                                    step_name    = step_name,
                                    invert_S1S2  = invert_S1S2)
                
                # The last binary is selected
                if any(logic) == True:
                    df_tmp = last_binary_df.loc[logic]
                    df_tmp['metallicity'] = met
                    indices_sel_met.extend(df_tmp.index.drop_duplicates().values)
                    df_tmp.index += index_shift
                    parsed_population_file.append('history',
                                                  df_tmp,
                                                  data_columns=True,
                                                  min_itemsize=history_min_itemsize)
                    del df_tmp
            
            # calculate population masses and parameters
            self.simulated_mass_per_met[met_index] = simulated_mass_for_met
            self.underlying_mass_per_met[met_index] = initial_total_underlying_mass(df=simulated_mass_for_met, **self.ini_params)[0] # This used to be init_kw
            self.selected_systems[met_index] = count-tmp
            self.total_systems[met_index] = total
            
            # get unique indicies
            #sel_met = df_sel_met.index.drop_duplicates()
            
            # Now select the oneline data
            # 1. Get indices from history
            unique_indices = np.unique(indices_sel_met).tolist()
            
            # 2. Select indices from oneline
            for i, df in enumerate(pd.read_hdf(file, 
                                               key='oneline',
                                               chunksize=self.chunksize,
                                               where="index in unique_indices")):
                df.index += index_shift
                df['metallicity'] = met
                parsed_population_file.append('oneline',
                                              df,
                                              data_columns=True,
                                              min_itemsize=oneline_min_itemsize)
            
            if self.verbose:
                print(f'in {file} are {count-tmp}')
            tmp = count
            # either the max of the population in the file or the total number of systems.
            index_shift += np.max([max_index, total])

        if self.verbose:
            print('Total binaries found are', count)
        
        # close the file first since the other functions open the file again
        parsed_population_file.close()    
        # store the information into the ParsedPopulation file
        self.io.save_per_metallicity_info(self)
        self.io.save_parse_kwargs(self)
        self.io.save_path_to_data(self)
        
    def head(self, n=10, type='history'):
        """Return the first n rows of the parsed population."""
        return pd.read_hdf(self.parsed_pop_file, key=type, start=0, stop=n)
    
    def iterator(self, key='history'):
        """Return an iterator over the parsed population."""
        return pd.read_hdf(self.parsed_pop_file, key=key, iterator=True)
    
    def read_binary(self, index, type='history', columns=None):
        """Return a specific binary and specified columns
        
        Parameters
        ----------
        index : int
            Index of the binary to read.
        type : str
            Type of data to read, either 'history' or 'oneline'.
        columns : list of str
            List of columns to read.
        
        Returns
        -------
        pd.DataFrame
            Dataframe containing the binary.
        """
        # No checks are done if the columns of indices are present in the file.
        return pd.read_hdf(self.parsed_pop_file, key=type, columns=columns, where=f'index=={index}')
        
    def apply_logic(self, df, S1_state=None, S2_state=None, binary_state=None,
                    binary_event=None, step_name=None, invert_S1S2=False,
                    warn=True):
        """Select binaries in a dataframe given some properties.

        Parameters
        ----------
        df : pd.DataFrame
            POSYDON binary population synthesis dataframe.
        S1_state : str
            Star1 stellar state.
        S2_state : str
            Star2 stellar state.
        binary_state : str
            Binary state.
        binary_event : str
            Binary event.
        step_name : str
            Name of posydon step.
        invert_S1S2 : bool
            If `True` isolated also sort S1_state=S2_state and S2_state=S1_state
            systems.

        Returns
        -------
        pd.DataFrame of bools
            List of binaries to select given the search parameters.

        """
        if not invert_S1S2 and S1_state != S2_state and warn:
            warnings.warn('Note that invert_S1S2=False, hence you are not parsing '
                          f'the dataset for {S1_state}-{S2_state} binaries and '
                          f'and not for for {S1_state}-{S2_state}. If this is '
                          'done on purpose, ignore this message!')

        sel_all = df['S1_state'].astype(bool)

        if S1_state is not None:
            S1_logic = (df['S1_state'] == S1_state)
        else:
            S1_logic = sel_all
        if S2_state is not None:
            S2_logic = (df['S2_state'] == S2_state)
        else:
            S2_logic = sel_all
        if binary_state is not None:
            binary_state_logic = (df['state'] == binary_state)
        else:
            binary_state_logic = sel_all
        if binary_event is not None:
            binary_event_logic = (df['event'] == binary_event)
        else:
            binary_event_logic = sel_all
        if step_name is not None:
            step_name_logic = (df['step_names'] == step_name)
        else:
            step_name_logic = sel_all

        if invert_S1S2:
            S1_logic_inverted = (df['S1_state'] == S2_state)
            S2_logic_inverted = (df['S2_state'] == S1_state)
            # find systems
            logic = ((S1_logic & S2_logic & binary_state_logic &
                        binary_event_logic & step_name_logic) |
                     (S1_logic_inverted & S2_logic_inverted &
                        binary_state_logic  & binary_event_logic &
                        step_name_logic))
        else:
            # find systems
            logic = (S1_logic & S2_logic & binary_state_logic &
                        binary_event_logic & step_name_logic)

        return logic
 
    @property
    def indices(self):
        """Return the indices of the parsed population."""
        if not hasattr(self, '_indices'):
            with pd.HDFStore(self.parsed_pop_file, mode='r') as store:
                self._indices = store.select_column('oneline', 'index').to_numpy()
        return self._indices

    def get_formation_channels(self, mt_history=False):
        """Get formation channel and add to df_oneline.
        
        Parameters
        ----------
        mt_history : bool
            If `True`, split the event oRLO1/oRLO2 into oRLO1-contact/oRLO2-contact,
            oRLO1-reverse/oRLO2-reverse and oRLO1/oRLO2. This is useful to
            identify binaries undergoing contact stable mass-transfer phases and 
            reverse mass-transfer phase .
        """
        
        # 1. Get the binary indices from the oneline df
        # 2. Read the specific columns from the history df and oneline for each binary
        # 3. Store them into a new dataframe and save it to the ParsedPopulation file
        
        if self.verbose: print('Calculating formation channels...')  
        
        
        # load the HMS-HMS interp class
        HMS_HMS_event_dict = {'stable_MT'   : 'oRLO1', 
                              'no_MT'       : 'None', 
                              'unstable_MT' : 'oCE1/oDoubleCE1'}
        
        unique_binary_indices = self.indices
        
        # check if formation channels already exist
        with pd.HDFStore(self.parsed_pop_file, mode='a') as store:
            if '/formation_channels' in store.keys():
                print('Formation channels already exist in the parsed population file!')
                print('Channels will be overwriten')
                del store['formation_channels']
        
        
        def get_events(group):
             # for now, only append information for RLO1; unstable_MT information already exists
            if 'oRLO1' in group['interp_class_HMS_HMS'].tolist():
                combined_events = group['event'].iloc[0] + '_' + group['interp_class_HMS_HMS'].iloc[0]
                tmp = [combined_events]
                tmp.extend(group['event'].iloc[1:])
                combined_events = '_'.join(tmp)
            else:
                combined_events = '_'.join(group['event'])
            return pd.Series({'channel_debug': combined_events})
        
        def mt_history(row):
            if pd.notna(row['mt_history_HMS_HMS']) and row['mt_history_HMS_HMS'] == 'Stable contact phase':
                return row['channel'].replace('oRLO1','oRLO1-contact')
            elif pd.notna(row['mt_history_HMS_HMS']) and row['mt_history_HMS_HMS'] == 'Stable reverse mass-transfer phase':
                return row['channel'].replace('oRLO1', 'oRLO1-reverse')
            else:
                return row['channel']
        
        if self.verbose:
            print('Getting history counts...')
        
        # get the length of each binary's history
        # oneline can be accessed using the counter (i) + chunksize
        # but each binary has a different history length
        with pd.HDFStore(self.parsed_pop_file, mode='r') as store:
            history_events = store.select_column('history', 'index')
        history_lengths = history_events.groupby(history_events).count()
        
        if self.verbose:
            print("Done")

        unique_binary_indices = self.indices
        
        if self.verbose: print('looping over the binaries now')
        
        previous = 0
        for i in tqdm(range(0,len(unique_binary_indices), self.chunksize), disable=not self.verbose):
            selection = unique_binary_indices[i:i+self.chunksize]
            
            # create the dataframe for the chunk
            df = pd.DataFrame(index=selection, columns=['channel_debug', 'channel'])
            end = previous + history_lengths[i:i+self.chunksize].sum()
            
            # get the history of chunk events and transform the interp_class_HMS_HMS
            with pd.HDFStore(self.parsed_pop_file, mode='r') as store:
                interp_class_HMS_HMS = store.select_column('oneline',
                                                            'interp_class_HMS_HMS',
                                                            start=i,
                                                            stop=i+self.chunksize)
                                
                interp_class_HMS_HMS.index = selection
                
                # mask the nan values
                mask = ~interp_class_HMS_HMS.isna()
                interp_class_HMS_HMS[mask] = interp_class_HMS_HMS[mask].apply(lambda x : HMS_HMS_event_dict[x])
                del mask
                
                hist = store.select_column('history',
                                            'event',
                                            start=previous,
                                            stop=end)
            
            hist.index = history_events.loc[previous:end-1]
            previous = end
            
            # combine based on the index, this allows for an easier apply later
            merged = pd.merge(hist.dropna(), interp_class_HMS_HMS, left_index=True, right_index=True)
            del hist, interp_class_HMS_HMS
            
            merged.index.name='binary_index'
            df['channel_debug'] = merged.groupby('binary_index').apply(get_events)
            del merged
            df['channel'] = df['channel_debug'].str.replace('_redirect', '').str.replace('_CO_contact', '')
            
            if mt_history:
                
                with pd.HDFStore(self.parsed_pop_file, mode='r') as store:
                    # check if mt_history in df
                    columns = store.select('oneline', start=0, stop=0).columns
                    if 'mt_history_HMS_HMS' not in columns:
                        raise ValueError('mt_history_HMS_HMS not saved in the oneline dataframe!')
                    else:
                        tmp_df = pd.DataFrame(index=selection, columns=['channel', 'mt_history_HMS_HMS'])
                        #print(df.loc[selection,'channel'])
                        tmp_df['channel'] = df['channel']
                        x = store.select_column('oneline',
                                                'mt_history_HMS_HMS',
                                                start=i,
                                                stop=i+self.chunksize)
                        
                        x.index = selection
                        tmp_df['mt_history_HMS_HMS'] = x
                        df['channel'] = tmp_df.apply(mt_history, axis=1)
                        del tmp_df
            
            self.io.write_formation_channels(self, df)
            del df
            
    @property
    def formation_channels(self):
        '''Return the formation channel if calculated
        '''
        
        with pd.HDFStore(self.parsed_pop_file, mode='r') as store:
            if '/formation_channels' in store.keys():
                return store['formation_channels']
            else:
                print('No formation channels found in the file')                
                return None
        
    def create_DCO_population(self, output_file=None, additional_oneline_cols=None):
        '''Create a DCO population from the parsed population.
        
        A DCO Population will contain one 'time' for each 'event',
        which represents the time after starburst that the event takes place.
        
        All DCO events in the population will be contained in the new population.
        If no filtering for the BH or NS formation is done, the DCO population
        will contain both.
        
        Parameters
        -----------
        oneline_cols : list of str
            List of columns to extract from the oneline dataframe and add to
            the DCO population dataframe.
            
        Returns
        --------
        DCO_synthetic_population : SyntheticPopulation
            A synthetic population containing DCOs with the time of formation
            and merger.
        '''
        
        # The user will have done an initial parse, when inputting the parsed
        # population data. Additional filtering can be done manually.
        # Using the data from the parsed population, 
        # we find the DCO at their formation.
        if output_file == None:
            DCO_synthetic_population = SyntheticPopulation(self.parsed_pop_file,
                                                           verbose=self.verbose)
            parsed_file = pd.HDFStore(self.parsed_pop_file, mode='a')
            output_store = parsed_file
            
        else:
            DCO_synthetic_population = SyntheticPopulation(output_file,
                                                           verbose=self.verbose)
            
            # write population data to the new file!
            DCO_synthetic_population.io.save_ini_params(self)
            DCO_synthetic_population.io.save_parse_kwargs(self)
            DCO_synthetic_population.io.save_per_metallicity_info(self)
            DCO_synthetic_population.io.save_path_to_data(self)
            DCO_synthetic_population.io.save_parsed_pop_path(self.parsed_pop_file)
            
            # load data into DCO class
            DCO_synthetic_population.io.load_ini_params(DCO_synthetic_population)
            DCO_synthetic_population.io.load_parse_kwargs(DCO_synthetic_population)
            DCO_synthetic_population.io.load_per_metallicity_info(DCO_synthetic_population)
            DCO_synthetic_population.io.load_path_to_data(DCO_synthetic_population)
            DCO_synthetic_population.io.load_parsed_pop_path(DCO_synthetic_population)
            
            # Open the file for writing the synthetic population
            parsed_file = pd.HDFStore(self.parsed_pop_file, mode='r')
            output_store = pd.HDFStore(DCO_synthetic_population.pop_file, mode='a')

        
        if '/synthetic' in output_store.keys():
            print('A synthetic population already exists in this file.\
                The current population will overwrite the existing one.!')                    
            output_store.remove('synthetic')
        
        DCO_synthetic_population.population_selection = 'DCO'
        DCO_synthetic_population._save_population_selection()
        
        where_BHNS = '((S1_state == "BH")'\
                    + ' & (S2_state == "NS")'\
                    + ' & (state == "detached")'\
                    + ' & (step_names == "step_SN"))'\
                    + ' | ((S1_state == "NS")' \
                    + ' & (S2_state == "BH")' \
                    + ' & (state == "detached")' \
                    + ' & (step_names == "step_SN"))'
        where_BNS = '((S1_state == "NS")'\
                    + ' & (S2_state == "NS")'\
                    + ' & (state == "detached")'\
                    + ' & (step_names == "step_SN"))'
        where_BBH = '((S1_state == "BH")'\
                    + ' & (S2_state == "BH")'\
                    + ' & (state == "detached")'\
                    + ' & (step_names == "step_SN"))'
        
        # get the history columns + min_itemsize
        history_cols = parsed_file.select('history',start=0, stop=0).columns
        history_min_itemsize = {key: val for key, val in
                                HISTORY_MIN_ITEMSIZE.items()
                                if key in history_cols}
        
        oneline_cols = parsed_file.select(key='oneline', start=0, stop=0).columns
        oneline_min_itemsize = {key: val for key, val in 
                                ONELINE_MIN_ITEMSIZE.items()
                                if key in oneline_cols}
        
        min_itemsize = history_min_itemsize
        min_itemsize.update(oneline_min_itemsize)
        min_itemsize.update({'channel': 100, 'channel_debug': 100}) 
        
        save_cols = ['S1_spin_orbit_tilt', 'S2_spin_orbit_tilt']
        if additional_oneline_cols is not None:
            for c in additional_oneline_cols:
                if c not in save_cols:
                    save_cols.append(c)

        # remove columnns from min_itemsize not in save_cols
        for c in min_itemsize.copy():
            if c not in save_cols:
                min_itemsize.pop(c)
        
        # never add time to the save_cols from the oneline dataframe
        # add based on the history data

        for selection in [where_BHNS, where_BNS, where_BBH]:
            if self.verbose:
                print(f'Parsing {selection}')
                
            # select BBH models and store them in the DCO population
            for df_synthetic in parsed_file.select('history', 
                                                where=selection,
                                                chunksize=self.chunksize):
                selected_indices = df_synthetic.index.values.tolist()
                if len(selected_indices) == 0:
                    continue
                
                # compute the inspiral timescale from the integrated orbit
                # this estimate is better than evaluating Peters approxiamtion
                time_contact = parsed_file.select('history',
                                                where='event == END & index in selected_indices',
                                                columns=['time'])
                # NOTE/TODO: we change the units of time in the dataframe to Myr
                # this might be confusion to the user? Note that Myr are convinient
                # when inspecting the data frame.
                df_synthetic['t_delay'] = (time_contact - df_synthetic[['time']])*1e-6 # Myr
                df_synthetic['time'] *= 1e-6 # Myr
                
                # If columns are not present, they're skipped.
                df_oneline = parsed_file.select('oneline',
                                                where='index in selected_indices',
                                                columns=save_cols)
                
                df_mt_channels = parsed_file.select('formation_channels',
                                                    where=f'index in selected_indices',
                                                    columns=save_cols)
               

                # add the columns to the synthetic population
                df_synthetic = pd.concat([df_synthetic, df_oneline], axis=1)
                df_synthetic = pd.concat([df_synthetic, df_mt_channels], axis=1)
        
                # store the synthetic populations
                output_store.append('synthetic',
                                df_synthetic,
                                format='table',
                                data_columns=True,
                                min_itemsize=min_itemsize
                                )
            
        output_store.close()
        parsed_file.close()
        return DCO_synthetic_population


    # TODO: add output_file to this function!!
    def create_GRB_population(self, output_file=None, GRB_properties={}, additional_oneline_cols=None):
        '''Create a GRB population from the parsed population.
        
        The Synthetic Population will contain a 'time for each 'event',
        which represents the time after starburst that the event takes place.
        
        Parameters
        ----------
        GRB_properties : dict
            Dictionary containing the GRB properties to use.
            The dictionary must contain
            
            - 'GRB_efficiency',
            - 'GRB_beaming',
            - 'E_GRB_iso_min'
        
        '''
        # The user will have done an initial parse, when inputting the parsed
        # population data. Additional filtering can be done manually.
        # Using the data from the parsed population, 
        # we find the DCO at their formation.
        if output_file is None:
            GRB_synthetic_population = SyntheticPopulation(self.parsed_pop_file,
                                                           verbose=self.verbose)
            parsed_store = pd.HDFStore(self.parsed_pop_file, mode='a')
            output_store = parsed_file
        else:
            GRB_synthetic_population = SyntheticPopulation(output_file,
                                                           verbose=self.verbose)
            
            # write population data to the new file!
            GRB_synthetic_population.io.save_ini_params(self)
            GRB_synthetic_population.io.save_parse_kwargs(self)
            GRB_synthetic_population.io.save_per_metallicity_info(self)
            GRB_synthetic_population.io.save_path_to_data(self)
            GRB_synthetic_population.io.save_parsed_pop_path(self.parsed_pop_file)
            
            # load data into GRB_data
            GRB_synthetic_population.io.load_ini_params(GRB_synthetic_population)
            GRB_synthetic_population.io.load_parse_kwargs(GRB_synthetic_population)
            GRB_synthetic_population.io.load_per_metallicity_info(GRB_synthetic_population)
            GRB_synthetic_population.io.load_path_to_data(GRB_synthetic_population)
            GRB_synthetic_population.io.load_parsed_pop_path(GRB_synthetic_population)
            
            # Open the file for writing the synthetic population
            parsed_store = pd.HDFStore(self.parsed_pop_file, mode='r')
            output_store = pd.HDFStore(GRB_synthetic_population.pop_file, mode='a')

        
        if 'synthetic' in output_store.keys():
            print('A synthetic population already exists in this file.\
                The current population will be removed!')      
            del output_store['synthetic']
        
        GRB_synthetic_population.population_selection = 'GRB'
        GRB_synthetic_population._save_population_selection()
        GRB_synthetic_population.model_parameters = GRB_properties  
        
        if ('GRB_efficiency' not in GRB_properties or
            GRB_properties['GRB_efficiency'] is None):
            raise ValueError('Missing GRB_efficiency variable in the MODEL!')
        if ('GRB_beaming' not in GRB_properties or
            GRB_properties['GRB_beaming'] is None):
            raise ValueError('Missing GRB_beaming variable in the MODEL!')
        if ('E_GRB_iso_min' not in GRB_properties or
            GRB_properties['E_GRB_iso_min'] is None):
            raise ValueError('Missing GRB_beaming variable in the MODEL!')
        
        # S1 and S2 mass are autmatically added to their respective columns
        # changing columns over the SN
        columns_pre_post = ['orbital_period', 'eccentricity']
        # unchanged columns
        columns = ['metallicity']
        # LGRB parameters
        oneline_columns = ['S1_spin_orbit_tilt']
        if additional_oneline_cols is not None:
            for c in additional_oneline_cols:
                if c not in oneline_columns:
                    oneline_columns.append(c)
                    
                    
        min_itemsize = {'channel': 100,}
        min_itemsize.update({key:val for key, val in 
                             HISTORY_MIN_ITEMSIZE.items() 
                             if key in columns_pre_post})
        min_itemsize.update({key:val for key, val in
                             HISTORY_MIN_ITEMSIZE.items()
                                if key in columns})
        min_itemsize.update({key:val for key, val in
                             ONELINE_MIN_ITEMSIZE.items()
                            if key in oneline_columns})
        
        def multi_columns_read(key, parsed_store, columns, previous, end, additional_oneline_cols=None):
            '''Read multiple columns from the history dataframe
            '''
            df = pd.DataFrame()
            for c in columns:
                df[c] = parsed_store.select_column(key, c, start=previous, stop=end)
            df.index = parsed_store.select_column(key, 'index', start=previous, stop=end)
            return df                  
        
        def GRB_data_store(tmp_df, parsed_store, previous, end, S1_S2='S1'):                    
            indices = tmp_df.index
            
            # Read history of the chunk
            hist = multi_columns_read('history',
                                    parsed_store, 
                                    ['S1_state', 'S2_state', 'event', 'step_names', 'time', 'S1_mass', 'S2_mass', 'metallicity', 'orbital_period', 'eccentricity'],
                                    previous,
                                    end)
            
            GRB_df_synthetic = pd.DataFrame()
            
            if S1_S2 == 'S1':
            # get the SN event of the GRB
            # This also select the second step_SN if S1_S2 goes SN first.
                logic1 = self.apply_logic(hist,
                                        S1_state='BH',
                                        S2_state=None,
                                        step_name='step_SN',
                                        invert_S1S2=False)
            elif S1_S2 == 'S2':
                logic1 = self.apply_logic(hist,
                                        S1_state=None,
                                        S2_state='BH',
                                        step_name='step_SN',
                                        invert_S1S2=False)   
        
            # select the previous row
            mask2 = logic1.shift(-1)
            mask2.iloc[-1]  = False
            
            # Select just the event where S1_S2 was not BH before undergoing the SN
            S1_SN = hist.loc[mask2].loc[indices][f'{S1_S2}_state'] != 'BH'
            post_sn = hist.loc[logic1].loc[indices][S1_SN]
            pre_sn = hist.loc[mask2].loc[indices][S1_SN]
            # write properties to the Synthetic Population
            GRB_df_synthetic = pd.DataFrame()
            if S1_S2 == 'S1':
                columns_pre_post.append('S1_mass')
                columns.append('S2_mass')
            elif S1_S2 == 'S2':
                columns_pre_post.append('S2_mass')
                columns.append('S1_mass')
                
            for c in columns_pre_post:
                GRB_df_synthetic['preSN_'+c] = pre_sn[c].values
                GRB_df_synthetic['postSN_'+c] = post_sn[c].values
            
            GRB_df_synthetic.index = post_sn.index
            GRB_df_synthetic['time'] = post_sn['time'].values * 1e-6 # Myr
            for c in columns:
                GRB_df_synthetic[c] = pre_sn[c].values
            

            # add oneline parameters
            df_oneline = multi_columns_read('oneline',
                                            parsed_store,
                                            oneline_columns,
                                            i,
                                            i+self.chunksize)
            
            for c in oneline_columns:
                GRB_df_synthetic[c] = df_oneline.loc[indices,c].values
            
            # Add GRB parameters
            for c in tmp_df.columns:
                GRB_df_synthetic[c] = tmp_df.loc[indices,c].values

            
            
            # add formation channels
            df_channel = parsed_store.select_column('formation_channels',
                                                    'channel', 
                                                    start = i,
                                                    stop=i+self.chunksize)
            
            df_channel.index = parsed_store.select_column('formation_channels',
                                                         'index',
                                                          start = i,
                                                          stop=i+self.chunksize)
            
            GRB_df_synthetic['channel'] = df_channel.loc[indices].values
            
            return GRB_df_synthetic
        
        history_events = parsed_store.select_column('history', 'index')
        history_lengths = history_events.groupby(history_events).count()
        del history_events
        
        unique_binary_indices = self.indices
            
        if self.verbose: print('looping over the binaries now')
        
        previous = 0
        for i in tqdm(range(0,len(unique_binary_indices), self.chunksize), disable=not self.verbose):
            
            selection = unique_binary_indices[i:i+self.chunksize]
            end = previous + history_lengths[i:i+self.chunksize].sum()
                        
            # Read oneline
            m_disk_radiated = parsed_store.select_column('oneline',
                                                         'S1_m_disk_radiated',
                                                         start = i,
                                                         stop=i+self.chunksize)
            m_disk_radiated = pd.concat([m_disk_radiated, 
                                         parsed_store.select_column('oneline',
                                                               'S2_m_disk_radiated',
                                                               start = i,
                                                               stop=i+self.chunksize)],
                                   axis=1)
            
            m_disk_radiated.index = parsed_store.select_column('oneline',
                                                              'index',
                                                              start = i,
                                                              stop=i+self.chunksize)
            
            tmp_df = get_GRB_properties(m_disk_radiated,
                                        GRB_properties['GRB_efficiency'],
                                        GRB_properties['GRB_beaming'],
                                        GRB_properties['E_GRB_iso_min']
                                        )
            
            # S1 GRBs
            S1_tmp_df = tmp_df[tmp_df['GRB1'] == True]
            S1_GRB_df_synthetic = GRB_data_store(S1_tmp_df, parsed_store, previous, end, 'S1')
            # Set all S2 columns to NaN
            for c in S1_GRB_df_synthetic.columns:
                if c in ['S2_eta', 'S2_E_GRB', 'S2_f_beaming', 'S2_E_GRB_iso', 'S2_L_GRB_iso', 'GRB2']:
                    S1_GRB_df_synthetic[c] = None
            
            
            # S2 GRBs
            S2_tmp_df = tmp_df[tmp_df['GRB2'] == True]
            S2_GRB_df_synthetic = GRB_data_store(S2_tmp_df, parsed_store, previous, end, 'S2')
            # Set all S1 columns to NaN
            for c in S2_GRB_df_synthetic.columns:
                if c in ['S1_eta', 'S1_E_GRB', 'S1_f_beaming', 'S1_E_GRB_iso', 'S1_L_GRB_iso', 'GRB1']:
                    S2_GRB_df_synthetic[c] = None
            
            out = pd.concat([S1_GRB_df_synthetic, S2_GRB_df_synthetic])
            
            out['GRB1'] = out['GRB1'].astype(bool)
            out['GRB2'] = out['GRB2'].astype(bool)
            # store the synthetic populations
            output_store.append('synthetic',
                                out,
                                format='table',
                                data_columns=True,
                                min_itemsize=min_itemsize
                                )
            previous = end
            
            
        output_store.close()
        parsed_store.close()
        return GRB_synthetic_population


    def plot_popsyn_over_grid_slice(self, grid_type, met_Zsun, **kwargs):
      # option for formation channel selection
    # needs the full history
    # This should be part of parsed population since it uses only df and df_oneline.
    # should not rely on a population class.
        """Plot popsyn over grid slice."""
        
        # TODO: does this fit into memory?
        df = pd.read_hdf(self.parsed_pop_file,
                               key='history',
                               columns=['S1_mass', 'S2_mass', 'metallicity', 'event'],)
        df_oneline = pd.read_hdf(self.parsed_pop_file,
                                 key='oneline',
                                 columns=['channel'])
        
        plot_pop.plot_popsyn_over_grid_slice(temp_pop, grid_type, met_Zsun, **kwargs)

class SyntheticPopulation():
    
    def __init__(self,
                 pop_file,
                 verbose=False,
                 
                 ):
        '''This class contains a synthetic population of transient events.

        You can calculate additional properties of the population, such as
        the formation channel, merger times, GRB properties, etc.
        
        pop_file : str
            Path to the synthetic population file.
        verbose : bool
            If `True`, print additional information.
        '''
        # initialise the parsed population class variables
        ParsedPopAttrs.__init__(self)
        self.io = ParsedPopulationIO(pop_file)
        self.pop_file = pop_file
        self.verbose = verbose
        self.population_selection = None
        if os.path.isfile(self.pop_file):
            store = pd.HDFStore(self.pop_file, mode='r')
            keys = store.keys()
            store.close()
            # Try to load the population parameters
            if '/ini_parameters' in keys:
                self.io.load_ini_params(self)
            if '/mass_per_met' in keys:
                self.io.load_per_metallicity_info(self)
            if '/parse_kwargs' in keys:
                self.io.load_parse_kwargs(self)            
            if '/path_to_data' in keys:
                self.io.load_path_to_data(self)
            if '/parsed_pop_file' in keys:
                self.io.load_parsed_pop_file(self)
            if '/parsed_pop_path' in keys:
                self.io.load_parsed_pop_path(self)
            if '/efficiency' in keys:
                self._load_efficiency()
            if '/population_selection' in keys:
                self._load_population_selection()
    
    def _load_efficiency(self):
        '''Loads the transient efficiency over metallicity from the population file'''
        with pd.HDFStore(self.pop_file, mode='r') as store:
            self.efficiency = store['efficiency'].to_dict()
            for key in self.efficiency:
                self.efficiency[key] = np.fromiter(self.efficiency[key].values(), dtype=float)
            self.met_efficiency = store['efficiency'].index.to_numpy()
    
    def _save_efficiency(self):
        '''Save the transient efficiency over metallicity to the population file.
        '''
        with pd.HDFStore(self.pop_file, mode='a') as store:
            # writing them as a fixed format to avoid problems with the
            # columns names
            store.put('efficiency',
                      pd.DataFrame(index=self.met_efficiency,
                                   data=self.efficiency),
                      format='fixed')

    def _load_population_selection(self):
        '''Loads the population selection from the synthetic population file'''
        with pd.HDFStore(self.pop_file, mode='r') as store:
            if '/population_selection' in store.keys():
                self.population_selection = store['population_selection'].values[0][0]
            
    def _save_population_selection(self):
        '''Stores the population selection to the synthetic population file'''
        with pd.HDFStore(self.pop_file, mode='a') as store:
            store.put('population_selection',
                      pd.DataFrame([self.population_selection]),
                      format='fixed')

    def export_synthetic_population(self, output_file, columns=None, chunksize=500000):
        '''Export the synthetic population to a separate file.'''
        
        synthetic_population = SyntheticPopulation(output_file)
        
        # write population data to the new file!
        synthetic_population.io.save_ini_params(self)
        synthetic_population.io.save_parse_kwargs(self)
        synthetic_population.io.save_per_metallicity_info(self)
        synthetic_population.io.save_path_to_data(self)
        synthetic_population.io.save_parsed_pop_path(self.parsed_pop_path)
        synthetic_population.population_selection = self.population_selection
        synthetic_population._save_population_selection()
        # Open the file for writing the synthetic population
        with pd.HDFStore(synthetic_population.pop_file, mode='a') as store:
            if '/synthetic' in store.keys():
                print('A synthetic population already exists in this file.\
                    The current population will be removed!')
                store.remove('synthetic')
                
            for df in pd.read_hdf(self.pop_file, key='synthetic', chunksize=chunksize):
                if columns is not None:
                    df = df[columns]
                store.append('synthetic',
                             df,
                             format='table',
                             data_columns=True)
                                       
    def head(self, n=10):
        """Return the first n rows of the synthetic population."""
        return pd.read_hdf(self.pop_file, key='synthetic', start=0, stop=n)
    
    def iloc(self, i,j):
        """Return the `i` to `j` entry of the synthetic population."""
        return pd.read_hdf(self.pop_file, key='synthetic', start=i, stop=j)
    
    def get_transient(self, index, columns=None):
        """Get a specific transient in the population."""
        return pd.read_hdf(self.pop_file,
                           key='synthetic',
                           columns=columns,
                           where=f'index=={index}')
    
    @property
    def synthetic_population(self):
        '''Returns the whole synthetic populaiton as a pandas dataframe.
        
        Warning: might be too big to load in memory!
        
        Returns
        -------
        pd.DataFrame
            Dataframe containing the synthetic population.
        '''
        return pd.read_hdf(self.pop_file, key='synthetic')
    
    @property
    def indices(self):
        if not hasattr(self, '_indices'):
            with pd.HDFStore(self.pop_file, mode='r') as store:
                self._indices = store.select_column('synthetic', 'index').to_numpy()
        return self._indices
    
    def plot_DTD(self, metallicity=None, ax=None, bins=100,):
        '''Plot the delay time distribution of the events in the population.
        
        Parameters
        ----------
        
        metallicity : float
            The metallicity in absolute metals. 
            Check `SyntheticPopulation.metallicities` for options.
            If `None`, all metallicities are combined
        ax : matplotlib.axes.Axes, optional
            The axes to plot on. If None given, a new figure will be created.
        bins : int or sequence, optional
            If bins is an int, it defines the number of equal-width bins in the
            given range (10, by default). If bins is a sequence, it defines the
            bin edges, including the rightmost edge, allowing for non-uniform
            bin widths. 
        '''
        if ax is None:
            fig, ax = plt.subplots()
        
        
        # plot all solar_metallicities
        if metallicity is None:
            # Values have to be copied to avoid overwriting the original data
            with pd.HDFStore(self.pop_file, mode='r') as store:
                # only load in 1 column at the time
                time = store.select_column('synthetic', 'time').to_numpy()
                if 't_delay' in store.select('synthetic', start=0, stop=0).columns:
                    t_delay = store.select_column('synthetic', 't_delay').to_numpy()
                    time += t_delay
            time *= 1e6 # yr
            h, bin_edges = np.histogram(time, bins=bins)
        
            h = h/np.diff(bin_edges)/self.underlying_mass_per_met.sum()
            
        else:
            if not any(np.isclose(metallicity, self.metallicities)):
                raise ValueError('The metallicity is not present in the population!')

            # get the time of the events
            # select_columns cannot be used with where
            # TODO: this might need to be rewritten to not go out of memory!
            df_tmp = pd.read_hdf(self.pop_file,
                                 key='synthetic',
                                 where='metallicity == metallicity',
                                 columns= ['time', 't_delay'])
            time = df_tmp['time'].to_numpy()
            
            # Add the delay time (GW inspiral time) if present
            if 't_delay' in df_tmp:
                time += df_tmp['t_delay'].to_numpy()
            time *= 1e6 # yr
                
            h, bin_edges = np.histogram(time, bins=bins)
            
            h = h/np.diff(bin_edges)/self.underlying_mass_per_met[np.isclose(metallicity, self.metallicities)]
        
        ax.step(bin_edges[:-1], h, where='post')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('Time [yr]')
        ax.set_ylabel('Number of events/Msun/yr')
     
    def get_efficiency_over_metallicity(self):
        """Compute the efficiency of events per Msun for each solar_metallicities."""
        
        if hasattr(self, 'efficiency'):
            print('Efficiencies already computed! Overwriting them!')
            
        efficiencies = []
        self.met_efficiency = sorted(self.metallicities, reverse=True)
        for met in self.met_efficiency:
            count = len(pd.read_hdf(self.pop_file, key='synthetic', 
                              where='metallicity == met',
                              columns=['index']))
            # just sums the number of events
            #count = np.sum(sel)
            underlying_stellar_mass = self.underlying_mass_per_met[np.isclose(met, self.metallicities)][0]
            eff = count/underlying_stellar_mass
            efficiencies.append(eff)
            print(f'Efficiency at Z={met:1.2E}: {eff:1.2E} Msun^-1')
        self.met_efficiency = np.array(self.met_efficiency)
        self.efficiency = {'total' : np.array(efficiencies)}
        # if the channel column is present compute the merger efficiency per channel
        if 'channel' in pd.read_hdf(self.pop_file,
                                              key='synthetic',
                                              start=0,
                                              stop=0).columns:   
        #if "channel" in self.df_synthetic:
            channels = np.unique(pd.read_hdf(self.pop_file,
                                             key='synthetic',
                                             columns=['channel']).values)
            for ch in channels:
                efficiencies = []
                for met in self.met_efficiency:
                    count = len(pd.read_hdf(self.pop_file,
                                          key='synthetic',
                                          where='metallicity == met & channel == ch',))

                    #count = self.df_synthetic[sel].shape[0]
                    if count > 0:
                        underlying_stellar_mass =self.underlying_mass_per_met[np.isclose(met, self.metallicities)][0]
                        eff = count/underlying_stellar_mass
                    else:
                        eff = np.nan
                    efficiencies.append(eff)
                    # print(f'Z={met:1.2E} {ch}: {eff:1.2E} Msun^-1')
                self.efficiency[ch] = np.array(efficiencies)
        
        # save the efficiency
        self._save_efficiency()
         
    def plot_met_efficiency(self, **kwargs):
        """Plot merger rate efficiency.
        
        Parameters
        ----------
        channel : bool
            plot the subchannels
        """
        if self.met_efficiency is None or self.efficiency is None:
            raise ValueError('First you need to compute the merger efficinty!')
        plot_pop.plot_merger_efficiency(self.met_efficiency, self.efficiency, **kwargs)
    
    
    def plot_hist_properties(self, var, intrinsic=False, observable=False, pop=None, **kwargs):
        """Plot histogram of intrinsic/observable properites.

        Parameters
        ----------
        var : str
            Property to plot stored in intrinsic/observable dataframe.
        intrinsic : bool
            `True` if you want to deplay the intrisc population.
        observable : bool
            `True` if you want to deplay the observable population.
        **kwargs : dict
            ploting arguments

        """
        if pop == 'DCO':
            if self.df_dco_intrinsic is None and self.df_dco_observable is None:
                raise ValueError('First you need to compute the merger rate density!')
            if intrinsic:
                df_intrinsic = self.df_dco_intrinsic
            else:
                df_intrinsic = None
            if observable:
                df_observable = self.df_dco_observable
            else:
                df_observable = None
        elif pop == 'GRB':
            if self.df_grb_intrinsic is None and self.df_grb_observable is None:
                raise ValueError('First you need to compute the merger rate density!')
            if intrinsic:
                df_intrinsic = self.df_grb_intrinsic
            else:
                df_intrinsic = None
            if observable:
                df_observable = self.df_grb_observable
            else:
                df_observable = None       
        else:
            raise ValueError('Population not recognized!')
        plot_pop.plot_hist_properties(var, df_intrinsic=df_intrinsic, df_observable=df_observable, pop=pop, **kwargs)

    def plot_rate_density(self, DCO=False, GRB=False, **kwargs):
        """Plot DCO and GRB rate densities."""
        if not DCO and not GRB:
            raise ValueError('You need to choose at least one population to plot!')
        if DCO:
            if self.dco_z_rate_density is None or self.dco_rate_density is None:
                raise ValueError('First you need to compute the merger rate density!')
            else:
                z_dco = self.dco_z_rate_density
                rate_dco =  self.dco_rate_density
        else:
            z_dco = None
            rate_dco = None
        if GRB:
            if self.grb_z_rate_density is None or self.grb_rate_density is None:
                raise ValueError('First you need to compute the GRB rate density!')
            else:
                z_grb = self.grb_z_rate_density
                rate_grb =  self.grb_rate_density
        else:
            z_grb = None
            rate_grb = None
        plot_pop.plot_rate_density(z_dco, rate_dco, z_grb, rate_grb, **kwargs)

    def plot_popsyn_over_grid_slice(self, grid_type, met_Zsun, **kwargs):
        """Plot popsyn over grid slice."""
        plot_pop.plot_popsyn_over_grid_slice(self, grid_type, met_Zsun, **kwargs)
    
    def compute_cosmological_weights(self, sensitivity, flag_pdet, working_dir, load_data, pop='DCO'):
        """Compute the GRB/DCO merger rate weights.

        Parameters
        ----------
        sensitivity : str
            GW detector sensitivity and network configuration you want to use, see arXiv:1304.0670v3
            detector sensitivities are taken from: https://dcc.ligo.org/LIGO-T2000012-v2/public
            available sensitivity keys (for Hanford, Livingston, Virgo network):
                'O3actual_H1L1V1' : aligo_O3actual_H1.txt, aligo_O3actual_L1.txt, avirgo_O3actual.txt
                'O4low_H1L1V1' : aligo_O4low.txt, aligo_O4low.txt, avirgo_O4high_NEW.txt
                'O4high_H1L1V1' : aligo_O4high.txt, aligo_O4high.txt, avirgo_O4high_NEW.txt
                'design_H1L1V1' : AplusDesign.txt, AplusDesign.txt, avirgo_O5high_NEW.txt
                'infinite': intrinsic merging DCO population, i.e. p_det = 1
        flag_pdet : bool
            `True` if you use sensitivity != 'infinite'.
        working_dir : str
            Working directory where the weights will be saved.
        load_data : bool
            `True` if you want to load the weights computed by this function
            in your working directory.
        **kwargs : dict
            Kwargs containing the model parameters of your rate calculation.
            See posydon/popsyn/rate_calculation.py

        Returns
        -------
        array doubles
            Return the cosmological weights, z_formation, z_merger and binary
            index k associated to each weighted binary.

        """

        # TODO: make the class inputs kwargs
        self.rates = Rates(self.df_synthetic, **GRB_properties)

        # compute DCO merger rate density
        if pop == 'DCO':
            if not load_data:
                self.rates.compute_merger_rate_weights(sensitivity=sensitivity, flag_pdet=flag_pdet, path_to_dir=working_dir)
            index, z_formation, z_merger, w_ijk = self.rates.load_merger_rate_weights(sensitivity, path_to_dir=working_dir)

            return index, z_formation, z_merger, w_ijk
        elif pop == 'GRB':
            if not load_data:
                self.rates.compute_GRB_rate_weights(sensitivity=sensitivity, path_to_dir=working_dir)
            index_1, z_formation_1, z_grb_1, w_ijk_1, \
                index_2, z_formation_2, z_grb_2, w_ijk_2 = self.rates.load_grb_rate_weights(sensitivity, path_to_dir=working_dir)

            return index_1, z_formation_1, z_grb_1, w_ijk_1, index_2, z_formation_2, z_grb_2, w_ijk_2
        else:
            raise ValueError('Population not recognized!')  
    
    def resample_synthetic_population(self, index, z_formation, z_event, w_ijk, export_cols=None, 
                                      pop='DCO', reset_grb_properties=None):
        """Resample synthetc population to obtain intrinsic/observable population.

        Parameters
        ----------
        index : array int
            Index k of each binary corresponding to the synthetc dataframe
            proeprties.
        z_formation : array float
            Redshift of formation of each binary.
        z_merger : array float
            Redshift of merger of each binary.
        w_ijk : array float
            Cosmological weights computed with Eq. B.8 of Bavera et at. (2020).
        export_cols : list str
            List of additional columns to save in the intrinsic/observable
            population.

        Returns
        -------
        pd.DataFrame
            Resampled synthetc population to intrinsic or detecatable
            population.

        """
        # compute GRB properties boolean
        compute_GRB_properties = (GRB_properties is not None and 
                                  "compute_GRB_properties" in GRB_properties and 
                                  GRB_properties["compute_GRB_properties"])

        # drop all zero weights to save memory
        sel = w_ijk > 0.
        index = index[sel]
        z_formation = z_formation[sel]
        z_event = z_event[sel]
        w_ijk = w_ijk[sel]

        # export results, we do a selections of columns else the dataframe
        # risks to go out of memory
        df = pd.DataFrame()
        df['weight'] = w_ijk
        df['z_formation'] = z_formation
        if pop == 'DCO':
            df['z_merger'] = z_event
        elif pop == 'GRB':
            df['z_grb'] = z_event
        else:
            raise ValueError('Population not recognized!')
        save_cols = ['metallicity','time','t_delay','S1_state','S2_state',
                     'S1_mass','S2_mass','S1_spin','S2_spin',
                     'orbital_period','eccentricity', 'q', 'm_tot',
                     'm_chirp', 'chi_eff']
        if compute_GRB_properties:
            save_cols += GRB_PROPERTIES
        if "channel" in self.rates.df:
            save_cols.append('channel')
        if export_cols is not None:
            for c in export_cols:
                if c not in save_cols:
                    save_cols.append(c)
        for c in save_cols:
            df[c] = self.rates.get_data(c, index)
        
        # the same binary system can emit two GRBs, we remove the
        # GRB properties of the other GRB to prevent mistakes
        if reset_grb_properties is not None:
            if reset_grb_properties == 'GRB1':
                for key in GRB_PROPERTIES:
                    if "S1" in key:
                        df[key] = np.nan
            elif reset_grb_properties == 'GRB2':
                for key in GRB_PROPERTIES:
                    if "S2" in key:
                        df[key] = np.nan
            else:
                raise ValueError(f'reset_grb_properties=={reset_grb_properties} key not recognized!')

        return df
    
    def get_dco_merger_rate_density(self, export_cols=None,  working_dir='./', load_data=False):
        """Compute the merger rate density as a function of redshift.

        Parameters
        ----------
        export_cols : list str
            List of additional columns to save in the intrinsic/observable
            population.
        working_dir : str
            Working directory where the weights will be saved.
        load_data : bool
            `True` if you want to load the weights computed by this function
            in your working directory.

        """
        if self.df_synthetic is None:
            raise ValueError('You first need to isolated the DCO synthetic population!')

        # compute cosmological weights (detection rate weights with infinite sensitivity)
        sensitivity='infinite'
        flag_pdet = False
        index, z_formation, z_merger, w_ijk = self.compute_cosmological_weights(sensitivity, flag_pdet, working_dir=working_dir, load_data=load_data, pop='DCO')
        # compute rate density weights
        self.dco_z_rate_density = self.rates.get_centers_redshift_bins()
        total_rate = self.rates.compute_rate_density(w_ijk, z_merger, observable='DCO', sensitivity=sensitivity)
        self.dco_rate_density = {'total' : total_rate}
        if "channel" in self.df_synthetic:
            channels = np.unique(self.df_synthetic['channel'])
            for ch in tqdm(channels):
                sel = (self.rates.get_data('channel', index) == ch)
                rate = self.rates.compute_rate_density(w_ijk[sel], z_merger[sel], observable='DCO', sensitivity=sensitivity)
                self.dco_rate_density[ch] = rate
            
        print(f'DCO merger rate density in the local Universe (z={self.dco_z_rate_density[0]:1.2f}): {round(total_rate[0],2)} Gpc^-3 yr^-1')

        # export the intrinsic DCO population
        self.df_dco_intrinsic = self.resample_synthetic_population(index, z_formation, z_merger, w_ijk, export_cols=export_cols, pop='DCO')

    def get_dco_detection_rate(self, sensitivity='design_H1L1V1', export_cols=None,  working_dir='./', load_data=False):
        """Compute the detection rate per yr.

        Parameters
        ----------
        sensitivity : str
            GW detector sensitivity and network configuration you want to use, see arXiv:1304.0670v3
            detector sensitivities are taken from: https://dcc.ligo.org/LIGO-T2000012-v2/public
            available sensitivity keys (for Hanford, Livingston, Virgo network):
                'O3actual_H1L1V1' : aligo_O3actual_H1.txt, aligo_O3actual_L1.txt, avirgo_O3actual.txt
                'O4low_H1L1V1' : aligo_O4low.txt, aligo_O4low.txt, avirgo_O4high_NEW.txt
                'O4high_H1L1V1' : aligo_O4high.txt, aligo_O4high.txt, avirgo_O4high_NEW.txt
                'design_H1L1V1' : AplusDesign.txt, AplusDesign.txt, avirgo_O5high_NEW.txt
        export_cols : list str
            List of additional columns to save in the intrinsic/observable
            population.
        working_dir : str
            Working directory where the weights will be saved.
        load_data : bool
            `True` if you want to load the weights computed by this function
            in your working directory.

        """
        if self.df_synthetic is None:
            raise ValueError('You first need to isolated the DCO synthetic population!')

        # compute detection rate weights
        flag_pdet = True
        index, z_formation, z_merger, w_ijk = self.compute_cosmological_weights(sensitivity, flag_pdet, working_dir=working_dir, load_data=load_data)
        print(f'DCO detection rate at {sensitivity} sensitivity: {sum(w_ijk):1.2f} yr^-1')

        # export the observable DCO population
        # TODO: store p_det
        self.df_dco_observable = self.resample_synthetic_population(index, z_formation, z_merger, w_ijk, export_cols=export_cols)


    def get_grb_rate_density(self, export_cols=None,  working_dir='./', load_data=False):
        """Compute the GRB density as a function of redshift.

        Parameters
        ----------
        export_cols : list str
            List of additional columns to save in the intrinsic/observable
            population.
        working_dir : str
            Working directory where the weights will be saved.
        load_data : bool
            `True` if you want to load the weights computed by this function
            in your working directory.

        """
        if self.df_synthetic is None:
            raise ValueError('You first need to isolated the DCO synthetic population!')

        # compute cosmological weights (detection rate weights with infinite sensitivity)
        sensitivity='infinite'
        flag_pdet = False
        index_1, z_formation_1, z_grb_1, w_ijk_1, \
            index_2, z_formation_2, z_grb_2, w_ijk_2 = self.compute_cosmological_weights(sensitivity, flag_pdet, working_dir=working_dir, load_data=load_data, pop='GRB')

        # compute beamed rate density weights
        # note we compute the beamed rate density for GRBs as this is what is often reported in the literature
        # as the beaming factor is not known a priori
        self.grb_z_rate_density = self.rates.get_centers_redshift_bins()
        total_GRB1 = self.rates.compute_rate_density(w_ijk_1, z_grb_1, observable='GRB1', sensitivity='beamed', index=index_1)
        total_GRB2 = self.rates.compute_rate_density(w_ijk_2, z_grb_2, observable='GRB2', sensitivity='beamed', index=index_2)
        total_rate = total_GRB1 + total_GRB2
        self.grb_rate_density = {'total' : total_rate, 'total_GRB1' : total_GRB1, 'total_GRB2': total_GRB2}
        if "channel" in self.df_synthetic:
            channels = np.unique(self.df_synthetic['channel'])
            for ch in tqdm(channels):
                sel1 = (self.rates.get_data('channel', index_1) == ch)
                if any(sel1):
                    self.grb_rate_density[ch+'_GRB1'] = self.rates.compute_rate_density(w_ijk_1[sel1], z_grb_1[sel1], observable='GRB1', sensitivity='beamed', index=index_1[sel1])
                else:
                    self.grb_rate_density[ch+'_GRB1'] = np.zeros(len(self.grb_z_rate_density))
                sel2 = (self.rates.get_data('channel', index_2) == ch)
                if any(sel2):
                    self.grb_rate_density[ch+'_GRB2'] = self.rates.compute_rate_density(w_ijk_2[sel2], z_grb_2[sel2], observable='GRB2', sensitivity='beamed', index=index_2[sel2])
                else:
                    self.grb_rate_density[ch+'_GRB2'] = np.zeros(len(self.grb_z_rate_density))
                self.grb_rate_density[ch] = self.grb_rate_density[ch+'_GRB1'] + self.grb_rate_density[ch+'_GRB2']
            
        print(f'GRB (beamed) rate density in the local Universe (z={self.grb_z_rate_density[0]:1.2f}): {round(total_rate[0],2)} Gpc^-3 yr^-1')

        # export the intrinsic grb intrisic population
        # TODO: instead of concatenating two dataframe and duplicating the information of the same
        # binary system we should combine the two datraframes where we have two columns for GRB1 and GRB2
        # such dataframe will have z_grb_1, z_grb_2, weight_1, weight_2
        # when this is addressed, remove reset_grb_properties feature
        self.df_grb_intrinsic = pd.DataFrame()
        df_grb_1 = self.resample_synthetic_population(index_1, z_formation_1, z_grb_1, w_ijk_1, export_cols=export_cols, pop='GRB', reset_grb_properties='GRB2')
        df_grb_2 = self.resample_synthetic_population(index_2, z_formation_2, z_grb_2, w_ijk_2, export_cols=export_cols, pop='GRB', reset_grb_properties='GRB1')
        self.df_grb_intrinsic = pd.concat([df_grb_1, df_grb_2], ignore_index=True, sort=False)
        # the observable population accounts for beaming
        self.df_grb_observable = self.df_grb_intrinsic.copy()
        for i in [1,2]:
            sel = self.df_grb_observable[f'S{i}_f_beaming'] > 0 
            self.df_grb_observable.loc[sel,'weight'] *= self.df_grb_observable.loc[sel,f'S{i}_f_beaming']
    
    def get_dco_detection_rate(self, sensitivity='design_H1L1V1', export_cols=None,  working_dir='./', load_data=False):
        """Compute the detection rate per yr.

        Parameters
        ----------
        sensitivity : str
            GW detector sensitivity and network configuration you want to use, see arXiv:1304.0670v3
            detector sensitivities are taken from: https://dcc.ligo.org/LIGO-T2000012-v2/public
            available sensitivity keys (for Hanford, Livingston, Virgo network):
                'O3actual_H1L1V1' : aligo_O3actual_H1.txt, aligo_O3actual_L1.txt, avirgo_O3actual.txt
                'O4low_H1L1V1' : aligo_O4low.txt, aligo_O4low.txt, avirgo_O4high_NEW.txt
                'O4high_H1L1V1' : aligo_O4high.txt, aligo_O4high.txt, avirgo_O4high_NEW.txt
                'design_H1L1V1' : AplusDesign.txt, AplusDesign.txt, avirgo_O5high_NEW.txt
        export_cols : list str
            List of additional columns to save in the intrinsic/observable
            population.
        working_dir : str
            Working directory where the weights will be saved.
        load_data : bool
            `True` if you want to load the weights computed by this function
            in your working directory.

        """
        if self.df_synthetic is None:
            raise ValueError('You first need to isolated the DCO synthetic population!')

        # compute detection rate weights
        flag_pdet = True
        index, z_formation, z_merger, w_ijk = self.compute_cosmological_weights(sensitivity, flag_pdet, working_dir=working_dir, load_data=load_data)
        print(f'DCO detection rate at {sensitivity} sensitivity: {sum(w_ijk):1.2f} yr^-1')

        # export the observable DCO population
        # TODO: store p_det
        self.df_dco_observable = self.resample_synthetic_population(index, z_formation, z_merger, w_ijk, export_cols=export_cols)
   
    def save(self, path):
        # overwrite the save function to also store the MODEL parameters
        return
        
    def load(self, path):
        # overwrite the save function.
        return 
                

##-------------------------------------------------------------------------####    

class SyntheticPopulation2:

    def __init__(self, path_to_ini, verbose=False, MODEL={}):
        """
        Parameters
        ----------

        path : str
            Path to the inifile to parse. You can supply a list in the
            metallicity parameter to evolve more than one population.
        """
        self.synthetic_pop_params = None
        self.solar_metallicities = None
        self.binary_populations = None

        self.verbose = verbose
        GRB_properties = MODEL
        # DCOs
        self.df = None
        self.df_oneline = None
        self.df_synthetic = None
        self.df_dco_intrinsic = None
        self.df_dco_observable = None
        self.met_efficiency = None
        self.efficiency = None
        self.dco_z_rate_density = None
        self.dco_rate_density = None
        # GRBs
        self.df_grb_intrinsic = None
        self.df_grb_observable = None
        self.grb_z_rate_density = None
        self.grb_rate_density = None

        if '.ini' not in path_to_ini:
            raise ValueError('You did not provide a valid path_to_ini!')
        else:
            self.synthetic_pop_params = binarypop_kwargs_from_ini(path_to_ini)
            self.solar_metallicities = self.synthetic_pop_params['metallicity']
            if not isinstance( self.solar_metallicities, list):
                self.solar_metallicities = [self.solar_metallicities]
            self.binary_populations = None
        
    def create_binary_populations(self):
        """Create a list of BinaryPopulation objects."""
        self.binary_populations = []
        for met in self.metallicities[::-1]:
            ini_kw = copy.deepcopy(self.synthetic_pop_params)
            ini_kw['metallicity'] = met
            ini_kw['temp_directory'] = self.create_met_prefix(met) + self.synthetic_pop_params['temp_directory']
            self.binary_populations.append(BinaryPopulation(**ini_kw))

    def get_ini_kw(self):
        return self.synthetic_pop_params.copy()

    def evolve(self):
        """Evolve population(s) at given Z(s)."""
        if self.binary_populations is None:
            self.create_binary_populations()
        while self.binary_populations:
            pop =  self.binary_populations.pop()
            
            if self.verbose:
                print(f'Z={pop.kwargs["metallicity"]:.2e} Z_sun')
                            
            process = mp.Process(target=pop.evolve)
            process.start()
            process.join()
            
            pop.close()
            del pop


    def merge_parallel_runs(self, path_to_batches):
        """
        Merge the folder or list of folders into a single file per metallicity.

        Parameters
        ----------
        path_to_batches : str or list of str
            Path to the folder(s) containing the batch folders.
        """
        
        if isinstance(path_to_batches, str):
            path_to_batches = [path_to_batches]
        # check if path_to_batches is the same length as the number of solar_metallicities
        if len(path_to_batches) != len(self.solar_metallicities):
            raise ValueError('The number of metallicity and batch directories do not match!')

        for met, path_to_batch in zip(self.solar_metallicities, path_to_batches):
            met_prefix = self.create_met_prefix(met)
            tmp_files = [os.path.join(path_to_batch, f)     \
                         for f in os.listdir(path_to_batch) \
                            if os.path.isfile(os.path.join(path_to_batch, f))]
            
            BinaryPopulation(**self.get_ini_kw()).combine_saved_files(met_prefix+ 'population.h5', tmp_files)
            print(f'Population at Z={met:.2e} Z_sun successfully merged!')
            if len(os.listdir(path_to_batch)) == 0:
                os.rmdir(path_to_batch)
            elif self.verbose:
                print(f'{path_to_batch} is not empty, it was not removed!')

    @staticmethod
    def create_met_prefix(met):
        """Append a prefix to the name of directories for batch saving."""
        return convert_metallicity_to_string(met) + '_Zsun_'

    def apply_logic(self, df, S1_state=None, S2_state=None, binary_state=None,
                    binary_event=None, step_name=None, invert_S1S2=False,
                    warn=True):
        """Select binaries in a dataframe given some properties.

        Parameters
        ----------
        df : pd.DataFrame
            POSYDON binary population synthesis dataframe.
        S1_state : str
            Star1 stellar state.
        S2_state : str
            Star2 stellar state.
        binary_state : str
            Binary state.
        binary_event : str
            Binary event.
        step_name : str
            Name of posydon step.
        invert_S1S2 : bool
            If `True` isolated also sort S1_state=S2_state and S2_state=S1_state
            systems.

        Returns
        -------
        pd.DataFrame of bools
            List of binaries to select given the search parameters.

        """
        if not invert_S1S2 and S1_state != S2_state and warn:
            warnings.warn('Note that invert_S1S2=False, hence you are not parsing '
                          f'the dataset for {S1_state}-{S2_state} binaries and '
                          f'and not for for {S1_state}-{S2_state}. If this is '
                          'done on purpose, ignore this message!')

        sel_all = df['S1_state'].astype(bool)

        if S1_state is not None:
            S1_logic = (df['S1_state'] == S1_state)
        else:
            S1_logic = sel_all
        if S2_state is not None:
            S2_logic = (df['S2_state'] == S2_state)
        else:
            S2_logic = sel_all
        if binary_state is not None:
            binary_state_logic = (df['state'] == binary_state)
        else:
            binary_state_logic = sel_all
        if binary_event is not None:
            binary_event_logic = (df['event'] == binary_event)
        else:
            binary_event_logic = sel_all
        if step_name is not None:
            step_name_logic = (df['step_names'] == step_name)
        else:
            step_name_logic = sel_all

        if invert_S1S2:
            S1_logic_inverted = (df['S1_state'] == S2_state)
            S2_logic_inverted = (df['S2_state'] == S1_state)
            # find systems
            logic = ((S1_logic & S2_logic & binary_state_logic &
                        binary_event_logic & step_name_logic) |
                     (S1_logic_inverted & S2_logic_inverted &
                        binary_state_logic  & binary_event_logic &
                        step_name_logic))
        else:
            # find systems
            logic = (S1_logic & S2_logic & binary_state_logic &
                        binary_event_logic & step_name_logic)

        return logic


    def parse(self, path_to_data, S1_state=None, S2_state=None, binary_state=None,
               binary_event=None, step_name=None, invert_S1S2=False, chunksize=500000):
        """Sort binaries of interests given some properties.
        
        It also stores the underlying stellar mass and 
        the initial simulated stellar mass for each metallicity.

        Parameters
        ----------
        S1_state : str
            Star1 stellar state.
        S2_state : str
            Star2 stellar state.
        binary_state : str
            Binary state.
        binary_event : str
            Binary event.
        step_name : str
            Name of posydon step.
        invert_S1S2 : bool
            If `True` isolated also sort S1_state=S2_state and S2_state=S1_state
            systems.
        chunksize : int
            Read the POSYDON binary population in chuncks to prevent OFM error.

        """
        
        # if the user provided a single string instead of a list of strings
        if type(path_to_data) is str and ('.h5' in path_to_data):
            path_to_data = [path_to_data]
        
        # catch the case where the user did not provide a path to data 
        if (isinstance(path_to_data, list)):
            for path in path_to_data:
                if os.path.splitext(path)[-1] != '.h5':
                    raise ValueError('You did not provide a valid path_to_data!')
        else:
            raise ValueError('You did not provide a valid path_to_data!')
        
        df_sel = pd.DataFrame()
        df_sel_oneline = pd.DataFrame()
        count = 0
        tmp = 0
        shift_index = 0
        if self.verbose:
            print('Binary count with (S1_state, S2_state, binary_state, binary_event, step_name) equal')
            print(f'to ({S1_state}, {S2_state}, {binary_state}, {binary_event}, {step_name})')
            if invert_S1S2:
                print(f'and ({S2_state}, {S1_state}, {binary_state}, {binary_event}, {step_name})')
        for k, file in enumerate(path_to_data):
            df_sel_met = pd.DataFrame()
            sel_met = []
            last_binary_df = None
            
            # read metallicity from path
            met = float(file.split('/')[-1].split('_Zsun')[0])*Zsun
            simulated_mass_for_met = 0.
            
            for i, df in enumerate(pd.read_hdf(file,  key='history', chunksize=chunksize)):
                
                df = pd.concat([last_binary_df, df])
                    
                last_binary_df = df.loc[[df.index[-1]]]
                df.drop(df.index[-1], inplace=True)
                
                logic = self.apply_logic(df, 
                                         S1_state     = S1_state,
                                         S2_state     = S2_state,
                                         binary_state = binary_state,
                                         binary_event = binary_event,
                                         step_name    = step_name,
                                         invert_S1S2  = invert_S1S2)
                
                # select systems
                # remove duplicate indicies, e.g. if selecting 'contact' state it appears twice
                # if no specific event is selected (the second time is from the copied END event)
                sel = df.loc[logic].index.drop_duplicates()
                
                # count systems
                count += len(np.unique(sel))

                # read the simulated ZAMS mass
                sel_ZAMS = df['event'] == 'ZAMS'
                mass = sum(df['S1_mass'][sel_ZAMS]+df['S2_mass'][sel_ZAMS])
                simulated_mass_for_met += mass

                # sort systems
                if any(sel):
                    df_tmp = pd.DataFrame()
                    df_tmp = df.loc[sel]
                    # store metallicity
                    df_tmp['metallicity'] = met
                    # concatenate results
                    df_sel_met = pd.concat([df_sel_met, df_tmp])
                    del df_tmp

            # check last binary if it should be included
            if last_binary_df is not None:
                logic = self.apply_logic(last_binary_df, 
                                    S1_state     = S1_state,
                                    S2_state     = S2_state,
                                    binary_state = binary_state,
                                    binary_event = binary_event,
                                    step_name    = step_name,
                                    invert_S1S2  = invert_S1S2)
                
                # The last binary is selected
                if any(logic) == True:
                    df_tmp = last_binary_df.loc[logic]
                    df_sel_met = pd.concat([df_sel_met, df_tmp])

            # get unique indicies
            sel_met = df_sel_met.index.drop_duplicates()
            
            if k > 0 and df_sel.shape[0] > 0:
                shift_index += max(np.unique(df.index)) + 1 
            else:
                shift_index = 0
            
            # store simulated and underlying stellar mass
            df_sel_met['simulated_mass_for_met'] = simulated_mass_for_met
            df_sel_met['underlying_mass_for_met'] = initial_total_underlying_mass(df=simulated_mass_for_met, **self.synthetic_pop_params)[0] # This used to be init_kw

            # concatenate results with shifted indices for each metallicity
            df_sel_met.index += shift_index
            df_sel = pd.concat([df_sel, df_sel_met])
            del df_sel_met

            # load, parse and store oneline dataframe
            # this dataframe is smaller, we can load it all at once
            df_sel_met_oneline = pd.read_hdf(file,  key='oneline')
            df_sel_met_oneline = df_sel_met_oneline.loc[sel_met]
            df_sel_met_oneline['metallicity'] = met

            df_sel_met_oneline.index += shift_index
            df_sel_oneline = pd.concat([df_sel_oneline, df_sel_met_oneline])

            if self.verbose:
                print(f'in {file} are {count-tmp}')
            tmp = count
            
            # shift the index for the next metallicity
            shift_index += max(np.unique(df.index)) + 1

        if self.verbose:
            print('Total binaries found are', count)

        # save parsed population as synthetic population
        if self.df is not None:
            warnings.warn('Overwriting the df population!')
        self.df = df_sel
        self.df_oneline = df_sel_oneline

    def save_pop(self, path='./selfulation.h5'):
        """Save parsed population.

        Parameters
        ----------
        path : str
            Path to file where you want to export the dataset.

        """
        if self.df is None:
            raise ValueError('Nothing to save! The population was not parsed.')
        elif self.df_oneline is None:
            raise ValueError('Missing oneline dataframe!')
        else:
            self.df.to_hdf(path, key='history')
            self.df_oneline.to_hdf(path, key='oneline')
            if self.verbose:
                print('Population successfully saved!')

    def load_pop(self, path):
        """Load parsed population.

        Parameters
        ----------
        path : str
            Path to dataset.

        """
        if self.df is None and self.df_oneline is None :
            self.df = pd.read_hdf(path, key='history')
            self.df_oneline = pd.read_hdf(path, key='oneline')
            if self.verbose:
                print('Population successfully loaded!')
        else:
            raise ValueError('You already have a population stored in memory!')

    def get_dco_at_formation(self, S1_state, S2_state, oneline_cols=None, formation_channels=False, mt_history=False):
        """Populates `df_synthetic` with DCOs at their formation.
        
        If `formation_channels` is `True` the `channel` column is added to the
        `df_synthetic` dataframe.

        if MODEL on class initialization is not None and
        "compute_GRB_properties" in MODEL. 
        If MODEL["compute_GRB_properties"] is `True` the following columns are
        in the MODEL:

            - 'GRB_efficiency',
            - 'GRB_beaming',
            - 'E_GRB_iso_min'

        The following columns are added to the `df_synthetic` dataframe:
            - S1_m_disk_radiated
            - S2_m_disk_radiated

        Note: by default this function looks for the symmetric state
        S1_state = S2_sate and S2_state = S1_sate.

        Parameters
        ----------
        S1_state : str
            Star1 stellar state.
        S2_state : str
            Star2 stellar state.
        oneline_cols : list str
            List of columns preset in the oneline dataframe you want to export
            into the synthetic population.
        formation_channels : bool
            Compute the formation channel, a string containing the binary
            event evolution.
        mt_history : bool
            If `True`, split the event oRLO1/oRLO2 into oRLO1-contact/oRLO2-contact,
            oRLO1-reverse/oRLO2-reverse and oRLO1/oRLO2. This is useful to
            identify binaries undergoing contact stable mass-transfer phases and
            reverse mass-transfer phase .

        """
        # compute GRB properties boolean
        compute_GRB_properties = (GRB_properties is not None and 
                                  "compute_GRB_properties" in GRB_properties and 
                                  GRB_properties["compute_GRB_properties"])

        # add channel column to oneline dataframe
        if formation_channels:
            if self.verbose:
                print('Computing formation channels...')
            self.get_formation_channels(mt_history=mt_history)

        # to avoid the user making mistake automatically check the inverse of
        # the stellar states, since the df is already parsed this will not
        # take too much extra time
        logic = self.apply_logic(self.df, S1_state=S1_state,
                                          S2_state=S2_state,
                                          binary_state='detached',
                                          step_name = 'step_SN',
                                          invert_S1S2=True)
        self.df_synthetic = self.df.loc[logic].copy()

        # compute the inspiral timescale from the integrated orbit
        # this estimate is better than evaluating Peters approxiamtion
        time_contact = self.df.loc[self.df['event'] == 'END',['time']]
        # NOTE/TODO: we change the units of time in the dataframe to Myr
        # this might be confusion to the user? Note that Myr are convinient
        # when inspecting the data frame.
        self.df_synthetic['t_delay'] = (time_contact - self.df_synthetic[['time']])*1e-6 # Myr
        self.df_synthetic['time'] *= 1e-6 # Myr

        # add properties of the oneline dataframe
        if self.df_oneline is not None:
            # TODO: add kicks as well by default?
            save_cols = ['S1_spin_orbit_tilt', 'S2_spin_orbit_tilt']
            if compute_GRB_properties:
                save_cols += ['S1_m_disk_radiated', 'S2_m_disk_radiated']
            if formation_channels:
                save_cols.append('channel')
            if oneline_cols is not None:
                for c in oneline_cols:
                    if c not in save_cols:
                        save_cols.append(c)
            for c in save_cols:
                if c in self.df_oneline:
                    self.df_synthetic[c] = self.df_oneline[c]
                else:
                    warnings.warn(f'The column {c} is not present in the '
                                   'oneline dataframe.')

        # compute GRB properties
        if compute_GRB_properties:
            if ('GRB_efficiency' not in GRB_properties or
                GRB_properties['GRB_efficiency'] is None):
                raise ValueError('Missing GRB_efficiency variable in the MODEL!')
            if ('GRB_beaming' not in GRB_properties or
                GRB_properties['GRB_beaming'] is None):
                raise ValueError('Missing GRB_beaming variable in the MODEL!')
            if ('E_GRB_iso_min' not in GRB_properties or
                GRB_properties['E_GRB_iso_min'] is None):
                raise ValueError('Missing GRB_beaming variable in the MODEL!')
            self.df_synthetic = get_GRB_properties(self.df_synthetic,
                                                   GRB_properties['GRB_efficiency'],
                                                   GRB_properties['GRB_beaming'],
                                                   GRB_properties['E_GRB_iso_min']
                                                   )
            # get time_CC1 and time_CC2 note that for GRB calculations we
            # do not use the "time" column indicating the time of formation
            # of the DCO systems as there might be cases where
            # CC2 might happen before CC1 due to mass ratio reversal
            DCO = [[S1_state, None], [None, S2_state]]
            events = ['CC1', 'CC2']
            for i, event in enumerate(events):
                logic = self.apply_logic(self.df, S1_state=DCO[i][0],
                            S2_state=DCO[i][1],
                            binary_state='detached',
                            step_name = 'step_SN',
                            invert_S1S2=False,
                            warn=False)
                last_index = -1
                time = []
                for index, value in self.df.loc[logic,'time'].items():
                    if index != last_index:
                        time.append(value)
                    last_index = index
                if len(time) != self.df_synthetic.shape[0]:
                    raise ValueError('Missing values in time_{event}!')
                self.df_synthetic[f'time_{event}'] = np.array(time)*1e-6 # Myr

        # for convinience reindex the DataFrame
        n_rows = len(self.df_synthetic.index)
        self.df_synthetic = self.df_synthetic.set_index(np.linspace(0,n_rows-1,n_rows,dtype=int))

    def save_synthetic_pop(self, path='./synthetic_population.h5'):
        """Save synthetc population.

        Parameters
        ----------
        path : str
            Path to dataset.

        """
        if self.df_synthetic is None:
            raise ValueError('Nothing to save!')
        else:
            self.df_synthetic.to_hdf(path, key='history')
            if self.verbose:
                print('Synthetic population successfully saved!')

    def load_synthetic_pop(self, path):
        """Load synthetc population.

        Parameters
        ----------
        path : str
            Path to dataset.

        """
        if self.df_synthetic is None:
            self.df_synthetic = pd.read_hdf(path, key='history')
            if self.verbose:
                print('Synthetic population successfully loaded!')
        else:
            raise ValueError('You already have an synthetic population stored in memory!')




    def compute_cosmological_weights(self, sensitivity, flag_pdet, working_dir, load_data, pop='DCO'):
        """Compute the GRB/DCO merger rate weights.

        Parameters
        ----------
        sensitivity : str
            GW detector sensitivity and network configuration you want to use, see arXiv:1304.0670v3
            detector sensitivities are taken from: https://dcc.ligo.org/LIGO-T2000012-v2/public
            available sensitivity keys (for Hanford, Livingston, Virgo network):
                'O3actual_H1L1V1' : aligo_O3actual_H1.txt, aligo_O3actual_L1.txt, avirgo_O3actual.txt
                'O4low_H1L1V1' : aligo_O4low.txt, aligo_O4low.txt, avirgo_O4high_NEW.txt
                'O4high_H1L1V1' : aligo_O4high.txt, aligo_O4high.txt, avirgo_O4high_NEW.txt
                'design_H1L1V1' : AplusDesign.txt, AplusDesign.txt, avirgo_O5high_NEW.txt
                'infinite': intrinsic merging DCO population, i.e. p_det = 1
        flag_pdet : bool
            `True` if you use sensitivity != 'infinite'.
        working_dir : str
            Working directory where the weights will be saved.
        load_data : bool
            `True` if you want to load the weights computed by this function
            in your working directory.
        **kwargs : dict
            Kwargs containing the model parameters of your rate calculation.
            See posydon/popsyn/rate_calculation.py

        Returns
        -------
        array doubles
            Return the cosmological weights, z_formation, z_merger and binary
            index k associated to each weighted binary.

        """

        # TODO: make the class inputs kwargs
        self.rates = Rates(self.df_synthetic, **GRB_properties)

        # compute DCO merger rate density
        if pop == 'DCO':
            if not load_data:
                self.rates.compute_merger_rate_weights(sensitivity=sensitivity, flag_pdet=flag_pdet, path_to_dir=working_dir)
            index, z_formation, z_merger, w_ijk = self.rates.load_merger_rate_weights(sensitivity, path_to_dir=working_dir)

            return index, z_formation, z_merger, w_ijk
        elif pop == 'GRB':
            if not load_data:
                self.rates.compute_GRB_rate_weights(sensitivity=sensitivity, path_to_dir=working_dir)
            index_1, z_formation_1, z_grb_1, w_ijk_1, \
                index_2, z_formation_2, z_grb_2, w_ijk_2 = self.rates.load_grb_rate_weights(sensitivity, path_to_dir=working_dir)

            return index_1, z_formation_1, z_grb_1, w_ijk_1, index_2, z_formation_2, z_grb_2, w_ijk_2
        else:
            raise ValueError('Population not recognized!')

    def resample_synthetic_population(self, index, z_formation, z_event, w_ijk, export_cols=None,
                                      pop='DCO', reset_grb_properties=None):
        """Resample synthetc population to obtain intrinsic/observable population.

        Parameters
        ----------
        index : array int
            Index k of each binary corresponding to the synthetc dataframe
            properties.
        z_formation : array float
            Redshift of formation of each binary.
        z_merger : array float
            Redshift of merger of each binary.
        w_ijk : array float
            Cosmological weights computed with Eq. B.8 of Bavera et at. (2020).
        export_cols : list str
            List of additional columns to save in the intrinsic/observable
            population.

        Returns
        -------
        pd.DataFrame
            Resampled synthetc population to intrinsic or detecatable
            population.

        """
        # compute GRB properties boolean
        compute_GRB_properties = (GRB_properties is not None and 
                                  "compute_GRB_properties" in GRB_properties and 
                                  GRB_properties["compute_GRB_properties"])

        # drop all zero weights to save memory
        sel = w_ijk > 0.
        index = index[sel]
        z_formation = z_formation[sel]
        z_event = z_event[sel]
        w_ijk = w_ijk[sel]

        # export results, we do a selections of columns else the dataframe
        # risks to go out of memory
        df = pd.DataFrame()
        df['weight'] = w_ijk
        df['z_formation'] = z_formation
        if pop == 'DCO':
            df['z_merger'] = z_event
        elif pop == 'GRB':
            df['z_grb'] = z_event
        else:
            raise ValueError('Population not recognized!')
        save_cols = ['metallicity','time','t_delay','S1_state','S2_state',
                     'S1_mass','S2_mass','S1_spin','S2_spin',
                     'orbital_period','eccentricity', 'q', 'm_tot',
                     'm_chirp', 'chi_eff']
        if compute_GRB_properties:
            save_cols += GRB_PROPERTIES
        if "channel" in self.rates.df:
            save_cols.append('channel')
        if export_cols is not None:
            for c in export_cols:
                if c not in save_cols:
                    save_cols.append(c)
        for c in save_cols:
            df[c] = self.rates.get_data(c, index)

        # the same binary system can emit two GRBs, we remove the
        # GRB properties of the other GRB to prevent mistakes
        if reset_grb_properties is not None:
            if reset_grb_properties == 'GRB1':
                for key in GRB_PROPERTIES:
                    if "S1" in key:
                        df[key] = np.nan
            elif reset_grb_properties == 'GRB2':
                for key in GRB_PROPERTIES:
                    if "S2" in key:
                        df[key] = np.nan
            else:
                raise ValueError(f'reset_grb_properties=={reset_grb_properties} key not recognized!')

        return df

    def get_dco_merger_rate_density(self, export_cols=None,  working_dir='./', load_data=False):
        """Compute the merger rate density as a function of redshift.

        Parameters
        ----------
        export_cols : list str
            List of additional columns to save in the intrinsic/observable
            population.
        working_dir : str
            Working directory where the weights will be saved.
        load_data : bool
            `True` if you want to load the weights computed by this function
            in your working directory.

        """
        if self.df_synthetic is None:
            raise ValueError('You first need to isolated the DCO synthetic population!')

        # compute cosmological weights (detection rate weights with infinite sensitivity)
        sensitivity='infinite'
        flag_pdet = False
        index, z_formation, z_merger, w_ijk = self.compute_cosmological_weights(sensitivity, flag_pdet, working_dir=working_dir, load_data=load_data, pop='DCO')
        # compute rate density weights
        self.dco_z_rate_density = self.rates.get_centers_redshift_bins()
        total_rate = self.rates.compute_rate_density(w_ijk, z_merger, observable='DCO', sensitivity=sensitivity)
        self.dco_rate_density = {'total' : total_rate}
        if "channel" in self.df_synthetic:
            channels = np.unique(self.df_synthetic['channel'])
            for ch in tqdm(channels):
                sel = (self.rates.get_data('channel', index) == ch)
                rate = self.rates.compute_rate_density(w_ijk[sel], z_merger[sel], observable='DCO', sensitivity=sensitivity)
                self.dco_rate_density[ch] = rate

        print(f'DCO merger rate density in the local Universe (z={self.dco_z_rate_density[0]:1.2f}): {round(total_rate[0],2)} Gpc^-3 yr^-1')

        # export the intrinsic DCO population
        self.df_dco_intrinsic = self.resample_synthetic_population(index, z_formation, z_merger, w_ijk, export_cols=export_cols, pop='DCO')


    def get_grb_rate_density(self, export_cols=None,  working_dir='./', load_data=False):
        """Compute the GRB density as a function of redshift.

        Parameters
        ----------
        export_cols : list str
            List of additional columns to save in the intrinsic/observable
            population.
        working_dir : str
            Working directory where the weights will be saved.
        load_data : bool
            `True` if you want to load the weights computed by this function
            in your working directory.

        """
        if self.df_synthetic is None:
            raise ValueError('You first need to isolated the DCO synthetic population!')

        # compute cosmological weights (detection rate weights with infinite sensitivity)
        sensitivity='infinite'
        flag_pdet = False
        index_1, z_formation_1, z_grb_1, w_ijk_1, \
            index_2, z_formation_2, z_grb_2, w_ijk_2 = self.compute_cosmological_weights(sensitivity, flag_pdet, working_dir=working_dir, load_data=load_data, pop='GRB')

        # compute beamed rate density weights
        # note we compute the beamed rate density for GRBs as this is what is often reported in the literature
        # as the beaming factor is not known a priori
        self.grb_z_rate_density = self.rates.get_centers_redshift_bins()
        total_GRB1 = self.rates.compute_rate_density(w_ijk_1, z_grb_1, observable='GRB1', sensitivity='beamed', index=index_1)
        total_GRB2 = self.rates.compute_rate_density(w_ijk_2, z_grb_2, observable='GRB2', sensitivity='beamed', index=index_2)
        total_rate = total_GRB1 + total_GRB2
        self.grb_rate_density = {'total' : total_rate, 'total_GRB1' : total_GRB1, 'total_GRB2': total_GRB2}
        if "channel" in self.df_synthetic:
            channels = np.unique(self.df_synthetic['channel'])
            for ch in tqdm(channels):
                sel1 = (self.rates.get_data('channel', index_1) == ch)
                if any(sel1):
                    self.grb_rate_density[ch+'_GRB1'] = self.rates.compute_rate_density(w_ijk_1[sel1], z_grb_1[sel1], observable='GRB1', sensitivity='beamed', index=index_1[sel1])
                else:
                    self.grb_rate_density[ch+'_GRB1'] = np.zeros(len(self.grb_z_rate_density))
                sel2 = (self.rates.get_data('channel', index_2) == ch)
                if any(sel2):
                    self.grb_rate_density[ch+'_GRB2'] = self.rates.compute_rate_density(w_ijk_2[sel2], z_grb_2[sel2], observable='GRB2', sensitivity='beamed', index=index_2[sel2])
                else:
                    self.grb_rate_density[ch+'_GRB2'] = np.zeros(len(self.grb_z_rate_density))
                self.grb_rate_density[ch] = self.grb_rate_density[ch+'_GRB1'] + self.grb_rate_density[ch+'_GRB2']

        print(f'GRB (beamed) rate density in the local Universe (z={self.grb_z_rate_density[0]:1.2f}): {round(total_rate[0],2)} Gpc^-3 yr^-1')

        # export the intrinsic grb intrisic population
        # TODO: instead of concatenating two dataframe and duplicating the information of the same
        # binary system we should combine the two datraframes where we have two columns for GRB1 and GRB2
        # such dataframe will have z_grb_1, z_grb_2, weight_1, weight_2
        # when this is addressed, remove reset_grb_properties feature
        self.df_grb_intrinsic = pd.DataFrame()
        df_grb_1 = self.resample_synthetic_population(index_1, z_formation_1, z_grb_1, w_ijk_1, export_cols=export_cols, pop='GRB', reset_grb_properties='GRB2')
        df_grb_2 = self.resample_synthetic_population(index_2, z_formation_2, z_grb_2, w_ijk_2, export_cols=export_cols, pop='GRB', reset_grb_properties='GRB1')
        self.df_grb_intrinsic = pd.concat([df_grb_1, df_grb_2], ignore_index=True, sort=False)
        # the observable population accounts for beaming
        self.df_grb_observable = self.df_grb_intrinsic.copy()
        for i in [1,2]:
            sel = self.df_grb_observable[f'S{i}_f_beaming'] > 0
            self.df_grb_observable.loc[sel,'weight'] *= self.df_grb_observable.loc[sel,f'S{i}_f_beaming']

    def get_dco_detection_rate(self, sensitivity='design_H1L1V1', export_cols=None,  working_dir='./', load_data=False):
        """Compute the detection rate per yr.

        Parameters
        ----------
        sensitivity : str
            GW detector sensitivity and network configuration you want to use, see arXiv:1304.0670v3
            detector sensitivities are taken from: https://dcc.ligo.org/LIGO-T2000012-v2/public
            available sensitivity keys (for Hanford, Livingston, Virgo network):
                'O3actual_H1L1V1' : aligo_O3actual_H1.txt, aligo_O3actual_L1.txt, avirgo_O3actual.txt
                'O4low_H1L1V1' : aligo_O4low.txt, aligo_O4low.txt, avirgo_O4high_NEW.txt
                'O4high_H1L1V1' : aligo_O4high.txt, aligo_O4high.txt, avirgo_O4high_NEW.txt
                'design_H1L1V1' : AplusDesign.txt, AplusDesign.txt, avirgo_O5high_NEW.txt
        export_cols : list str
            List of additional columns to save in the intrinsic/observable
            population.
        working_dir : str
            Working directory where the weights will be saved.
        load_data : bool
            `True` if you want to load the weights computed by this function
            in your working directory.

        """
        if self.df_synthetic is None:
            raise ValueError('You first need to isolated the DCO synthetic population!')

        # compute detection rate weights
        flag_pdet = True
        index, z_formation, z_merger, w_ijk = self.compute_cosmological_weights(sensitivity, flag_pdet, working_dir=working_dir, load_data=load_data)
        print(f'DCO detection rate at {sensitivity} sensitivity: {sum(w_ijk):1.2f} yr^-1')

        # export the observable DCO population
        # TODO: store p_det
        self.df_dco_observable = self.resample_synthetic_population(index, z_formation, z_merger, w_ijk, export_cols=export_cols)
    
    def save_intrinsic_pop(self, path='./intrinsic_population_type.h5', pop='DCO'):
        """Save intrinsic population.

        Parameters
        ----------
        path : str
            Path to dataset.

        """
        if pop == 'DCO':
            if self.df_dco_intrinsic is None:
                raise ValueError('Nothing to save!')
            else:
                self.df_dco_intrinsic.to_hdf(path.replace('type', pop), key='history')
                if self.verbose:
                    print('Intrinsic population successfully saved!')
        elif pop == 'GRB':
            if self.df_grb_intrinsic is None:
                raise ValueError('Nothing to save!')
            else:
                self.df_grb_intrinsic.to_hdf(path.replace('type', pop), key='history')
                if self.verbose:
                    print('Intrinsic population successfully saved!')
        else:
            raise ValueError('Population not recognized!')

    def load_intrinsic_pop(self, path, pop='DCO'):
        """Load intrinsic population.

        Parameters
        ----------
        path : str
            Path to dataset.

        """
        if pop == 'DCO':
            if self.df_dco_intrinsic is None:
                self.df_dco_intrinsic = pd.read_hdf(path, key='history')
                if self.verbose:
                    print('Intrinsic population successfully loaded!')
            else:
                raise ValueError('You already have an intrinsic population stored in memory!')
        elif pop == 'GRB':
            if self.df_grb_intrinsic is None:
                self.df_grb_intrinsic = pd.read_hdf(path, key='history')
                if self.verbose:
                    print('Intrinsic population successfully loaded!')
            else:
                raise ValueError('You already have an intrinsic population stored in memory!')

    def save_observable_pop(self, path='./observable_population_type.h5', pop='DCO'):
        """Save observable population.

        Parameters
        ----------
        path : str
            Path to dataset.

        """
        if pop == 'DCO':
            if self.df_dco_observable is None:
                raise ValueError('Nothing to save!')
            else:
                self.df_dco_observable.to_hdf(path.replace('type', pop), key='history')
                if self.verbose:
                    print('observable population successfully saved!')
        elif pop == 'GRB':
            if self.df_grb_observable is None:
                raise ValueError('Nothing to save!')
            else:
                self.df_grb_observable.to_hdf(path.replace('type', pop), key='history')
                if self.verbose:
                    print('observable population successfully saved!')
        else:
            raise ValueError('Population not recognized!')

    def load_observable_pop(self, path, pop='DCO'):
        """Load observable population.

        Parameters
        ----------
        path : str
            Path to dataset.

        """
        if pop == 'DCO':
            if self.df_dco_observable is None:
                self.df_dco_observable = pd.read_hdf(path, key='history')
                if self.verbose:
                    print('observable population successfully loaded!')
            else:
                raise ValueError('You already have an observable population stored in memory!')
        elif pop == 'GRB':
            if self.df_grb_observable is None:
                self.df_grb_observable = pd.read_hdf(path, key='history')
                if self.verbose:
                    print('observable population successfully loaded!')
            else:
                raise ValueError('You already have an observable population stored in memory!')
        else:
            raise ValueError('Population not recognized!')

    def plot_hist_properties(self, var, intrinsic=False, observable=False, pop=None, **kwargs):
        """Plot histogram of intrinsic/observable properites.

        Parameters
        ----------
        var : str
            Property to plot stored in intrinsic/observable dataframe.
        intrinsic : bool
            `True` if you want to deplay the intrisc population.
        observable : bool
            `True` if you want to deplay the observable population.
        **kwargs : dict
            ploting arguments

        """
        if pop == 'DCO':
            if self.df_dco_intrinsic is None and self.df_dco_observable is None:
                raise ValueError('First you need to compute the merger rate density!')
            if intrinsic:
                df_intrinsic = self.df_dco_intrinsic
            else:
                df_intrinsic = None
            if observable:
                df_observable = self.df_dco_observable
            else:
                df_observable = None
        elif pop == 'GRB':
            if self.df_grb_intrinsic is None and self.df_grb_observable is None:
                raise ValueError('First you need to compute the merger rate density!')
            if intrinsic:
                df_intrinsic = self.df_grb_intrinsic
            else:
                df_intrinsic = None
            if observable:
                df_observable = self.df_grb_observable
            else:
                df_observable = None
        else:
            raise ValueError('Population not recognized!')
        plot_pop.plot_hist_properties(var, df_intrinsic=df_intrinsic, df_observable=df_observable, pop=pop, **kwargs)

    def plot_rate_density(self, DCO=False, GRB=False, **kwargs):
        """Plot DCO and GRB rate densities."""
        if not DCO and not GRB:
            raise ValueError('You need to choose at least one population to plot!')
        if DCO:
            if self.dco_z_rate_density is None or self.dco_rate_density is None:
                raise ValueError('First you need to compute the merger rate density!')
            else:
                z_dco = self.dco_z_rate_density
                rate_dco =  self.dco_rate_density
        else:
            z_dco = None
            rate_dco = None
        if GRB:
            if self.grb_z_rate_density is None or self.grb_rate_density is None:
                raise ValueError('First you need to compute the GRB rate density!')
            else:
                z_grb = self.grb_z_rate_density
                rate_grb =  self.grb_rate_density
        else:
            z_grb = None
            rate_grb = None
        plot_pop.plot_rate_density(z_dco, rate_dco, z_grb, rate_grb, **kwargs)

    def plot_popsyn_over_grid_slice(self, grid_type, met_Zsun, **kwargs):
        """Plot popsyn over grid slice."""
        plot_pop.plot_popsyn_over_grid_slice(self, grid_type, met_Zsun, **kwargs)
