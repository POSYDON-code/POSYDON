"""Processing of population files.

This module contains classes and functions to process population files.
Population files are HDF5 files containing the history and oneline dataframes
of a population of binary systems. The history dataframe contains the detailed
evolution of each binary system, while the oneline dataframe contains the
final state of each binary system.

Classes
-------
Population
    A class to handle population files.
History
    A class to handle the history dataframe of a population file.
Oneline
    A class to handle the oneline dataframe of a population file.
    
TransientPopulation
    A class to handle transient populations.

Rates
    A class to handle the cosmic rates of in a population file.
"""

__authors__ = [
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Kyle Akira Rocha <kylerocha2024@u.northwestern.edu>",
    "Monica Gallegos-Garcia <monicagallegosgarcia2024@u.northwestern.edu>",
    "Max Briel <max.briel@unige.ch>",
]

import warnings
import numpy as np
import pandas as pd
from tqdm import tqdm
import os
from matplotlib import pyplot as plt

from posydon.utils.constants import Zsun
from posydon.popsyn.io import binarypop_kwargs_from_ini
from posydon.popsyn.normalized_pop_mass import initial_total_underlying_mass
import posydon.visualization.plot_pop as plot_pop
from posydon.utils.common_functions import convert_metallicity_to_string

from astropy.cosmology import Planck15 as cosmology
from astropy.cosmology import z_at_value
from scipy.interpolate import interp1d
from astropy import units as u
from astropy import constants as const
import scipy as sp


from posydon.popsyn.rate_calculation import DEFAULT_MODEL

from posydon.popsyn.star_formation_history import star_formation_rate, SFR_Z_fraction_at_given_redshift

from posydon.popsyn.binarypopulation import (BinaryPopulation,
                                             HISTORY_MIN_ITEMSIZE,
                                             ONELINE_MIN_ITEMSIZE)
    

###############################################################################

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

class PopulationRunner():
    
    def __init__(self, path_to_ini, verbose=False):
        '''initialise the binary populations'''
        if '.ini' not in path_to_ini:
            raise ValueError('You did not provide a valid path_to_ini!')
        else:
            self.synthetic_pop_params = binarypop_kwargs_from_ini(path_to_ini)
            self.solar_metallicities = self.synthetic_pop_params['metallicity']
            self.verbose = verbose
            if not isinstance( self.solar_metallicities, list):
                self.solar_metallicities = [self.solar_metallicities]
                
            self.binary_populations = []
            for MET in self.solar_metallicities:
                ini_kw = binarypop_kwargs_from_ini(path_to_ini)
                ini_kw['metallicity'] = MET
                ini_kw['temp_directory'] = convert_metallicity_to_string(MET) + "_Zsun_" + ini_kw['temp_directory']
                self.binary_populations.append(BinaryPopulation(**ini_kw))
                
    def evolve(self):
        '''evolve the binary populations'''
        for pop in self.binary_populations:
            pop.evolve()
            # get the files in the batches
            if pop.comm is None:
                self.merge_parallel_runs(pop)
                    
    def merge_parallel_runs(self, pop):
        path_to_batch = pop.kwargs['temp_directory']
        tmp_files = [os.path.join(path_to_batch, f)    \
                             for f in os.listdir(path_to_batch) \
                                if os.path.isfile(os.path.join(path_to_batch, f))]
        pop.combine_saved_files(convert_metallicity_to_string(pop.metallicity) + '_Zsun_population.h5', tmp_files)
        # remove files
        if len(os.listdir(path_to_batch)) == 0:
            os.rmdir(path_to_batch)
            

# TODO: I don't know if the merge_populations function is still needed
def merge_populations(populations, filename, verbose=False):
    '''merges multiple populations into a single population file'''
    
    
    history_cols = populations[0].history.columns
    oneline_cols = populations[0].oneline.columns
    
    history_min_itemsize = {key: val for key, val in
                                HISTORY_MIN_ITEMSIZE.items()
                                if key in history_cols}
    oneline_min_itemsize = {key: val for key, val in
                                ONELINE_MIN_ITEMSIZE.items()
                                if key in oneline_cols}
    
    # open new file
    with pd.HDFStore(filename, mode='a') as store:
        
        prev = 0            
        for pop in populations:

            # write history
            for i in tqdm(range(0, len(pop), 10000), total=len(pop)//10000, disable=not verbose):
                tmp_df = pop.history[i:i+10000]
                tmp_df.index = tmp_df.index + prev
                store.append('history', tmp_df, format='table', data_columns=True, min_itemsize=history_min_itemsize)    
                
            if verbose:
                print(f'History for {pop.filename} written to population file!')
                
            # write oneline
            for i in tqdm(range(0, len(pop), 10000), total=len(pop)//10000, disable=not verbose):
                tmp_df = pop.oneline[i:i+10000]
                tmp_df.index = tmp_df.index + prev
                store.append('oneline', tmp_df, format='table', data_columns=True, min_itemsize=oneline_min_itemsize)
            
            if verbose:
                print(f'Oneline for {pop.filename} written to population file!')
                
            prev += len(pop)

        # merge mass_per_met
        tmp_df = pd.concat([pop.mass_per_met for pop in populations])
        mass_per_met = tmp_df.groupby(tmp_df.index).sum()
        
        store.put('mass_per_met', mass_per_met)
        
    # add ini parameters from the first population
    populations[0]._save_ini_params(filename)

    print(f'Populations merged into {filename}!')
    
##################
# Helper classes #
################## 

class History():
    
    def __init__(self, filename, verbose=False, chunksize=10000):
        self.filename = filename
        self.verbose = verbose
        self.chunksize = chunksize
        self.lengths = None
        self.number_of_systems = None
        self.columns = None
        
        # add history_lengths
        with pd.HDFStore(filename, mode='a') as store:
            # get the history lengths from the file
            if '/history_lengths' in store.keys():
                self.lengths = store['history_lengths']
            else:
                if self.verbose:
                    print('history_lengths not found in population file. Calculating history lengths...')
                history_events = store.select_column('history', 'index')
                self.lengths = pd.DataFrame(history_events.groupby(history_events).count())
                if self.verbose:
                    print('Storing history lengths in population file!')
                store.put('history_lengths', pd.DataFrame(self.lengths), format='table')
                
                del history_events
            
            self.columns = store.select('history', start=0, stop=0).columns
        
        self.indices = self.lengths.index.to_numpy()
        self.number_of_systems = len(self.lengths)
        
    def __getitem__(self, key):
        '''Return the history table'''
        if isinstance(key, slice):
            if key.start is None:
                pre = 0
            else: 
                pre = self.lengths.loc[:key.start-1].sum()
            if key.stop is None:
                chunk = self.lengths.loc[key.start:].sum()
            else:
                chunk = self.lengths.loc[key.start:key.stop-1].sum()
            return pd.read_hdf(self.filename, key='history', start=pre, stop=pre+chunk)
        elif isinstance(key, int):
            return pd.read_hdf(self.filename, where='index == key', key='history')
        elif isinstance(key, list) and all(isinstance(x, int) for x in key):
            with pd.HDFStore(self.filename, mode='r') as store:
                return store.select('history', where='index in key')
        elif isinstance(key, np.ndarray) and (key.dtype == int):
            indices = self.lengths[key].index.tolist()
            with pd.HDFStore(self.filename, mode='r') as store:
                return store.select('history', where='index in indices')
        elif isinstance(key, np.ndarray) and key.dtype == bool:
            out_df = pd.DataFrame()
            for i in range(0,len(key), self.chunksize):
                tmp_df = pd.read_hdf(self.filename, key='history', start=i, stop=i+self.chunksize)[key[i:i+self.chunksize]]
                out_df = pd.concat([out_df, tmp_df])
            return out_df
            
        elif isinstance(key, str):
            if key in self.columns:
                return pd.read_hdf(self.filename, key='history', columns=[key])
            else:
                raise ValueError(f'{key} is not a valid column name!')
        elif isinstance(key, list) and all(isinstance(x, str) for x in key):
            if all(x in self.columns for x in key):
                return pd.read_hdf(self.filename, key='history', columns=key)
            else:
                raise ValueError(f'Not all columns in {key} are valid column names!')
        else:
            raise ValueError('Invalid key type!')            
        
    def __len__(self):
        return np.sum(self.lengths)
    
    def head(self, n=10):
        '''Return the first n rows of the history table'''
        return pd.read_hdf(self.filename, key='history', start=0, stop=n)
        
    def tail(self, n=10):
        '''Return the last n rows of the history table'''
        return pd.read_hdf(self.filename, key='history', start=-n)
    
    def __repr__(self):
        return pd.read_hdf(self.filename, key='history').__repr__()
    
    def _repr_html_(self):
        return pd.read_hdf(self.filename, key='history')._repr_html_()
    
    def iterator(self):
        return pd.read_hdf(self.filename, key='history', iterator=True, chunksize=self.chunksize)
        
    def iterbinaries(self):
        for i in self.indices:
            yield pd.read_hdf(self.filename, key='history', where='index == i')   
         
    def select(self, where=None, start=None, stop=None, columns=None):
        with pd.HDFStore(self.filename, mode='r') as store:
            return store.select('history', where=where, start=start, stop=stop, columns=columns)

class Oneline():
    
    def __init__(self, filename, verbose=False, chunksize=10000):
        self.filename = filename
        self.verbose = verbose
        self.chunksize = chunksize
        self.number_of_systems = None
        self.indices = None
                
        with pd.HDFStore(filename, mode='r') as store:
            self.indices = store.select_column('oneline', 'index')
            self.columns = store.select('oneline', start=0, stop=0).columns
            
        self.number_of_systems = len(self.indices)
        
    def __getitem__(self, key):
        '''Return the oneline table'''
        if isinstance(key, slice):
            if key.start is None:
                pre = 0
            else:
                pre = key.start
            if key.stop is None:
                chunk = self.number_of_systems
            else:
                chunk = key.stop - key.start
            return pd.read_hdf(self.filename, key='oneline', start=pre, stop=pre+chunk)
        elif isinstance(key, int):
            return pd.read_hdf(self.filename, where='index == key', key='oneline')
        elif isinstance(key, list) and all(isinstance(x, int) for x in key):
            return pd.read_hdf(self.filename, where='index in key', key='oneline')
        elif isinstance(key, np.ndarray) and key.dtype == bool:
            indices = self.lengths[key].index
            return pd.read_hdf(self.filename, where='index in indices', key='oneline')
        elif isinstance(key, str):
            if key in self.columns:
                return pd.read_hdf(self.filename, key='oneline', columns=[key])
            else:
                raise ValueError(f'{key} is not a valid column se!')
            
        elif isinstance(key, list) and all(isinstance(x, str) for x in key):
            if all(x in self.columns for x in key):
                return pd.read_hdf(self.filename, key='oneline', columns=key)
            else:
                raise ValueError(f'Not all columns in {key} are valid column names!')
        else:
            raise ValueError('Invalid key type!')
        
    def __len__(self):
        return self.number_of_systems
    
    def head(self, n=10):
        '''Return the first n rows of the oneline table'''
        return pd.read_hdf(self.filename, key='oneline', start=0, stop=n)
    
    def tail(self, n=10):
        '''Return the last n rows of the oneline table'''
        return pd.read_hdf(self.filename, key='oneline', start=-n)
    
    def __repr__(self):
        return pd.read_hdf(self.filename, key='oneline').__repr__()
    
    def _repr_html_(self):
        return pd.read_hdf(self.filename, key='oneline')._repr_html_()

    def iterator(self):
        return pd.read_hdf(self.filename, key='oneline', iterator=True, chunksize=self.chunksize)
    
    def iterbinaries(self):
        for i in self.indices:
            yield pd.read_hdf(self.filename, key='oneline', where='index == i')

    def select(self,where=None, start=None, stop=None, columns=None):
        with pd.HDFStore(self.filename, mode='r') as store:
            return store.select('oneline', where=where, start=start, stop=stop, columns=columns)


class PopulationIO():
    def __init__(self):
        pass
    
    def _load_metadata(self, filename):
        '''Load the metadata from the file'''
        if '.h5' not in filename:
            raise ValueError(f'{filename} does not contain .h5 in the se.\n Is this a valid population file?')
        
        self._load_ini_params(filename)
        self._load_mass_per_met(filename)
    
    def _save_mass_per_met(self, filename):
        with pd.HDFStore(filename, mode='a') as store:
            store.put('mass_per_met', self.mass_per_met)
            if self.verbose:
                print('mass_per_met table written to population file!')
        
    def _load_mass_per_met(self, filename):
        with pd.HDFStore(filename, mode='r') as store:
            self.mass_per_met = store['mass_per_met']
            if self.verbose:
                print('mass_per_met table read from population file!')

    def _save_ini_params(self, filename):
        '''
        Store the ini parameters to the ParsedPopulation file'''
        with pd.HDFStore(filename, mode='a') as store:
            # write ini parameters to file
            tmp_df = pd.DataFrame()
            for c in parameter_array:
                tmp_df[c] = [self.ini_params[c]]
            store.put('ini_parameters', tmp_df)
            
    def _load_ini_params(self, filename):
        # load ini parameters
        with pd.HDFStore(filename, mode='r', ) as store:
            tmp_df = store['ini_parameters']
            self.ini_params = {}
            for c in parameter_array:
                self.ini_params[c] = tmp_df[c][0]

##########################
# Main interface classes #
##########################

class Population(PopulationIO):
    
    def __init__(self, filename, metallicity=None, ini_file=None, verbose=False, chunksize=1000):
        '''A Population file.
        
        Parameters
        -----------
        filename : str
            The path to the population file
        verbose : bool
            If `True`, print additional information.
        chunksize : int
            The chunksize to use when reading the population file.

        Optional parameters for adding additional information to old population files.            
        metallicity : float
            The metallicity of the population in solar units.
        ini_file : str
            The path to the ini file used to create the population. 
        
        '''
        
        self.filename = filename
        self.verbose = verbose
        self.chunksize = chunksize
        
        self.mass_per_met = None
        self.number_of_systems = None
        self.history_lengths = None

        # if the user provided a single string instead of a list of strings
        if not ('.h5' in filename):
            raise ValueError(f'{filename} does not contain .h5 in the se.\n Is this a valid population file?')
        
        # read the population file
        with pd.HDFStore(filename, mode='r') as store:
            keys = store.keys()
        
        # check if pop contains history
        if '/history' not in keys:
            raise ValueError(f'{filename} does not contain a history table!')
        else:
            self.history = History(filename, self.verbose, self.chunksize)
        
        # check if pop contains oneline
        if '/oneline' not in keys:
            raise ValueError(f'{filename} does not contain an oneline table!')
        else:
            self.oneline = Oneline(filename, self.verbose, self.chunksize)
        
        # check if formation channels are present
        if '/formation_channels' not in keys:
            if self.verbose:
                warnings.warn(f'{filename} does not contain formation channels!')
            self._formation_channels = None
        else:
            self._formation_channels = pd.read_hdf(self.filename, key='formation_channels')
        
        # if an ini file is given, read the parameters from the ini file
        if ini_file is not None:
            self.ini_params = binarypop_kwargs_from_ini(ini_file)
            self._save_ini_params(filename)
            self._load_ini_params(filename)
        else:
            if '/ini_parameters' not in keys:
                raise ValueError(f'{filename} does not contain an ini_parameters table!')
            else:
                self._load_ini_params(filename)
        
        # check if pop contains mass_per_met table
        if '/mass_per_met' in keys and metallicity is None:
            self._load_mass_per_met(filename)
            self.solar_metallicities = self.mass_per_met.index.to_numpy()
            self.metallicities = self.solar_metallicities * Zsun
        elif metallicity is None:
            raise ValueError(f'{filename} does not contain a mass_per_met table and no metallicity for the file was given!')
        
        # calculate the metallicity information. This assumes the metallicity is for the whole file!    
        if metallicity is not None and ini_file is not None:
            if '/mass_per_met' in keys:
                warnings.warn(f'{filename} already contains a mass_per_met table. Overwriting the table!')
                
            simulated_mass = np.sum(self.oneline[['S1_mass_i', 'S2_mass_i']].to_numpy())
            underlying_mass = initial_total_underlying_mass(df=simulated_mass, **self.ini_params)[0]
            self.mass_per_met = pd.DataFrame(index=[metallicity], 
                                              data={'simulated_mass':simulated_mass,
                                                    'underlying_mass':underlying_mass,
                                                    'number_of_systems':len(self.oneline)})
            
            self._save_mass_per_met(filename)
            self.solar_metallicities = self.mass_per_met.index.to_numpy()
            self.metallicities = self.solar_metallicities * Zsun
            
        elif metallicity is not None and ini_file is None:
            raise ValueError(f'{filename} does not contain a mass_per_met table and no ini file was given!')        
        
        # add number of systems
        self.history_lengths = self.history.lengths
        self.number_of_systems = self.oneline.number_of_systems
        self.indices = self.history.indices
        
    def export_selection(self, selection, filename, chunksize=1):
        '''Export a selection of the population to a new file
        
        Parameters
        -----------
        selection  : list of int
            The indices of the systems to export    
        filename   : str
            The se of the export file to create or append to
        '''
        
        if not ('.h5' in filename):
            raise ValueError(f'{filename} does not contain .h5 in the se.\n Is this a valid population file?')
        
        history_cols = self.history.columns
        oneline_cols = self.oneline.columns
    
        history_min_itemsize = {key: val for key, val in
                                HISTORY_MIN_ITEMSIZE.items()
                                if key in history_cols}
        oneline_min_itemsize = {key: val for key, val in
                                ONELINE_MIN_ITEMSIZE.items()
                                if key in oneline_cols}
        
        
        with pd.HDFStore(filename, mode='a') as store:
            # shift all new indices by the current length of data in the file
            last_index_in_file = 0
            
            if '/oneline' in store.keys():
                last_index_in_file = np.sort(store['oneline'].index)[-1]
            elif '/history' in store.keys():
                last_index_in_file = np.sort(store['history'].index)[-1]
            
            if '/history' in store.keys() and self.verbose:
                print('history in file. Appending to file')
                
            if '/oneline' in store.keys() and self.verbose:
                print('oneline in file. Appending to file')
                
            if '/formation_channels' in store.keys() and self.verbose:
                print('formation_channels in file. Appending to file')
                
            if '/history_lengths' in store.keys() and self.verbose:
                print('history_lengths in file. Appending to file')
            
            # TODO: I need to shift the indices of the binaries or should I reindex them?
            # since I'm storing the information, reindexing them should be fine. 
            
            if last_index_in_file == 0:
                reindex = {i:j for i,j in zip(selection, 
                                          np.arange(last_index_in_file, last_index_in_file+len(selection), 1)
                                   )}
            else:
                reindex = {i:j for i,j in zip(selection, 
                                          np.arange(last_index_in_file+1, last_index_in_file+len(selection)+1, 1)
                                   )}
                       
            if 'metallicity' not in self.oneline.columns:
                warnings.warn('No metallicity column in oneline dataframe! Using the metallicity of the population file and adding it to the oneline.')
                if len(self.metallicities) > 1:
                    raise ValueError('The population file contains multiple metallicities. Please add a metallicity column to the oneline dataframe!')
            
            if self.verbose:
                print('Writing selected systems to population file...')
            
            # write oneline of selected systems
            for i in tqdm(range(0, len(selection), 1000), total=len(selection)//1000, disable=not self.verbose):
                tmp_df = self.oneline[selection[i:i+1000]]
                if 'metallicity' in tmp_df.columns:
                    tmp_df['metallicity'] = tmp_df['metallicity'].astype('float')
                else:
                    tmp_df['metallicity'] = self.metallicities[0]
                tmp_df.rename(index=reindex, inplace=True)
                store.append('oneline', tmp_df, format='table', data_columns=True, min_itemsize=oneline_min_itemsize, index=False)
            
            if self.verbose:
                print('Oneline: Done')

            # write history of selected systems    
            for i in tqdm(range(0,len(selection), chunksize), total=len(selection)//chunksize, disable=not self.verbose):
                tmp_df = self.history[selection[i:i+chunksize]]
                tmp_df.rename(index=reindex, inplace=True)
                store.append('history', tmp_df, format='table', data_columns=True, min_itemsize=history_min_itemsize, index=False)
        
            if self.verbose:
                print('History: Done')
        
            # write formation channels of selected systems
            if self.formation_channels is not None:
                for i in tqdm(range(0, len(selection), 10000), total=len(selection)//10000, disable=not self.verbose):
                    tmp_df = self.formation_channels.loc[selection[i:i+10000]]
                    tmp_df.rename(index=reindex, inplace=True)
                    store.append('formation_channels', tmp_df, format='table', data_columns=True, min_itemsize={'channel_debug': 100, 'channel': 100})
            
            
            ## METADATA
            
            # write the history lengths
            for i in tqdm(range(0, len(selection), 10000), total=len(selection)//10000, disable=not self.verbose):
                tmp_df = self.history.lengths.loc[selection[i:i+10000]]
                tmp_df.rename(index=reindex, inplace=True)
                store.append('history_lengths', pd.DataFrame(tmp_df), format='table', index=False)
            
            # write mass_per_met
            if '/mass_per_met' in store.keys():
                self_mass = self.mass_per_met
                self_mass['number_of_systems'] = len(selection)
                tmp_df = pd.concat([store['mass_per_met'], self_mass])
                mass_per_met = tmp_df.groupby(tmp_df.index).sum()
                store.put('mass_per_met', mass_per_met)
        
            else:
                self_mass = self.mass_per_met
                self_mass['number_of_systems'] = len(selection)
                store.put('mass_per_met', self_mass)
        
        # write ini parameters
        self._save_ini_params(filename)
    
    @property
    def formation_channels(self):
        '''Return the formation channels of the population'''
        if self._formation_channels is None:
            with pd.HDFStore(self.filename, mode='r') as store:
                if '/formation_channels' in store.keys():
                    self._formation_channels = pd.read_hdf(self.filename, key='formation_channels')
                else:
                    warnings.warn('No formation channels in the population file!')
                    self._formation_channels = None
            
        return self._formation_channels

    def calculate_formation_channels(self, mt_history=False):
        
        if self.verbose: print('Calculating formation channels...')  
        
        # load the HMS-HMS interp class
        HMS_HMS_event_dict = {'initial_MT'  : 'initial_MT',
                              'stable_MT'   : 'oRLO1', 
                              'no_MT'       : 'None', 
                              'unstable_MT' : 'oCE1/oDoubleCE1'}
        
        unique_binary_indices = self.indices
        
        # check if formation channels already exist
        with pd.HDFStore(self.filename, mode='a') as store:
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
    
        previous = 0
        
        for i in tqdm(range(0,len(unique_binary_indices), self.chunksize), disable=not self.verbose):
            selection = unique_binary_indices[i:i+self.chunksize]
            
            # create the dataframe for the chunk
            df = pd.DataFrame(index=selection, columns=['channel_debug', 'channel'])
            end = previous + self.history_lengths.iloc[i:i+self.chunksize].sum().iloc[0]

            # get the history of chunk events and transform the interp_class_HMS_HMS
            interp_class_HMS_HMS = self.oneline.select(start=i, stop=i+self.chunksize, columns=['interp_class_HMS_HMS'])
            events = self.history.select(start=previous, stop=end, columns=['event'])
            
            mask = ~ pd.isna(interp_class_HMS_HMS['interp_class_HMS_HMS'].values)
            interp_class_HMS_HMS.loc[mask, 'interp_class_HMS_HMS'] = interp_class_HMS_HMS[mask].apply(lambda x: HMS_HMS_event_dict[x['interp_class_HMS_HMS']], axis=1).values
            del mask
            
            previous = end
            # combine based on the index, this allows for an easier apply later
            merged = pd.merge(events.dropna(), interp_class_HMS_HMS, left_index=True, right_index=True)
            del events, interp_class_HMS_HMS
            
            merged.index.name ='binary_index'
            df['channel_debug'] = merged.groupby('binary_index').apply(get_events)
            del merged
            df['channel'] = df['channel_debug'].str.replace('_redirect', '').str.replace('_CO_contact', '')

            if mt_history:
                columns = self.oneline.columns
                if 'mt_history_HMS_HMS' not in columns:
                    raise ValueError('mt_history_HMS_HMS not saved in the oneline dataframe!')
                else:
                    tmp_df = pd.DataFrame(index=selection, columns=['channel', 'mt_history_HMS_HMS'])
                    tmp_df['channel'] = df['channel']
                    x = self.oneline.select(start=i, stop=i+self.chunksize, columns=['mt_history_HMS_HMS'])
                    tmp_df['mt_history_HMS_HMS'] = x
                    df['channel'] = tmp_df.apply(mt_history, axis=1)
                    del tmp_df
                    del x
                            
            self._write_formation_channels(self.filename, df)
            del df

    def _write_formation_channels(self, filename, df):
        '''Write the formation channels to the population file'''
        with pd.HDFStore(filename, mode='a') as store:
            store.append('formation_channels', df, format='table', data_columns=True, min_itemsize={'channel_debug': 100, 'channel': 100})
            if self.verbose:
                print('formation_channels written to population file!')

    def __len__(self):
        return self.number_of_systems

    @property
    def columns(self):
        return {'history':self.history.columns,
                'oneline':self.oneline.columns}
        
    def create_transient_population(self, func, transient_name, oneline_cols=None, hist_cols=None):
        '''Given a function, create a synthetic population
        
        Parameters
        ----------
        func : function
            Function to apply to the parsed population to create the synthetic population.
            The function needs to take 3 arguments:
                - history_chunk : pd.DataFrame
                - oneline_chunk : pd.DataFrame
                - formation_channels_chunk : pd.DataFrame
                and return a pd.DataFrame containing the synthetic population, which needs to contain a column 'time'.
                
        oneline_cols : list of str
            Columns to extract from the oneline dataframe. default is all columns.
        hist_cols : list of str
            Columns to extract from the history dataframe. default is all columns.
            
        Returns
        -------
        SyntheticPopulation
            A synthetic population containing the synthetic population. 
        '''
        
        with pd.HDFStore(self.filename, mode='a') as store:
            if f'/transients/{transient_name}' in store.keys():
                print('overwriting transient population')
                del store['transients/'+transient_name]

        min_itemsize = {'channel': 100,}
        if hist_cols is not None:
            if 'time' not in hist_cols:
                raise ValueError('The transient population requires a time column!')
        
            min_itemsize.update({key:val for key, val in 
                                HISTORY_MIN_ITEMSIZE.items() 
                                if key in hist_cols})
        else:
            hist_cols = self.history.columns
            min_itemsize.update(HISTORY_MIN_ITEMSIZE)
        
        if oneline_cols is not None:
            min_itemsize.update({key:val for key, val in
                                ONELINE_MIN_ITEMSIZE.items()
                            if key in oneline_cols})
        else:
            oneline_cols = self.oneline.columns
            min_itemsize.update(ONELINE_MIN_ITEMSIZE)

        # setup a mapping to the size of each history colummn
        history_lengths = self.history_lengths
        unique_binary_indices = self.indices
        
        previous = 0
        for i in tqdm(range(0,len(unique_binary_indices), self.chunksize), disable=not self.verbose):

            end = previous + history_lengths[i:i+self.chunksize].sum().iloc[0]
            
            oneline_chunk = self.oneline.select(start=i,
                                        stop=i+self.chunksize,
                                        columns=oneline_cols)
            
            history_chunk = self.history.select(start=previous,
                                                stop=end, 
                                                columns=hist_cols)
            
            if self.formation_channels is not None:
                formation_channels_chunk = self.formation_channels[i:i+self.chunksize]
            else:
                formation_channels_chunk = None
            
            syn_df = func(history_chunk, oneline_chunk, formation_channels_chunk)

            if len(syn_df.columns) != len(syn_df.columns.unique()):
                raise ValueError('Synthetic population contains duplicate columns!')
            
            # filter out the columns in min_itemsize that are not in the dataframe
            min_itemsize = {key:val for key, val in min_itemsize.items() if key in syn_df.columns}
            
            with pd.HDFStore(self.filename, mode='a') as store:
                store.append('transients/'+transient_name,
                             syn_df,
                             format='table',
                             data_columns=True,
                             min_itemsize=min_itemsize
                             )
            
            previous = end
        
        synth_pop = TransientPopulation(self.filename, transient_name, verbose=self.verbose)
        return synth_pop

    
    def plot_binary_evolution(self, index):
        '''Plot the binary evolution of a system'''
        pass

class TransientPopulation(Population):
    
    def __init__(self, filename, transient_name, verbose=False):
        '''This class contains a synthetic population of transient events.

        You can calculate additional properties of the population, such as
        the formation channel, merger times, GRB properties, etc.
        
        pop_file : str
            Path to the synthetic population file.
        verbose : bool
            If `True`, print additional information.
        ''' 
        super().__init__(filename, verbose=verbose)
        
        with pd.HDFStore(self.filename, mode='r') as store:
            if '/transients/'+transient_name not in store.keys():
                raise ValueError(f'{transient_name} is not a valid transient population in {filename}!')

            self.transient_name = transient_name        
            if '/transients/'+transient_name+'/efficiencies' in store.keys():
                self._load_efficiency(filename)
                
    @property
    def population(self):
        '''Returns the whole synthetic populaiton as a pandas dataframe.
        
        Warning: might be too big to load in memory!
        
        Returns
        -------
        pd.DataFrame
            Dataframe containing the synthetic population.
        '''
        return pd.read_hdf(self.filename, key='transients/'+self.transient_name)        

    def _load_efficiency(self, filename):
        '''Load the efficiency from the file'''
        with pd.HDFStore(filename, mode='r') as store:
            self.efficiency = store['transients/'+self.transient_name+'/efficiencies']
            if self.verbose:
                print('Efficiency table read from population file!')
                
    def _save_efficiency(self, filename):
        '''Save the efficiency to the file'''
        with pd.HDFStore(filename, mode='a') as store:
            store.put('transients/'+self.transient_name+'/efficiencies', self.efficiency)
    
    @property
    def columns(self):
        '''Return the columns of the transient population'''
        if not hasattr(self, '_columns'):
            with pd.HDFStore(self.filename, mode='r') as store:
                self._columns = store.select('transients/'+self.transient_name, start=0, stop=0).columns
        return self._columns
    
        
    def select(self, where=None, start=None, stop=None, columns=None):
        '''Select a subset of the transient population'''
        return pd.read_hdf(self.filename, key='transients/'+self.transient_name, where=where, start=start, stop=stop, columns=columns)
        
    def get_efficiency_over_metallicity(self):
        """Compute the efficiency of events per Msun for each solar_metallicities."""
        
        if hasattr(self, 'efficiency'):
            print('Efficiencies already computed! Overwriting them!')
            with pd.HDFStore(self.filename, mode='a') as store:
                if '/transients/'+self.transient_name+'/efficiencies' in store.keys():
                    del store['transients/'+self.transient_name+'/efficiencies']
                
        metallicities = self.mass_per_met.index.to_numpy()
        efficiencies = []
        for met in metallicities:
            count = self.select(where='metallicity == {}'.format(met)).shape[0]
            
            # just sums the number of events
            underlying_stellar_mass = self.mass_per_met['underlying_mass'][met]
            
            eff = count/underlying_stellar_mass
            efficiencies.append(eff)
            print(f'Efficiency at Z={met:1.2E}: {eff:1.2E} Msun^-1')
        
        self.efficiency = pd.DataFrame(index=metallicities, data={'total':np.array(efficiencies)})

        # if the channel column is present compute the merger efficiency per channel
        if 'channel' in self.columns:
            channels = np.unique(self.select(columns=['channel']).values)
            for ch in channels:
                efficiencies = []
                for met in metallicities:
                    count = np.sum((self.select(where='(metallicity == {})'.format(met), columns=['channel']) == ch).values)
                    if count > 0:
                        underlying_stellar_mass = self.mass_per_met['underlying_mass'][met]
                        eff = count/underlying_stellar_mass
                    else:
                        eff = 0
                    efficiencies.append(eff)
                self.efficiency[ch] = np.array(efficiencies)
        
        # save the efficiency
        self._save_efficiency(self.filename)


    def calculate_cosmic_weights(self, SFH_identifier, MODEL_in=None):
        '''Calculate the cosmic weights of the transient population'''
        
        # Set model to DEFAULT or provided MODEL parameters
        # Allows for partial model specification
        if MODEL_in is None:
            MODEL = DEFAULT_MODEL
        else:
            for key in MODEL_in:
                if key not in DEFAULT_MODEL:
                    raise ValueError(key + " is not a valid parameter name!")
            
            # write the DEFAULT_MODEL with updates parameters to self.MODEL.
            MODEL = DEFAULT_MODEL
            MODEL.update(MODEL_in)
        
        path_in_file = '/transients/'+self.transient_name+'/rates/'+SFH_identifier+'/'
        
        with pd.HDFStore(self.filename, mode='a') as store:
            if path_in_file+'MODEL' in store.keys():
                store.remove(path_in_file+'MODEL')
                if self.verbose:
                    print('Cosmic weights already computed! Overwriting them!')
                if path_in_file+'weights' in store.keys():
                    store.remove(path_in_file+'weights')
                if path_in_file+'z_events' in store.keys():
                    store.remove(path_in_file+'z_events')
                if path_in_file+'birth' in store.keys():
                    store.remove(path_in_file+'birth')
                    
        self._write_MODEL_data(self.filename, path_in_file, MODEL)
        
        rates = Rates(self.filename,
                      self.transient_name,
                      SFH_identifier,
                      verbose=self.verbose)
        
        z_birth = rates.centers_redshift_bins
        t_birth = rates.get_cosmic_time_from_redshift(z_birth)
        nr_of_birth_bins = len(z_birth)
        # write birth to the population file
        with pd.HDFStore(self.filename, mode='a') as store:
            store.put(path_in_file+'birth', pd.DataFrame(data={'z':z_birth,
                                                               't':t_birth}))
        
        get_redshift_from_time_cosmic_time = rates.redshift_from_cosmic_time_interpolator
        indices = self.indices

        # sample the SFH for only the events that are within the Hubble time
        # I only need to sample the SFH at each metallicity and z_birth
        # Not for every event!
        SFR_at_z_birth = star_formation_rate(rates.MODEL['SFR'], z_birth)
        # get metallicity bin edges
        met_edges = rates.edges_metallicity_bins

        # get the fractional SFR at each metallicity and z_birth
        fSFR = SFR_Z_fraction_at_given_redshift(z_birth,
                                                rates.MODEL['SFR'],
                                                rates.MODEL['sigma_SFR'],
                                                met_edges,
                                                rates.MODEL['Z_max'],
                                                rates.MODEL['select_one_met'])
        
        # simulated mass per given metallicity corrected for the unmodeled
        # single and binary stellar mass
        M_model = rates.mass_per_met.loc[rates.centers_metallicity_bins/Zsun]['underlying_mass'].values
        
        # speed of light
        c = const.c.to('Mpc/yr').value  # Mpc/yr
        
        # delta cosmic time bin
        deltaT = rates.MODEL['delta_t'] * 10 ** 6  # yr
    
        for i in tqdm(range(0, len(indices), self.chunksize), desc='event loop', disable=not self.verbose):
        
            selected_indices = self.select(start=i,
                                           stop=i+self.chunksize,
                                           columns=['index']).index.to_numpy().flatten()

            #selected_indices = indices[i:i+self.chunksize]
            delay_time = self.select(start=i,
                                     stop=i+self.chunksize,
                                     columns=['time']).to_numpy() *1e-3 # Gyr
            
            t_events = t_birth + delay_time
            hubble_time_mask = t_events  <= cosmology.age(1e-08).value*0.9999999
            
            # get the redshift of the events
            z_events = np.full(t_events.shape, np.nan)
            z_events[hubble_time_mask] = get_redshift_from_time_cosmic_time(t_events[hubble_time_mask])
            
            
            D_c = rates.get_comoving_distance_from_redshift(z_events)  # Mpc
            
            # the events have to be in solar metallicity
            met_events = self.select(start=i,
                                      stop=i+self.chunksize,
                                      columns=['metallicity']).to_numpy().flatten() * Zsun

            weights = np.zeros((len(met_events), nr_of_birth_bins))
            for i, met in enumerate(rates.centers_metallicity_bins):
                mask = met_events == met
                weights[mask,:] = 4.*np.pi * c * D_c[mask]**2 * deltaT * (fSFR[:,i]*SFR_at_z_birth) / M_model[i] # yr^-1
            
            
            with pd.HDFStore(self.filename, mode='a') as store:
                store.append(path_in_file+'weights',
                             pd.DataFrame(data=weights, index=selected_indices),
                             format='table')
                store.append(path_in_file+'z_events',
                             pd.DataFrame(data=z_events, index=selected_indices),
                             format='table')
        return rates
    
    def plot_efficiency_over_metallicity(self, **kwargs):
        '''plot the efficiency over metallicity
        
        Parameters
        ----------
        channel : bool
            plot the subchannels
        '''
        if not hasattr(self, 'efficiency'):
            raise ValueError('First you need to compute the efficiency over metallicity!')
        plot_pop.plot_merger_efficiency(self.efficiency.index.to_numpy()*Zsun, self.efficiency, **kwargs)

    def plot_delay_time_distribution(self, metallicity=None, ax=None, bins=100,):
        '''Plot the delay time distribution of the transient population'''
        
        if ax is None:
            fig, ax = plt.subplots()
            
        if metallicity is None:
            time = self.select(columns=['time']).values
            time = time*1e6 # yr
            h, bin_edges = np.histogram(time, bins=bins)
            h = h/np.diff(bin_edges)/self.mass_per_met['underlying_mass'].sum()
            
        else:
            if not any(np.isclose(metallicity, self.solar_metallicities)):
                raise ValueError('The metallicity is not present in the population!')
        
            time = self.select(where='metallicity == {}'.format(metallicity), columns=['time']).values
            time = time*1e6 # yr
            h, bin_edges = np.histogram(time, bins=bins)
            h = h/np.diff(bin_edges)/self.mass_per_met['underlying_mass'][metallicity]
                    
        ax.step(bin_edges[:-1], h, where='post', color='black')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('Time [yr]')
        ax.set_ylabel('Number of events/Msun/yr')
    
    def plot_popsyn_over_grid_slice(self, grid_type, met_Zsun, **kwargs):
        '''Plot the transients over the grid slice'''
        
        plot_pop.plot_popsyn_over_grid_slice(pop=self,
                                             grid_type=grid_type, 
                                             met_Zsun=met_Zsun,
                                             **kwargs)
    
    def _write_MODEL_data(self, filename, path_in_file, MODEL):
        with pd.HDFStore(filename, mode='a') as store:
            if MODEL['dlogZ'] is not None:
                store.put(path_in_file+'MODEL', pd.DataFrame(MODEL))
            else:
                store.put(path_in_file+'MODEL', pd.DataFrame(MODEL, index=[0]))
            if self.verbose:
                print('MODEL written to population file!')
    
class Rates(TransientPopulation):
    
    def __init__(self, filename, transient_name, SFH_identifier, verbose=False):
        
        super().__init__(filename, transient_name, verbose=verbose)
        self.SFH_identifier = SFH_identifier
        
        self.base_path = '/transients/'+self.transient_name+'/rates/'+self.SFH_identifier+'/'
        
        with pd.HDFStore(self.filename, mode='r') as store:
            if '/transients/'+self.transient_name+'/rates/'+self.SFH_identifier+'/MODEL' not in store.keys():
                raise ValueError(f'{self.SFH_identifier} is not a valid SFH_identifier in {filename}!')
            
        # load in the SFH_model
        self._read_MODEL_data(self.filename)
        
    def _read_MODEL_data(self, filename):
        with pd.HDFStore(filename, mode='r') as store:
            tmp_df = store[self.base_path+'MODEL']
            if len(tmp_df) > 1:
                self.MODEL = tmp_df.iloc[0].to_dict()
                self.MODEL['dlogZ'] = [tmp_df['dlogZ'].min(), tmp_df['dlogZ'].max()]
            else:
                self.MODEL = tmp_df.iloc[0].to_dict()    
            
            if self.verbose:
                print('MODEL read from population file!')    
    
    @property
    def weights(self):
        with pd.HDFStore(self.filename, mode='r') as store:
            return store[self.base_path+'weights']
           
    @property
    def z_birth(self):
        '''Get the redshift of the birth bins. Should return the same as centers_redshift_bins'''
        with pd.HDFStore(self.filename, mode='r') as store:
            return store[self.base_path+'birth']
        
    @property
    def z_events(self):
        with pd.HDFStore(self.filename, mode='r') as store:
            return store[self.base_path+'z_events']
    
    def select_rate_slice(self, key, start=None, stop=None):
        '''select a slice of the rates'''
        if key not in ['weights', 'z_events', 'birth']:
            raise ValueError('key not in [weights, z_events, birth]')
        
        with pd.HDFStore(self.filename, mode='r') as store:
            return store.select(self.base_path+key, start=start, stop=stop)
        
    
    def calculate_intrinsic_rate_density(self, mt_channels=False):
        '''Compute the intrinsic rate density of the transient population'''
        
        
        z_events = self.z_events.to_numpy()
        weights = self.weights.to_numpy()
        z_horizon = self.edges_redshift_bins
        n = len(z_horizon)
        
        if mt_channels:
            channels = self.select(columns=['channel'])
            unique_channels = np.unique(channels)
        else:
            unique_channels = []
            
        intrinsic_rate_density = pd.DataFrame(index=z_horizon[:-1], columns=['total'])
        
        normalisation = np.zeros(n-1)
        
        for i in tqdm(range(1, n), total=n-1, disable=not self.verbose):
            normalisation[i-1] = self.get_shell_comoving_volume(z_horizon[i-1],z_horizon[i],'infinite')
        
        for i in tqdm(range(1, n), total=n-1, disable=not self.verbose):
            mask = (z_events > z_horizon[i-1]) & (z_events <= z_horizon[i])
            for ch in unique_channels:
                mask_ch = channels.to_numpy() == ch
                intrinsic_rate_density.loc[z_horizon[i-1], ch] = np.nansum(weights[mask & mask_ch])/normalisation[i-1]
            
            intrinsic_rate_density.loc[z_horizon[i-1], 'total'] = np.nansum(weights[mask])/normalisation[i-1]
                
        with pd.HDFStore(self.filename, mode='a') as store:
            store.put(self.base_path+'intrinsic_rate_density', intrinsic_rate_density)
        
        return intrinsic_rate_density                               

    def calculate_observable_population(self, observable_func, observable_name):
        '''Calculate the observable population'''
        
        # file management and creation happens here
        
        # 1. recalculate the weights based on the observability function.
        # 2. the observability function takes the TransientPopulation as input
        
        # 3. We loop over the transient population (events) and recalulate the weights.
        # 4. the observability function takes the a transient population, the weights, and the z_events as input.
        # It outputs a new weights as an output
        # 5. It does this chunkwise
        
        with pd.HDFStore(self.filename, mode='a') as store:
            # remove the observable population if it already exists
            if '/transients/'+self.transient_name+'/rates/observable/'+observable_name in store.keys():
                if self.verbose:
                    print('Overwriting observable population!')
                del store['transients/'+self.transient_name+'/rates/observable/'+observable_name]
            
            
        # loop over the transient population and calculate the new weights, while writing to the file
        for i in tqdm(range(0, len(self), self.chunksize), total=len(self)//self.chunksize, disable=not self.verbose):
            transient_pop_chunk = self.select(start=i, stop=i+self.chunksize)
            weights_chunk = self.select_rate_slice('weights', start=i, stop=i+self.chunksize)
            z_events_chunk = self.select_rate_slice('z_events', start=i, stop=i+self.chunksize)
            new_weights = observable_func(transient_pop_chunk, z_events_chunk, weights_chunk)
            
            with pd.HDFStore(self.filename, mode='a') as store:
                store.append('transients/'+self.transient_name+'/rates/observable/'+observable_name,
                            new_weights,
                            format='table')
                
    def observable_population(self, observable_name):
        '''Return the observable population'''
        with pd.HDFStore(self.filename, mode='r') as store:
            if '/transients/'+self.transient_name+'/rates/observable/'+observable_name not in store.keys():
                raise ValueError(f'{observable_name} is not a valid observable population!')
            else:
                return store['transients/'+self.transient_name+'/rates/observable/'+observable_name]

    @property
    def observable_population_names(self):
        '''Return the names of the observable populations'''
        with pd.HDFStore(self.filename, mode='r') as store:
            return [key.split('/')[-1] for key in store.keys() if '/transients/'+self.transient_name+'/rates/observable/' in key]

    @property        
    def intrinsic_rate_density(self):
        with pd.HDFStore(self.filename, mode='r') as store:
            if self.base_path+'intrinsic_rate_density' not in store.keys():
                raise ValueError('First you need to compute the intrinsic rate density!')
            else:
                return store[self.base_path+'intrinsic_rate_density']
       
    def plot_hist_properties(self, prop, intrinsic=True, observable=None, bins=50, channel=None, **kwargs):
        '''plot a histogram of a given property available in the transient population'''
        
        if prop not in self.columns:
            raise ValueError(f'{prop} is not a valid property in the transient population!')
       
        # get the property and its associated weights in the population.
        
        df = self.select(columns=[prop])
        df['property'] = df[prop]
        del df[prop]
        if intrinsic:
            df['intrinsic'] = np.sum(self.weights, axis=1)       

        if observable is not None:
            with pd.HDFStore(self.filename, mode='r') as store:
                if '/transients/'+self.transient_name+'/rates/observable/'+observable not in store.keys():
                    raise ValueError(f'{observable} is not a valid observable population!')
                else:
                    df['observable'] = np.sum(store['transients/'+self.transient_name+'/rates/observable/'+observable], axis=1)
        if channel is not None:
            df['channel'] = self.select(columns=['channel'])
            df = df[df['channel'] == channel]
            if len(df) == 0:
                raise ValueError(f'{channel} is not present in the transient population!')
            plot_pop.plot_hist_properties(df, bins=bins, **kwargs)
            
        else:
            # plot the histogram using plot_pop.plot_hist_properties
            plot_pop.plot_hist_properties(df, bins=bins, **kwargs)
        
        
        
    
                                                    
    #### cosmolgy ####
    ##################
    
    def get_shell_comoving_volume(self, z_hor_i, z_hor_f, sensitivity='infinite'):
        """Compute comoving volume corresponding to a redshift shell.

        Parameters
        ----------
        z_hor_i : double
            Cosmological redshift. Lower bound of the integration.
        z_hor_f : double
            Cosmological redshift. Upper bound of the integration.
        sensitivity : string
            hoose which GW detector sensitivity you want to use. At the moment
            only 'infinite' is available, i.e. p_det = 1.

        Returns
        -------
        double
            Retruns the comoving volume between the two shells z_hor_i
            and z_hor_f in Gpc^3.

        """
        c = const.c.to('Gpc/yr').value  # Gpc/yr
        H_0 = cosmology.H(0).to('1/yr').value # km/Gpc*s
        def E(z):
            Omega_m = cosmology.Om0
            Omega_L = 1-cosmology.Om0
            return np.sqrt(Omega_m*(1.+z)**3+Omega_L)
        def f(z,sensitivity):
            if sensitivity=='infinite':
                return (1./(1.+z) * 4*np.pi*c / H_0
                        * (self.get_comoving_distance_from_redshift(z)
                        * 10**(-3.))**2. / E(z))
            else:
                # TODO: peanut-shaped antenna patter comoving volume calculation
                raise ValueError('Sensitivity not supported!')
        return sp.integrate.quad(f, z_hor_i, z_hor_f, args=(sensitivity))[0] # Gpc^3
        
    def get_comoving_distance_from_redshift(self, z):
        """Compute the comoving distance from redshift.

        Parameters
        ----------
        z : double
            Cosmological redshift.

        Returns
        -------
        double
            Comoving distance in Mpc corresponding to the redhisft z.

        """
        return cosmology.comoving_distance(z).value  # Mpc    
    
    ### metallicity bins ###
    ########################
    @property
    def edges_metallicity_bins(self):
        """Return the edges of the metallicity bins.

        Returns
        -------
        array double
            Returns the edges of all metallicity bins. We assume metallicities
            were binned in log-space.

        """
        met_val = np.log10(self.centers_metallicity_bins)
        bin_met = np.zeros(len(met_val)+1)
        # if more than one metallicty bin
        if len(met_val) > 1 :
            bin_met[0] = met_val[0] - (met_val[1] - met_val[0]) / 2.
            bin_met[-1] = met_val[-1] + (met_val[-1] - met_val[-2]) / 2.
            bin_met[1:-1] = met_val[:-1] + (met_val[1:] - met_val[:-1]) / 2.
        # one metallicty bin
        elif len(met_val) == 1 :
            if isinstance(self.MODEL['dlogZ'], float):
                bin_met[0] = met_val[0] - self.MODEL['dlogZ'] / 2.
                bin_met[-1] = met_val[0] + self.MODEL['dlogZ'] / 2.
            elif isinstance(self.MODEL['dlogZ'], list) or isinstance(self.MODEL['dlogZ'], np.array):
                bin_met[0] = self.MODEL['dlogZ'][0]
                bin_met[-1] = self.MODEL['dlogZ'][1]
        
        return 10**bin_met
    
    @property
    def centers_metallicity_bins(self):
        """Return the centers of the metallicity bins.

        Returns
        -------
        array double
            Returns sampled metallicities of the population. This corresponds
            to the center of each metallicity bin.
        """
        return np.sort(self.metallicities)
    
    ### redshift bins ###
    #####################
    @property
    def edges_redshift_bins(self):
        """Compute redshift bin edges.

        Returns
        -------
        array doubles
            We devide the cosmic time history of the Universe in equally spaced
            bins of cosmic time of self.delta_t (100 Myr default) an compute the
            redshift corresponding to edges of these bins.

        """
        self.n_redshift_bin_centers = int(cosmology.age(0).to('Myr').value/self.MODEL['delta_t'])
        # generate t_birth at the middle of each self.delta_t bin
        t_birth_bin = [cosmology.age(0.).value]
        for i in range(self.n_redshift_bin_centers+1):
            t_birth_bin.append(t_birth_bin[i] - self.MODEL['delta_t']*1e-3) # Gyr
        # compute the redshift
        z_birth_bin = []
        for i in range(self.n_redshift_bin_centers):
            # do not count first edge, we add z=0. later
            z_birth_bin.append(z_at_value(cosmology.age, t_birth_bin[i+1] * u.Gyr))
        # add the first and last bin edge at z=0. and z=inf.=100
        z_birth_bin = np.array([0.]+z_birth_bin+[100.])

        return z_birth_bin
    
    @property
    def centers_redshift_bins(self):
        """Compute redshift bin centers.

        Returns
        -------
        array doubles
            We devide the cosmic time history of the Universe in equally spaced
            bins of cosmic time of self.delta_t (100 Myr default) an compute the
            redshift corresponding to center of these bins.

        """
        # generate t_birth at the middle of each self.delta_t bin
        t_birth_bin = [cosmology.age(0.).value]
        t_birth = []
        # devide the
        self.n_redshift_bin_centers = int(cosmology.age(0).to('Myr').value/self.MODEL['delta_t'])
        for i in range(self.n_redshift_bin_centers+1):
            t_birth.append(t_birth_bin[i] - self.MODEL['delta_t']*1e-3/2.) # Gyr
            t_birth_bin.append(t_birth_bin[i] - self.MODEL['delta_t']*1e-3) # Gyr
        t_birth = np.array(t_birth)
        # compute the redshift
        z_birth = []
        for i in range(self.n_redshift_bin_centers+1):
            # z_at_value is from astopy.cosmology
            z_birth.append(z_at_value(cosmology.age, t_birth[i] * u.Gyr))
        z_birth = np.array(z_birth)

        return z_birth
    
    ###############################
    ### cosmic time to redshift ###
    ###############################
    @property
    def redshift_from_cosmic_time_interpolator(self):
        """Interpolator to compute the cosmological redshift given the cosmic time.

        Returns
        -------
        object
            Returns the trained interpolator object.

        """
        # astropy z_at_value method is too slow to compute z_mergers efficinty
        # we must implement interpolation
        t = np.linspace(1e-2, cosmology.age(1e-08).value*0.9999999, 1000)
        z = np.zeros(1000)
        for i in range(1000):
            z[i] = z_at_value(cosmology.age, t[i] * u.Gyr)
        f_z_m = interp1d(t, z, kind='cubic')
        return f_z_m
    
    def get_redshift_from_cosmic_time(self, t_cosm):
        """Compute the cosmological redshift given the cosmic time..

        Parameters
        ----------
        t_cosm : array doubles
            Cosmic time to which you want to know the redhisft.

        Returns
        -------
        array doubles
            Cosmolgocial redshift corresponding to the cosmic time.

        """
        return self.redshift_from_cosmic_time_interpolator(t_cosm)
    
    ###############################
    ### redshift to cosmic time ###
    ###############################
    
    def get_cosmic_time_from_redshift(self, z):
        """Compute the cosmic time from redshift.

        Parameters
        ----------
        z : double
            Cosmological redshift.

        Returns
        -------
        double
            Return age of the cosmic time in Gyr given the redshift z.

        """
        return cosmology.age(z).value  # Gyr