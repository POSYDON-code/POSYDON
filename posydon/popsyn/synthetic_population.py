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
#from posydon.popsyn.rate_calculation import Rates
import posydon.visualization.plot_pop as plot_pop
from posydon.popsyn.GRB import get_GRB_properties, GRB_PROPERTIES

# TODO: temp import, remove after TF2 classification is implemented in pop synth
from posydon.interpolation.IF_interpolation import IFInterpolator
from posydon.binary_evol.binarystar import BinaryStar
from posydon.binary_evol.singlestar import SingleStar
from posydon.utils.common_functions import convert_metallicity_to_string


from astropy.cosmology import Planck15 as cosmology
from astropy.cosmology import z_at_value
from scipy.interpolate import interp1d
from astropy import units as u
from astropy import constants as const
import scipy as sp


from posydon.popsyn.star_formation_history import star_formation_rate,fractional_SFR_at_given_redshift, SFR_Z_fraction_at_given_redshift


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
        """Append a prefix to the se of directories for batch saving."""
        return convert_metallicity_to_string(met) + '_Zsun_'
    
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


DEFAULT_MODEL = {
    'delta_t' : 100, # Myr
    'SFR' : 'IllustrisTNG',
    'sigma_SFR' : None,
    'Z_max' : 1.,
    'select_one_met' : False,
    'dlogZ' : None, # e.g, [np.log10(0.0142/2),np.log10(0.0142*2)]
    'Zsun' : Zsun,
    'compute_GRB_properties' : False,
    'GRB_beaming' : 1., # e.g., 0.5, 'Goldstein+15'
    'GRB_efficiency' : 0., # e.g., 0.01
    'E_GRB_iso_min' : 0., # e.g., 1e51 erg 
}

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
            print(previous, end)
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
        
            #store.put(path_in_file+'MODEL', pd.DataFrame(MODEL))
        
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
            store.put(path_in_file+'MODEL', pd.DataFrame(MODEL))
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

    @property        
    def intrinsic_rate_density(self):
        with pd.HDFStore(self.filename, mode='r') as store:
            return store[self.base_path+'intrinsic_rate_density']
        
                                                    
        
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
        
    def get_cosmic_time_at_event(self, z_birth, event_time):
        """Get cosmic time at DCO merger.

        Parameters
        ----------
        z_birth : double
            Cosmological redshift of formation of the DCO system
            (must be the same for every binary).
        event_time : double 
            Time at which the DCO system merges in Myr.

        Returns
        -------
        double
            Cosmic time in Gyr at DCO merger for all binaries born at z_birth.

        """
        t_birth = self.get_cosmic_time_from_redshift(z_birth) # Gyr
        # add the cosmic time at birth
        return t_birth + (event_time * 10 ** (-3))  # Gyr
    
    def observable_population(self):
        pass
    
    def resampled_population(self):
        pass
        



class PopulationOld():
    
    def __init__(self, population_files, ini_file=None, verbose=False, chunksize=500000):

        # initialise the population attributes
        self.verbose = verbose
        self.chunksize = chunksize
        
        # if the user provided a single string instead of a list of strings
        if isinstance(population_files, str):
            population_files = [population_files]
        
        self.population_files = population_files
        self.populations = []
        
        if ini_file is not None:
            # read the metallicities from the ini file and match them with the population files
            self.ini_params = binarypop_kwargs_from_ini(ini_file)
            self.solar_metallicities = np.array(self.ini_params['metallicity'])
            print(self.solar_metallicities)
            self.metallicities = self.solar_metallicities * Zsun
            if len(self.solar_metallicities) != len(self.population_files):
                raise ValueError('The number of metallicities in the ini file does not match the number of population files!')
            
            for file, met in zip(self.population_files, self.solar_metallicities):
                if not os.path.isfile(file):
                    raise ValueError(f'{file} does not exist!')
                print('Loading population file: ', file)
                print('attributing Metallicity: ', met)
                self.populations.append(PopulationFile(file, metallicity=met, ini_file=ini_file, verbose=self.verbose, chunksize=self.chunksize))
            
        else:
            # load the population files
            for file in self.population_files:
                if not os.path.isfile(file):
                    raise ValueError(f'{file} does not exist!')
                self.populations.append(PopulationFile(file, verbose=self.verbose, chunksize=self.chunksize))
        
        self.number_of_systems = np.array([pop.number_of_systems for pop in self.populations])
        self.total_systems = np.sum(self.number_of_systems)
        tmp_df = pd.concat([pop.mass_per_met for pop in self.populations])
        self.mass_per_met = tmp_df.groupby(tmp_df.index).sum()
        self.solar_metallicities = np.array([pop.metallicities for pop in self.populations]).flatten()
        self.metallicities = self.solar_metallicities * Zsun
        

        self.indices = []
        prev = 0
        for pop in self.populations:
            self.indices.extend(pop.indices+prev)
            prev += pop.number_of_systems
        self._cumsum = np.cumsum(self.number_of_systems) - self.number_of_systems[0]

    def history(self, key):
            
        if isinstance(key, int):
            if key > self.total_systems:
                raise ValueError('Index out of range!')
            # first find in which popfile the key is
            
            pop_idx = np.searchsorted(self._cumsum, key, side='right') 
            id = int(key-self._cumsum[pop_idx-1])
            df = self.populations[pop_idx-1].history[id]
            df.index = [key] * len(df)
            return df
        elif isinstance(key, str):
            df = pd.DataFrame()
            for i, pop in enumerate(self.populations):
                tmp_df = pop.history[key]
                tmp_df.rename(index={x:y for x,y in zip(pop.indices, np.arange(self._cumsum[i-1], self._cumsum[i]))}, inplace=True)
                df = pd.concat([df, tmp_df])
            return df
        
        elif (isinstance(key, list) or isinstance(key, np.ndarray)) and all(isinstance(x, int) for x in key):
            key = np.sort(key)
            pop_idx = np.searchsorted(self._cumsum, key, side='right')
            df = pd.DataFrame()
            
            for i in np.unique(pop_idx):
                idx = pop_idx == i
                index = key[idx] - self._cumsum[i-1]
                tmp_df = self.populations[i-1].history[index]
                tmp_df.rename(index={x:y for x,y in zip(index, key[idx])}, inplace=True)
                df = pd.concat([df, tmp_df])
            return df
        
            
        elif isinstance(key, np.ndarray) and key.dtype == bool:
            # split at the boundaries of the populations
            # length should be the total number of systems
            df = pd.DataFrame()
            for i in range(0,len(self.populations)):
                tmp_df = self.populations[0].history[key[self._cumsum[i-1]:self._cumsum[i]]]
                tmp_df.rename(index={x:y for x,y in zip(tmp_df.index, np.arange(self._cumsum[i-1], self._cumsum[i]))}, inplace=True)
                df = pd.concat([df, tmp_df])
            return df
        elif isinstance(key, slice):
            raise ValueError('Slicing is not supported!')
        
        else:
            raise ValueError('Invalid key type!')
        
    def oneline(self, key):
        '''access the oneline from multiple files'''
        if isinstance(key, int):
            if key > self.total_systems:
                raise ValueError('Index out of range!')
            
            pop_idx = np.searchsorted(self._cumsum, key, side='right')
            id = int(key-self._cumsum[pop_idx-1])
            df = self.populations[pop_idx-1].oneline[id]
            df.index = key
            return df
            

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
    
    def save_filename(self, parsed_pop_instance):
        with pd.HDFStore(self.parsed_pop_file, mode='a') as store:
            tmp_df = pd.Series(parsed_pop_instance.filename)
            store.put('filename', tmp_df)
        
    def load_per_metallicity_info(self, parsed_pop_instance):
        """Load the underlying mass, simulated mass, selected systems, 
        and total_systems"""
        with pd.HDFStore(self.parsed_pop_file, mode='r') as store:
            tmp_df = store['mass_per_met']
            parsed_pop_instance.underlying_mass_per_met = tmp_df['underlying_mass'].to_numpy()
            parsed_pop_instance.simulated_mass_per_met = tmp_df['simulated_mass'].to_numpy()
            parsed_pop_instance.selected_systems = tmp_df['selected_systems'].to_numpy()
            parsed_pop_instance.total_systems = tmp_df['total_systems'].to_numpy()
            
    def load_filename(self, parsed_pop_instance):
        """Load the filename array
        
        Contains which population files were used to create the
        Parsed Population file"""
        with pd.HDFStore(self.parsed_pop_file, model='r') as store:
            parsed_pop_instance.filename = store['filename'].to_numpy()

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
        
    # for the synthetic population only?
    def load_parsed_pop_path(self, synthetic_pop_instance):
        '''Load the path to the parsed population file'''
        with pd.HDFStore(self.parsed_pop_file, model='r') as store:
            synthetic_pop_instance.parsed_pop_path = store['parsed_pop_path'].to_numpy()[0]
    
    # for the synthetic population only?
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
        
        filename : str or list of str
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
            self.io.load_filename(self)
        
    def parse(self, filename, S1_state=None, S2_state=None, binary_state=None,
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
        if type(filename) is str and ('.h5' in filename):
            filename = [filename]
        
        # The ini file and filename metallities should match.
        if len(filename) != len(self.solar_metallicities):
            raise ValueError('The number of metallicities and data files do not match!')        
        
        # catch the case where the user did not provide a path to data 
        if (isinstance(filename, list)):
            for path in filename:
                if os.path.splitext(path)[-1] != '.h5':
                    raise ValueError('You did not provide a valid filename!')
        else:
            raise ValueError('You did not provide a valid filename!')
        
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

        self.filename = filename
        self.underlying_mass_per_met = np.zeros(len(self.metallicities))
        self.simulated_mass_per_met = np.zeros(len(self.metallicities))
        self.selected_systems = np.zeros(len(self.metallicities))
        self.total_systems = np.zeros(len(self.metallicities))
        
        if self.verbose:
            print('Binary count with (S1_state, S2_state, binary_state, binary_event, step_name) equal')
            print(f'to ({S1_state}, {S2_state}, {binary_state}, {binary_event}, {step_name})')
            if invert_S1S2:
                print(f'and ({S2_state}, {S1_state}, {binary_state}, {binary_event}, {step_name})')
        
        if ((S1_state == None) & (S2_state == None) & (binary_state == None) & (binary_event == None) & (step_name == None)) and self.verbose:
            print('You did not specify any conditions to select binaries!')
            print('All binaries will be selected!')
            
            for k, file in enumerate(filename):
                if self.verbose:
                    print(f'Parsing {file}...')
                    
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
                simulated_mass_for_met = 0.

                with pd.HDFStore(file, mode='r') as store:
                    hist_total = store.get_storer('history').nrows

                for df in tqdm(pd.read_hdf(file, key='history', chunksize=self.chunksize), total=hist_total/self.chunksize , disable=not self.verbose, desc='history'):
                    df.index += index_shift
                    simulated_mass_for_met += sum(df['S1_mass'][df['event'] == 'ZAMS']) + sum(df['S2_mass'][df['event'] == 'ZAMS'])
                    
                    # write history to file
                    parsed_population_file.append('history',
                                                df,
                                                data_columns=True,
                                                min_itemsize=history_min_itemsize)
                    
                
                with pd.HDFStore(file, mode='r') as store:
                    total = store.get_storer('oneline').nrows
                    
                for df in tqdm(pd.read_hdf(file, key='oneline', chunksize=self.chunksize), total=total/self.chunksize, disable=not self.verbose, desc='oneline'):
                    df.index += index_shift
                    parsed_population_file.append('oneline',
                                                df,
                                                data_columns=True,
                                                min_itemsize=oneline_min_itemsize)
                    max_index = df.index.max()
    
                self.simulated_mass_per_met[met_index] = simulated_mass_for_met        
                self.underlying_mass_per_met[met_index] = initial_total_underlying_mass(df=simulated_mass_for_met, **self.ini_params)[0]
                self.total_systems[met_index] = total
                self.selected_systems[met_index] = total
                
                index_shift += np.max([max_index, total])
                    
        else:
            for k, file in enumerate(filename):
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
                
                # Loop over the binaries using a history length index.
                # Store that history index too!
                for i, df in enumerate(pd.read_hdf(file,  key='history', chunksize=self.chunksize)):
                    
                    df = pd.concat([last_binary_df, df])
                        
                    last_binary_df = df.loc[[df.index[-1]]]
                    df.drop(df.index[-1], inplace=True)
                    total += len(df.index.drop_duplicates())
                    max_index = df.index.max()
                    
                    # TODO: skip the logic step
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
        self.io.save_filename(self)
        
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
            
            merged.index.se='binary_index'
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
            DCO_synthetic_population.io.save_filename(self)
            DCO_synthetic_population.io.save_parsed_pop_path(self.parsed_pop_file)
            
            # load data into DCO class
            DCO_synthetic_population.io.load_ini_params(DCO_synthetic_population)
            DCO_synthetic_population.io.load_parse_kwargs(DCO_synthetic_population)
            DCO_synthetic_population.io.load_per_metallicity_info(DCO_synthetic_population)
            DCO_synthetic_population.io.load_filename(DCO_synthetic_population)
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

    @property
    def history_lengths(self):
        '''Return the length of the history for each binary'''
        
        if not hasattr(self, '_history_lengths'):
            with pd.HDFStore(self.parsed_pop_file, mode='r') as store:
                history_events = store.select_column('history', 'index')
                self._history_lengths = history_events.groupby(history_events).count()
                del history_events
                
        return self._history_lengths
    
    def select(self, key, start, stop, columns=None):
        
        with pd.HDFStore(self.parsed_pop_file, mode='r') as parsed_store:
            return parsed_store.select(key,
                                start=start,
                                stop=stop,
                                columns=columns)

    
    @property
    def df(self):
        with pd.HDFStore(self.parsed_pop_file, mode='r') as store:
            return store.select('history')
    
    @property
    def df_oneline(self):
        with pd.HDFStore(self.parsed_pop_file, mode='r') as store:
            return store.select('oneline')

    def create_synpop(self, func, output_file, oneline_cols=None, hist_cols=None):
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
        synth_pop = SyntheticPopulation(output_file, verbose=self.verbose) 
        # write population data to the new file!
        synth_pop.io.save_ini_params(self)
        synth_pop.io.save_parse_kwargs(self)
        synth_pop.io.save_per_metallicity_info(self)
        synth_pop.io.save_filename(self)
        synth_pop.io.save_parsed_pop_path(self.parsed_pop_file)
            
        # load data into DCO class
        synth_pop.io.load_ini_params(synth_pop)
        synth_pop.io.load_parse_kwargs(synth_pop)
        synth_pop.io.load_per_metallicity_info(synth_pop)
        synth_pop.io.load_filename(synth_pop)
        synth_pop.io.load_parsed_pop_path(synth_pop)
                
        if '/synthetic' in synth_pop.keys():
            print('overwriting synthetic population')
            synth_pop.remove('synthetic')
            
        min_itemsize = {'channel': 100,}
        if hist_cols is not None:
            if 'time' not in hist_cols:
                raise ValueError('The synthetic population requires a time column!')
        
            min_itemsize.update({key:val for key, val in 
                                HISTORY_MIN_ITEMSIZE.items() 
                                if key in hist_cols})
        
        if oneline_cols is not None:
            min_itemsize.update({key:val for key, val in
                                ONELINE_MIN_ITEMSIZE.items()
                            if key in oneline_cols})

        # setup a mapping to the size of each history colummn
        history_lengths = self.history_lengths
        unique_binary_indices = self.indices
        
        previous = 0
        for i in tqdm(range(0,len(unique_binary_indices), self.chunksize), disable=not self.verbose):
            end = previous + history_lengths[i:i+self.chunksize].sum()
            
            oneline_chunk = self.select('oneline',
                                        start=i,
                                        stop=i+self.chunksize,
                                        columns=oneline_cols)
            
            history_chunk = self.select('history',
                                                start=previous,
                                                stop=end, 
                                                columns=hist_cols)
            
            formation_channels_chunk = self.select('formation_channels', 
                                                    start = i,
                                                    stop=i+self.chunksize,
                                                    columns=['channel'])
            
            syn_df = func(history_chunk, oneline_chunk, formation_channels_chunk)
            
            # filter out the columns in min_itemsize that are not in the dataframe
            min_itemsize = {key:val for key, val in min_itemsize.items() if key in syn_df.columns}
                
            synth_pop.append('synthetic',
                                syn_df,
                                format='table',
                                data_columns=True,
                                min_itemsize=min_itemsize
                                )
            
            previous = end
            
        return synth_pop

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
            output_store = parsed_store
        else:
            GRB_synthetic_population = SyntheticPopulation(output_file,
                                                           verbose=self.verbose)
            
            # write population data to the new file!
            GRB_synthetic_population.io.save_ini_params(self)
            GRB_synthetic_population.io.save_parse_kwargs(self)
            GRB_synthetic_population.io.save_per_metallicity_info(self)
            GRB_synthetic_population.io.save_filename(self)
            GRB_synthetic_population.io.save_parsed_pop_path(self.parsed_pop_file)
            
            # load data into GRB_data
            GRB_synthetic_population.io.load_ini_params(GRB_synthetic_population)
            GRB_synthetic_population.io.load_parse_kwargs(GRB_synthetic_population)
            GRB_synthetic_population.io.load_per_metallicity_info(GRB_synthetic_population)
            GRB_synthetic_population.io.load_filename(GRB_synthetic_population)
            GRB_synthetic_population.io.load_parsed_pop_path(GRB_synthetic_population)
            
            # Open the file for writing the synthetic population
            parsed_store = pd.HDFStore(self.parsed_pop_file, mode='r')
            output_store = pd.HDFStore(GRB_synthetic_population.pop_file, mode='a')

        
        if '/synthetic' in output_store.keys():
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
                GRB_df_synthetic['preSN_'+c] = pre_sn[c].valuesS
                
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

class SyntheticPopulation2():
    
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
            if '/filename' in keys:
                self.io.load_filename(self)
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
    
    def select(self, key, start, stop, columns=None):
        with pd.HDFStore(self.pop_file, mode='r') as store:
            return store.select(key, start=start, stop=stop, columns=columns)

    def remove(self, key):
        with pd.HDFStore(self.pop_file, mode='a') as store:
            store.remove(key)
    
    def append(self, key, df, format='table', data_columns=True, min_itemsize=None):
        with pd.HDFStore(self.pop_file, mode='a') as store:
            store.append(key,
                         df,
                         format=format,
                         data_columns=data_columns,
                         min_itemsize=min_itemsize
                         )
    
    def keys(self):
        with pd.HDFStore(self.pop_file, mode='r') as store:
            return store.keys()

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
        synthetic_population.io.save_filename(self)
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
        # Create a 
        if self.parsed_pop_path is None:
            raise ValueError('No parsed population file found!')
        try:
            tmp_pop = ParsedPopulation(self.parsed_pop_path)
        except FileNotFoundError:
            pass
        
        class temp_pop:
            def __init__(self, hist, oneline):
                self.df = hist
                self.df_oneline = oneline
        
        with pd.HDFStore(tmp_pop.parsed_pop_file, mode='r') as store:
            indices = self.indices.tolist()
            df = store.select('history', where='index in indices')
            df_oneline = store.select('oneline', where='index in indices')
            df_formation = store.select('formation_channels', where='index in indices')
            df_oneline['channel'] = df_formation['channel']
        
        tmp_pop = temp_pop(df, df_oneline)
        
        plot_pop.plot_popsyn_over_grid_slice(tmp_pop, grid_type, met_Zsun, **kwargs)
    
    
    
    
