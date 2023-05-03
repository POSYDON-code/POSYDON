"""Evolve multiple BinaryPopulations together.

e.g. with multiple metallicities
"""



__authors__ = [
    "Kyle Akira Rocha <kylerocha2024@u.northwestern.edu>",
    "Simone Bavera <Simone.Bavera@unige.ch>",
]

import warnings
import numpy as np
import pandas as pd
from posydon.popsyn.io import parse_inifile, binarypop_kwargs_from_ini
from posydon.popsyn.binarypopulation import BinaryPopulation
from posydon.utils.common_functions import convert_metallicity_to_string
from posydon.popsyn.normalized_pop_mass import initial_total_underlying_mass
from posydon.utils.common_functions import inspiral_timescale_from_orbital_period


class SyntheticPopulation:

    def __init__(self, path=None, verbose=False):
        """
        Parameters
        ----------

        path : str
            Path to the inifile to parse. You can supply a list in the
            metallicity parameter to evolve more than one population.

        """

        self.verbose = verbose
        self.df_synthetic = None

        if path is None:
            return
        elif '.ini' in path:
            synthetic_pop_params = binarypop_kwargs_from_ini(path)

            self.metallicity = synthetic_pop_params['metallicity']

            if not isinstance( self.metallicity, list):
                self.metallicity = [self.metallicity]

            self.binary_populations = []
            for met in self.metallicity:
                ini_kw = binarypop_kwargs_from_ini(path)
                ini_kw['metallicity'] = met
                ini_kw['temp_directory'] = self.create_met_prefix(met) + ini_kw['temp_directory']
                self.binary_populations.append(BinaryPopulation(**ini_kw))
        elif (isinstance(path, list) and '.h5' in path[0]) or ('.h5' in path):
            self.path_to_results = path

    def evolve(self):
        """Evolve population(s) at given Z(s)."""
        for ind, pop in enumerate( self.binary_populations ):
            print( f'Z={pop.kwargs["metallicity"]:.2e} Z_sun' )
            pop.evolve()
            met_prefix = f'{pop.kwargs["metallicity"]:.2e}_Zsun_'
            pop.save( met_prefix + 'population.h5' )

    @staticmethod
    def create_met_prefix(met):
        """Append a prefix to the name of directories for batch saving."""
        return convert_metallicity_to_string(met) + '_Zsun_'

    def apply_logic(self, df, S1_state=None, S2_state=None, binary_state=None,
                    binary_event=None, step_name=None, invert_S1S2=False):
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


    def parse(self, S1_state=None, S2_state=None, binary_state=None,
               binary_event=None, invert_S1S2=False, chunksize=500000,
               Zsun = 0.0142):

        # TODO: we should also export the oneline dataframe

        df_sel = pd.DataFrame()
        count = 0
        tmp = 0
        if self.verbose:
            print('Binary count with (S1_state, S2_state, binary_state, binary_event) equal')
            print(f'to ({S1_state}, {S2_state}, {binary_state}, {binary_event})')
            if invert_S1S2:
                print(f'and ({S1_state}, {S2_state}, {binary_state}, {binary_event})')
        for k, file in enumerate(self.path_to_results):
            # read metallicity from path
            met = float(file.split('/')[-1].split('_Zsun')[0])*Zsun
            for i, df in enumerate(pd.read_hdf(file,  key='history', chunksize=chunksize)):

                logic = self.apply_logic(df, S1_state=S1_state,
                                    S2_state=S2_state,
                                    binary_state=binary_state,
                                    binary_event=binary_event,
                                    invert_S1S2=invert_S1S2)

                # select systems
                # remove duplicate indicies, e.g. if selecting 'contact' state it appears twice
                # if no specific event is selected (the second time is from the copied END event)
                sel = df.loc[logic].index.drop_duplicates()

                # count systems
                count += len(np.unique(sel))

                # sort systems
                if any(sel):
                    df_tmp = pd.DataFrame()
                    df_tmp = df.loc[sel]
                    df_tmp['metallicity'] = met
                    df_sel = pd.concat([df_sel,df_tmp])

            if self.verbose:
                print(f'in {file} are {count-tmp}')
            tmp = count

        if self.verbose:
            print('Total binaries found are', count)

        # save parsed population as synthetic population
        if self.df_synthetic is not None:
            warnings.warn('Overwriting the df_synthetic population!')
        self.df_synthetic = df_sel

    def save_synthetic_pop(self, path='./parsed_population.h5'):
        if self.df_synthetic is None:
            raise ValueError('Nothing to save! The population was not parsed.')
        else:
            self.df_synthetic.to_hdf(path, key='history')
            if self.verbose:
                print('Population succesfully saved!')

    def load_synthetic_pop(self, path):
        if self.df_synthetic is None:
            self.df_synthetic = pd.read_hdf(path, key='history')
            if self.verbose:
                print('Population succesfully loaded!')
        else:
            raise ValueError('You already have a population stored in memory!')
