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

    def __init__(self, path_to_ini, path_to_data=None, verbose=False):
        """
        Parameters
        ----------

        path : str
            Path to the inifile to parse. You can supply a list in the
            metallicity parameter to evolve more than one population.

        """

        self.verbose = verbose
        self.df = None
        self.df_synthetic = None

        if '.ini' not in path_to_ini:
            raise ValueError('You did not provide a valid path_to_ini!')
        else:
            synthetic_pop_params = binarypop_kwargs_from_ini(path_to_ini)

            self.metallicity = synthetic_pop_params['metallicity']

            if not isinstance( self.metallicity, list):
                self.metallicity = [self.metallicity]

            self.binary_populations = []
            for met in self.metallicity:
                self.ini_kw = binarypop_kwargs_from_ini(path_to_ini)
                self.ini_kw['metallicity'] = met
                self.ini_kw['temp_directory'] = self.create_met_prefix(met) + self.ini_kw['temp_directory']
                self.binary_populations.append(BinaryPopulation(**self.ini_kw))

        if path_to_data is None:
            return
        elif (isinstance(path_to_data, list) and '.h5' in path_to_data[0]) or ('.h5' in path_to_data):
            self.path_to_data = path_to_data

    def get_ini_kw(self):
        return self.ini_kw

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
                print(f'and ({S2_state}, {S1_state}, {binary_state}, {binary_event})')
        for k, file in enumerate(self.path_to_data):
            df_sel_met = pd.DataFrame()
            # read metallicity from path
            met = float(file.split('/')[-1].split('_Zsun')[0])*Zsun
            simulated_mass_for_met = 0.
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

            # store simulated and underlying stellar mass
            df_sel_met['simulated_mass_for_met'] = simulated_mass_for_met
            df_sel_met['underlying_mass_for_met'] = initial_total_underlying_mass(df=simulated_mass_for_met, **self.ini_kw)[0]

            # concatenate results
            df_sel = pd.concat([df_sel, df_sel_met])
            del df_sel_met

            if self.verbose:
                print(f'in {file} are {count-tmp}')
            tmp = count

        if self.verbose:
            print('Total binaries found are', count)

        # save parsed population as synthetic population
        if self.df is not None:
            warnings.warn('Overwriting the df population!')
        self.df = df_sel

    def save_pop(self, path='./parsed_population.h5'):
        if self.df is None:
            raise ValueError('Nothing to save! The population was not parsed.')
        else:
            self.df.to_hdf(path, key='history')
            if self.verbose:
                print('Population successfully saved!')

    def load_pop(self, path):
        if self.df is None:
            self.df = pd.read_hdf(path, key='history')
            if self.verbose:
                print('Population successfully loaded!')
        else:
            raise ValueError('You already have a population stored in memory!')

    def get_dco_at_formation(self, S1_state, S2_state):

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
        self.df_synthetic['t_inspiral'] = (time_contact - self.df_synthetic[['time']])*1e-6 # Myr
        self.df_synthetic['time'] *= 1e-6 # Myr

        # for convinience reindex the DataFrame
        n_rows = len(self.df_synthetic.index)
        self.df_synthetic = self.df_synthetic.set_index(np.linspace(0,n_rows-1,n_rows,dtype=int))
