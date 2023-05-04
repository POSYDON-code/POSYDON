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
from posydon.popsyn.cosmology import DCOrates
import posydon.visualization.plot_dco as plot_dco


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
        self.df_intrinsic = None
        self.df_detectable = None
        self.met_merger_efficiency = None
        self.merger_efficiency = None
        self.z_rate_density = None
        self.rate_density = None

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
        if not invert_S1S2 and S1_state != S2_state:
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
            # TODO: handle binaries at the edge case of the chuncks
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
        self.df_synthetic['t_delay'] = (time_contact - self.df_synthetic[['time']])*1e-6 # Myr
        self.df_synthetic['time'] *= 1e-6 # Myr

        # for convinience reindex the DataFrame
        n_rows = len(self.df_synthetic.index)
        self.df_synthetic = self.df_synthetic.set_index(np.linspace(0,n_rows-1,n_rows,dtype=int))

    def save_synthetic_pop(self, path='./synthetic_population.h5'):
        if self.df_synthetic is None:
            raise ValueError('Nothing to save!')
        else:
            self.df_synthetic.to_hdf(path, key='history')
            if self.verbose:
                print('Synthetic population successfully saved!')

    def load_synthetic_pop(self, path):
        if self.df_synthetic is None:
            self.df_synthetic = pd.read_hdf(path, key='history')
            if self.verbose:
                print('Synthetic population successfully loaded!')
        else:
            raise ValueError('You already have an synthetic population stored in memory!')

    def get_dco_merger_efficiency(self):
        metallicities = np.unique(self.df_synthetic['metallicity'])
        efficiencies = []
        self.met_merger_efficiency = sorted(metallicities)[::-1]
        for met in self.met_merger_efficiency:
            sel = (self.df_synthetic['metallicity'] == met)
            count = self.df_synthetic[sel].shape[0]
            underlying_stellar_mass = self.df_synthetic.loc[sel,'underlying_mass_for_met'].values[0]
            eff = count/underlying_stellar_mass
            efficiencies.append(eff)
            print(f'DCO merger efficiency at Z={met:1.2E}: {eff:1.2E} Msun^-1')
        self.met_merger_efficiency = np.array(self.met_merger_efficiency)
        self.merger_efficiency = np.array(efficiencies)


    def compute_cosmological_weights(self, sensitivity, flag_pdet, working_dir, load_data):

        # TODO: make the class inputs kwargs
        self.cosmology = DCOrates(self.df_synthetic)

        # compute DCO merger rate density
        if not load_data:
            self.cosmology.RunDCOsSimulation(sensitivity=sensitivity, flag_pdet=flag_pdet, path_to_dir=working_dir)
        index, z_formation, z_merger, w_ijk = self.cosmology.loadDCOsSimulation(sensitivity, path_to_dir=working_dir)

        return index, z_formation, z_merger, w_ijk

    def resample_synthetic_population(self, index, z_formation, z_merger, w_ijk, export_cols=None):

        # drop all zero weights to save memory
        sel = w_ijk > 0.
        index = index[sel]
        z_formation = z_formation[sel]
        z_merger = z_merger[sel]
        w_ijk = w_ijk[sel]

        # export results, we do a selections of columns else the dataframe
        # risks to go out of memory
        df = pd.DataFrame()
        df['weight'] = w_ijk
        df['z_formation'] = z_formation
        df['z_merger'] = z_merger
        save_cols = ['metallicity','time','t_delay','S1_state','S2_state',
                     'S1_mass','S2_mass','S1_spin','S2_spin',
                     'orbital_period','eccentricity', 'q', 'm_tot',
                     'm_chirp', 'chi_eff']
        if export_cols is not None:
            for c in export_cols:
                if c not in save_cols:
                    save_cols.append(c)
        for c in save_cols:
            df[c] = self.cosmology.get_data(c, index)

        return df

    def get_dco_merger_rate_density(self, export_cols=None,  working_dir='./', load_data=False):
        if self.df_synthetic is None:
            raise ValueError('You first need to isolated the DCO synthetic population!')

        # compute cosmological weights (detection rate weights with infinite sensitivity)
        sensitivity='infinite'
        flag_pdet = False
        index, z_formation, z_merger, w_ijk = self.compute_cosmological_weights(sensitivity, flag_pdet, working_dir=working_dir, load_data=load_data)

        # compute rate density weights
        self.z_rate_density = self.cosmology.getRedshiftBinCenter()
        self.rate_density = self.cosmology.RateDensity(w_ijk, z_merger, Type='DCOs', sensitivity=sensitivity)
        print(f'DCO merger rate density in the local Universe ({self.z_rate_density[0]:1.2f}): {round(self.rate_density[0],2)} Gpc^-3 yr^-1')

        # export the intrinsic DCO population
        self.df_intrinsic = self.resample_synthetic_population(index, z_formation, z_merger, w_ijk, export_cols=export_cols)

    def get_dco_detection_rate(self, sensitivity='design', export_cols=None,  working_dir='./', load_data=False):
        if self.df_synthetic is None:
            raise ValueError('You first need to isolated the DCO synthetic population!')

        # compute detection rate weights
        flag_pdet = True
        index, z_formation, z_merger, w_ijk = self.compute_cosmological_weights(sensitivity, flag_pdet, working_dir=working_dir, load_data=load_data)
        print(f'DCO detection rate at {sensitivity} sensitivity: {sum(w_ijk):1.2f} yr^-1')

        # export the detectable DCO population
        # TODO: store p_det
        self.df_detectable = self.resample_synthetic_population(index, z_formation, z_merger, w_ijk, export_cols=export_cols)

    def save_intrinsic_pop(self, path='./intrinsic_population.h5'):
        if self.df_intrinsic is None:
            raise ValueError('Nothing to save!')
        else:
            self.df_intrinsic.to_hdf(path, key='history')
            if self.verbose:
                print('Intrinsic population successfully saved!')

    def load_intrinsic_pop(self, path):
        if self.df_intrinsic is None:
            self.df_intrinsic = pd.read_hdf(path, key='history')
            if self.verbose:
                print('Intrinsic population successfully loaded!')
        else:
            raise ValueError('You already have an intrinsic population stored in memory!')

    def save_detectable_pop(self, path='./detectable_population.h5'):
        if self.df_detectable is None:
            raise ValueError('Nothing to save!')
        else:
            self.df_detectable.to_hdf(path, key='history')
            if self.verbose:
                print('Detectable population successfully saved!')

    def load_detectable_pop(self, path):
        if self.df_detectable is None:
            self.df_detectable = pd.read_hdf(path, key='history')
            if self.verbose:
                print('Detectable population successfully loaded!')
        else:
            raise ValueError('You already have an detectable population stored in memory!')

    def plot_merger_efficiency(self, **kwargs):
        if self.met_merger_efficiency is None or self.merger_efficiency is None:
            raise ValueError('First you need to compute the merger efficinty!')
        plot_dco.plot_merger_efficiency(self.met_merger_efficiency, self.merger_efficiency, **kwargs)

    def plot_merger_rate_density(self, **kwargs):
        if self.z_rate_density is None or self.rate_density is None:
            raise ValueError('First you need to compute the merger rate density!')
        plot_dco.plot_merger_rate_density(self.z_rate_density, self.rate_density, **kwargs)

    def plot_hist_dco_properties(self, var, intrinsic=False, detectable=False, **kwargs):
        if self.z_rate_density is None or self.rate_density is None:
            raise ValueError('First you need to compute the merger rate density!')
        if intrinsic:
            intrinsic = self.df_intrinsic
        if detectable:
            detectable = self.df_detectable
        plot_dco.plot_hist_dco_properties(var, df_intrinsic=intrinsic, df_detectable=detectable, **kwargs)
