"""Evolve multiple BinaryPopulations together.

e.g. with multiple metallicities
"""

__authors__ = [
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Kyle Akira Rocha <kylerocha2024@u.northwestern.edu>",
]

import warnings
import numpy as np
import pandas as pd
from posydon.utils.constants import Zsun
from posydon.popsyn.io import parse_inifile, binarypop_kwargs_from_ini
from posydon.popsyn.binarypopulation import BinaryPopulation
from posydon.utils.common_functions import convert_metallicity_to_string
from posydon.popsyn.normalized_pop_mass import initial_total_underlying_mass
from posydon.utils.common_functions import inspiral_timescale_from_orbital_period
from posydon.popsyn.rate_calculation import Rates
import posydon.visualization.plot_dco as plot_dco


class SyntheticPopulation:

    def __init__(self, path_to_ini, path_to_data=None, verbose=False, MODEL=None):
        """
        Parameters
        ----------

        path : str
            Path to the inifile to parse. You can supply a list in the
            metallicity parameter to evolve more than one population.

        """

        self.verbose = verbose
        self.MODEL = MODEL
        self.df = None
        self.df_oneline = None
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
               binary_event=None, invert_S1S2=False, chunksize=500000):
        """Sort binaries of interests given some properties.

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

        df_sel = pd.DataFrame()
        df_sel_oneline = pd.DataFrame()
        count = 0
        tmp = 0
        if self.verbose:
            print('Binary count with (S1_state, S2_state, binary_state, binary_event) equal')
            print(f'to ({S1_state}, {S2_state}, {binary_state}, {binary_event})')
            if invert_S1S2:
                print(f'and ({S2_state}, {S1_state}, {binary_state}, {binary_event})')
        for k, file in enumerate(self.path_to_data):
            df_sel_met = pd.DataFrame()
            sel_met = []
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
                    sel_met.extend(sel)
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

            # load, parse and store oneline dataframe
            # this dataframe is smaller, we can load it all at once
            df_sel_met_oneline = pd.read_hdf(file,  key='oneline')
            df_sel_met_oneline = df_sel_met_oneline.loc[sel_met]
            df_sel_met_oneline['metallicity'] = met
            df_sel_oneline = pd.concat([df_sel_oneline, df_sel_met_oneline])

            if self.verbose:
                print(f'in {file} are {count-tmp}')
            tmp = count

        if self.verbose:
            print('Total binaries found are', count)

        # save parsed population as synthetic population
        if self.df is not None:
            warnings.warn('Overwriting the df population!')
        self.df = df_sel
        self.df_oneline = df_sel_oneline

    def save_pop(self, path='./parsed_population.h5'):
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

    def get_dco_at_formation(self, S1_state, S2_state, oneline_cols=None):
        """Sort synthetic population, i.e. DCO at formation.

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

        """

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

        # add properties of the oneline dataframe
        if self.df_oneline is not None:
            # TODO: add kicks as well by default?
            save_cols = ['S1_spin_orbit_tilt', 'S2_spin_orbit_tilt']
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

    def get_dco_merger_efficiency(self):
        """Compute the DCO merger efficinty per Msun for each metallicities."""
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
        """Compute the DCO merger rate weights.

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
        self.rates = Rates(self.df_synthetic, **self.MODEL)

        # compute DCO merger rate density
        if not load_data:
            self.rates.compute_merger_rate_weights(sensitivity=sensitivity, flag_pdet=flag_pdet, path_to_dir=working_dir)
        index, z_formation, z_merger, w_ijk = self.rates.load_merger_rate_weights(sensitivity, path_to_dir=working_dir)

        return index, z_formation, z_merger, w_ijk

    def resample_synthetic_population(self, index, z_formation, z_merger, w_ijk, export_cols=None):
        """Resample synthetc population to obtain intrinsic/detectable population.

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
            List of additional columns to save in the underlying/detectable
            population.

        Returns
        -------
        pd.DataFrame
            Resampled synthetc population to intrinsic or detecatable
            population.

        """

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
            df[c] = self.rates.get_data(c, index)

        return df

    def get_dco_merger_rate_density(self, export_cols=None,  working_dir='./', load_data=False):
        """Compute the merger rate density as a function of redshift.

        Parameters
        ----------
        export_cols : list str
            List of additional columns to save in the underlying/detectable
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
        index, z_formation, z_merger, w_ijk = self.compute_cosmological_weights(sensitivity, flag_pdet, working_dir=working_dir, load_data=load_data)

        # compute rate density weights
        self.z_rate_density = self.rates.get_centers_redshift_bins()
        self.rate_density = self.rates.compute_merger_rate_density(w_ijk, z_merger, observable='DCOs', sensitivity=sensitivity)
        print(f'DCO merger rate density in the local Universe (z={self.z_rate_density[0]:1.2f}): {round(self.rate_density[0],2)} Gpc^-3 yr^-1')

        # export the intrinsic DCO population
        self.df_intrinsic = self.resample_synthetic_population(index, z_formation, z_merger, w_ijk, export_cols=export_cols)

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
            List of additional columns to save in the underlying/detectable
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

        # export the detectable DCO population
        # TODO: store p_det
        self.df_detectable = self.resample_synthetic_population(index, z_formation, z_merger, w_ijk, export_cols=export_cols)

    def save_intrinsic_pop(self, path='./intrinsic_population.h5'):
        """Save intrinsic population.

        Parameters
        ----------
        path : str
            Path to dataset.

        """
        if self.df_intrinsic is None:
            raise ValueError('Nothing to save!')
        else:
            self.df_intrinsic.to_hdf(path, key='history')
            if self.verbose:
                print('Intrinsic population successfully saved!')

    def load_intrinsic_pop(self, path):
        """Load intrinsic population.

        Parameters
        ----------
        path : str
            Path to dataset.

        """
        if self.df_intrinsic is None:
            self.df_intrinsic = pd.read_hdf(path, key='history')
            if self.verbose:
                print('Intrinsic population successfully loaded!')
        else:
            raise ValueError('You already have an intrinsic population stored in memory!')

    def save_detectable_pop(self, path='./detectable_population.h5'):
        """Save detectable population.

        Parameters
        ----------
        path : str
            Path to dataset.

        """
        if self.df_detectable is None:
            raise ValueError('Nothing to save!')
        else:
            self.df_detectable.to_hdf(path, key='history')
            if self.verbose:
                print('Detectable population successfully saved!')

    def load_detectable_pop(self, path):
        """Load detectable population.

        Parameters
        ----------
        path : str
            Path to dataset.

        """
        if self.df_detectable is None:
            self.df_detectable = pd.read_hdf(path, key='history')
            if self.verbose:
                print('Detectable population successfully loaded!')
        else:
            raise ValueError('You already have an detectable population stored in memory!')

    def plot_merger_efficiency(self, **kwargs):
        """Plot merger rate efficinty."""
        if self.met_merger_efficiency is None or self.merger_efficiency is None:
            raise ValueError('First you need to compute the merger efficinty!')
        plot_dco.plot_merger_efficiency(self.met_merger_efficiency, self.merger_efficiency, **kwargs)

    def plot_merger_rate_density(self, **kwargs):
        """Plot merger rate density."""
        if self.z_rate_density is None or self.rate_density is None:
            raise ValueError('First you need to compute the merger rate density!')
        plot_dco.plot_merger_rate_density(self.z_rate_density, self.rate_density, **kwargs)

    def plot_hist_dco_properties(self, var, intrinsic=False, detectable=False, **kwargs):
        """Plot histogram of intrinsic/detectable properites.

        Parameters
        ----------
        var : str
            Property to plot stored in intrinsic/detectable dataframe.
        intrinsic : bool
            `True` if you want to deplay the intrisc population.
        detectable : bool
            `True` if you want to deplay the detectable population.
        **kwargs : dict
            ploting arguments

        """
        if self.z_rate_density is None or self.rate_density is None:
            raise ValueError('First you need to compute the merger rate density!')
        if intrinsic:
            intrinsic = self.df_intrinsic
        if detectable:
            detectable = self.df_detectable
        plot_dco.plot_hist_dco_properties(var, df_intrinsic=intrinsic, df_detectable=detectable, **kwargs)
