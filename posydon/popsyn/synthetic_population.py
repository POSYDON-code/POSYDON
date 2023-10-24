"""Evolve multiple BinaryPopulations together.

e.g. with multiple metallicities
"""

__authors__ = [
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Kyle Akira Rocha <kylerocha2024@u.northwestern.edu>",
    "Monica Gallegos-Garcia <monicagallegosgarcia2024@u.northwestern.edu>",
]

import warnings
import numpy as np
import pandas as pd
from tqdm import tqdm
import os
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

class SyntheticPopulation:

    def __init__(self, path_to_ini, verbose=False, MODEL={}):
        """
        Parameters
        ----------

        path : str
            Path to the inifile to parse. You can supply a list in the
            metallicity parameter to evolve more than one population.
        """
        self.synthetic_pop_params = None
        self.metallicity = None
        self.binary_populations = None

        self.verbose = verbose
        self.MODEL = MODEL
        # DCOs
        self.df = None
        self.df_oneline = None
        self.df_synthetic = None
        self.df_dco_intrinsic = None
        self.df_dco_observable = None
        self.met_merger_efficiency = None
        self.merger_efficiency = None
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
            self.metallicities = self.synthetic_pop_params['metallicity']
            if not isinstance( self.metallicity, list):
                self.metallicity = [self.metallicity]
            self.binary_populations = None
        
    def create_binary_populations(self):
        self.binary_populations = []
        ini_kw = self.synthetic_pop_params.copy()
        for met in self.metallicity[::-1]:
            ini_kw['metallicity'] = met
            ini_kw['temp_directory'] = self.create_met_prefix(met) + self.synthetic_pop_params['temp_directory']
            self.binary_populations.append(BinaryPopulation(**ini_kw))

    def get_ini_kw(self):
        return self.synthetic_pop_params

    def evolve(self):
        """Evolve population(s) at given Z(s)."""
        if self.binary_populations is None:
            self.create_binary_populations()
        while self.binary_populations:
            pop =  self.binary_populations.pop()
            print( f'Z={pop.kwargs["metallicity"]:.2e} Z_sun' )
            pop.evolve()
            #met_prefix = f'{pop.kwargs["metallicity"]:.2e}_Zsun_'
            #pop.save( met_prefix + 'population.h5' )
            del pop


    def merge_parallel_run(self, path_to_batches):
        """
        Merge the folder or list of folders into a single file per metallicity.

        Parameters
        ----------
        path_to_batches : str or list of str
            Path to the folder(s) containing the batch folders.
        """
        
        if isinstance(path_to_batches, str):
            path_to_batches = [path_to_batches]
        # check if path_to_batches is the same length as the number of metallicities
        if len(path_to_batches) != len(self.metallicity):
            raise ValueError('The number of metallicity and batch directories do not match!')

        for met, path_to_batch in zip(self.metallicity, path_to_batches):
            met_prefix = self.create_met_prefix(met)
            tmp_files = [f for f in os.listdir(path_to_batch) if os.isfile(os.join(path_to_batch, f))]
        

            BinaryPopulation(**self.ini_kw).combine_saved_files( met_prefix + 'population.h5', tmp_files)


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
        
        # catch the case where the user did not provide a path to data
        if not ((isinstance(path_to_data, list) and '.h5' in path_to_data[0]) or ('.h5' in path_to_data)):
            raise ValueError('You did not provide a valid path_to_data!')
        
        df_sel = pd.DataFrame()
        df_sel_oneline = pd.DataFrame()
        count = 0
        tmp = 0
        if self.verbose:
            print('Binary count with (S1_state, S2_state, binary_state, binary_event, step_name) equal')
            print(f'to ({S1_state}, {S2_state}, {binary_state}, {binary_event}, {step_name})')
            if invert_S1S2:
                print(f'and ({S2_state}, {S1_state}, {binary_state}, {binary_event}, {step_name})')
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
                                    step_name=step_name,
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

            if k > 0 and df_sel.shape[0] > 0:
                shift_index = max(np.unique(df_sel.index)) + 1 
                df_sel_met.index += shift_index

            # store simulated and underlying stellar mass
            df_sel_met['simulated_mass_for_met'] = simulated_mass_for_met
            df_sel_met['underlying_mass_for_met'] = initial_total_underlying_mass(df=simulated_mass_for_met, **self.synthetic_pop_params)[0] # This used to be init_kw

            # concatenate results
            df_sel = pd.concat([df_sel, df_sel_met])
            del df_sel_met

            # load, parse and store oneline dataframe
            # this dataframe is smaller, we can load it all at once
            df_sel_met_oneline = pd.read_hdf(file,  key='oneline')
            df_sel_met_oneline = df_sel_met_oneline.loc[sel_met]
            df_sel_met_oneline['metallicity'] = met

            if k > 0 and df_sel_oneline.shape[0] > 0:
                shift_index = max(np.unique(df_sel_oneline.index)) + 1
                df_sel_met_oneline.index += shift_index 

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

    def get_dco_at_formation(self, S1_state, S2_state, oneline_cols=None, formation_channels=False, mt_history=False):
        """Populates `df_synthetic` with DCOs at their formation.
        If `formation_channels` is `True` the `channel` column is added to the
        `df_synthetic` dataframe.

        if MODEL on class initialization is not None and
        "compute_GRB_properties" in MODEL. 
        If MODEL["compute_GRB_properties"] is `True` the following columns are
        in the MODEL:

            - ['GRB_efficiency'],
            - self.MODEL['GRB_beaming'],
            - self.MODEL['E_GRB_iso_min']

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
        compute_GRB_properties = (self.MODEL is not None and 
                                  "compute_GRB_properties" in self.MODEL and 
                                  self.MODEL["compute_GRB_properties"])

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
            if ('GRB_efficiency' not in self.MODEL or
                self.MODEL['GRB_efficiency'] is None):
                raise ValueError('Missing GRB_efficiency variable in the MODEL!')
            if ('GRB_beaming' not in self.MODEL or
                self.MODEL['GRB_beaming'] is None):
                raise ValueError('Missing GRB_beaming variable in the MODEL!')
            if ('E_GRB_iso_min' not in self.MODEL or
                self.MODEL['E_GRB_iso_min'] is None):
                raise ValueError('Missing GRB_beaming variable in the MODEL!')
            self.df_synthetic = get_GRB_properties(self.df_synthetic,
                                                   self.MODEL['GRB_efficiency'],
                                                   self.MODEL['GRB_beaming'],
                                                   self.MODEL['E_GRB_iso_min']
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

    def get_dco_merger_efficiency(self):
        """Compute the DCO merger efficinty per Msun for each metallicities."""
        metallicities = np.unique(self.df_synthetic['metallicity'])
        efficiencies = []
        self.met_merger_efficiency = sorted(metallicities)[::-1]
        for met in self.met_merger_efficiency:
            sel = (self.df_synthetic['metallicity'] == met)
            count = self.df_synthetic[sel].shape[0]
            underlying_stellar_mass = self.df_synthetic.loc[sel,'underlying_mss_for_met'].values[0]
            eff = count/underlying_stellar_mass
            efficiencies.append(eff)
            print(f'DCO merger efficiency at Z={met:1.2E}: {eff:1.2E} Msun^-1')
        self.met_merger_efficiency = np.array(self.met_merger_efficiency)
        self.merger_efficiency = {'total' : np.array(efficiencies)}
        # if the channel column is present compute the merger efficiency per channel
        if "channel" in self.df_synthetic:
            channels = np.unique(self.df_synthetic['channel'])
            for ch in channels:
                efficiencies = []
                for met in self.met_merger_efficiency:
                    sel = (self.df_synthetic['metallicity'] == met) & (self.df_synthetic['channel'] == ch)
                    count = self.df_synthetic[sel].shape[0]
                    if count > 0:
                        underlying_stellar_mass = self.df_synthetic.loc[sel,'underlying_mass_for_met'].values[0]
                        eff = count/underlying_stellar_mass
                    else:
                        eff = np.nan
                    efficiencies.append(eff)
                    # print(f'Z={met:1.2E} {ch}: {eff:1.2E} Msun^-1')
                self.merger_efficiency[ch] = np.array(efficiencies)


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
        self.rates = Rates(self.df_synthetic, **self.MODEL)

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
        compute_GRB_properties = (self.MODEL is not None and 
                                  "compute_GRB_properties" in self.MODEL and 
                                  self.MODEL["compute_GRB_properties"])

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

    def get_formation_channels(self, mt_history):
        """Get formation channel and add to df and df_oneline."""
        
        # loop through each binary
        unique_binary_index = np.unique(self.df.index)
        for index in unique_binary_index:

            # get event column and information from interpolated classes
            df_binary = self.df.loc[index,['event']].dropna()
            intep_cls = [key for key in self.df_oneline.keys() if 'interp_class' in key]
            df_binary_online = self.df_oneline.loc[index, intep_cls].dropna()
            event_array = df_binary['event'].values.tolist()
        
            # make interpolated class information consistent with event column
            HMS_HMS_event_dict = {'stable_MT':'oRLO1', 'no_MT':'None', 'unstable_MT':'oCE1/oDoubleCE1'}
            event_HMS_HMS = HMS_HMS_event_dict[df_binary_online['interp_class_HMS_HMS']]
            
            # for now, only append information for RLO1; unstable_MT information already exists
            if event_HMS_HMS == 'oRLO1':
                event_array.insert(1, event_HMS_HMS)
                formation_channel = "_".join(event_array)
            else:
                formation_channel = "_".join(event_array)
                        
            # TODO: drop the envent CO_contact
            # TODO: once we trust the redirection to the detached step
            #       drop also the redirect event
            
            # TODO: for debugging purposes we keep the full channel
            self.df_oneline.loc[index,'channel_debug'] = formation_channel
            # clean the redirect and CO_contact events
            formation_channel = formation_channel.replace('_redirect', '')
            formation_channel = formation_channel.replace('_CO_contact', '')
            self.df_oneline.loc[index,'channel'] = formation_channel
            
        if mt_history and 'mt_history_HMS_HMS' not in self.df_oneline:
            raise ValueError('mt_history_HMS_HMS not saved in the oneline dataframe!')
        else:   
            # split oRLO1 into oRLO1, oRLO1-contact and oRLO1-reverse
            sel = ((self.df_oneline['mt_history_HMS_HMS'] == 'Stable contact phase') & 
                    self.df_oneline['channel'].str.contains('oRLO1'))
            self.df_oneline.loc[sel, 'channel'] = self.df_oneline.loc[sel, 'channel'].apply(lambda x: x.replace('oRLO1', 'oRLO1-contact'))
            sel = ((self.df_oneline['mt_history_HMS_HMS'] == 'Stable reverse mass-transfer phase') & 
                    self.df_oneline['channel'].str.contains('oRLO1'))
            self.df_oneline.loc[sel, 'channel'] = self.df_oneline.loc[sel, 'channel'].apply(lambda x: x.replace('oRLO1', 'oRLO1-reverse'))

            # TODO: do the above split for unstable MT as well
            
    
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

    def plot_merger_efficiency(self, **kwargs):
        """Plot merger rate efficinty."""
        if self.met_merger_efficiency is None or self.merger_efficiency is None:
            raise ValueError('First you need to compute the merger efficinty!')
        plot_pop.plot_merger_efficiency(self.met_merger_efficiency, self.merger_efficiency, **kwargs)

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