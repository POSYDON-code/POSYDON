
__author__ = ['Simone Bavera <Simone.Bavera@unige.ch>']

import os
import warnings
import numpy as np
import pandas as pd
from tqdm import tqdm
import scipy as sp
from scipy.interpolate import interp1d
from astropy import units as u
from astropy import constants as const # TODO: replace this in favour of POSYDON constats module
from astropy.cosmology import Planck15 as cosmology # TODO: update to Planck18 once we update astropy module
from astropy.cosmology import z_at_value
import posydon.popsyn.selection_effects as selection_effects
from posydon.utils.constants import Zsun
from posydon.popsyn.star_formation_history import (star_formation_rate,
                                                   mean_metallicity,
                                                   std_log_metallicity_dist,
                                                   get_illustrisTNG_data,
                                                   fractional_SFR_at_given_redshift)
from posydon.utils.data_download import PATH_TO_POSYDON_DATA
from posydon.popsyn.GRB import GRB_PROPERTIES, get_GRB_properties
import posydon.visualization.plot_pop as plot_pop

PATH_TO_PDET_GRID = os.path.join(PATH_TO_POSYDON_DATA, 'selection_effects/pdet_grid.hdf5')


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

from astropy.cosmology import Planck15 as cosmology
from astropy.cosmology import z_at_value
from scipy.interpolate import interp1d
from tqdm import tqdm
import numpy as np
import os
from astropy import units as u
from posydon.popsyn.star_formation_history import (star_formation_rate,
                                                   mean_metallicity,
                                                   std_log_metallicity_dist,
                                                   get_illustrisTNG_data,
                                                   fractional_SFR_at_given_redshift)
from astropy import constants as const
import warnings
from posydon.utils.constants import Zsun
import scipy as sp
import posydon.popsyn.selection_effects as selection_effects


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


class Rates(object):
    """Object for calculating rates from a synthetic population and a model.
    
    """
    
    def __init__(self, synthetic_population, MODEL=None, rates_file=None, verbose=False, chunksize=1000000):
        '''Initialize the Rates object.
        
        Parameters
        ----------
        synthetic_population : SyntheticPopulation
            The synthetic population object.
        MODEL : dict
            a dictionary describing the model
            delta_t : float
                Cosmic time bin size used to bin distribute the synthetic binary
                population.
            SFR : string
                Star-formation-rate history you want to use:
                'Madau+Dickinson14', 'Madau+Fragos17', 'Neijssel+19', 'IllustrisTNG'
            sigma_SFR : string
                Standard deviation of the trucated log-normal metallicity distribution
                assumed for the star-formation rate history in the case you select
                'Madau+Dickinson14', 'Madau+Fragos17', 'Neijssel+19',
            Z_max : float
                The absolute maximum metallicity of the star formation rate history.
                NOTE: This should be used when assuming a truncated log-normal
                distribution of metallcities to not overpredict 1, also it can be
                used with the IllustrisTNG SFR to prevent unphysically high
                metallicities.
            select_one_met : bool
                If `True` your synthetic binary population should contain only one
                descrete metallicity.
            dlogZ : float or range
                Used when `select_one_met=True` to select the metallicity range
                around which you want to integrate the SFR history. If float value
                value is provided we assume a symmetric range dlogZ/2 around the
                provided metallcity else we consider the range provided
                [Z_min, Z_max] which could also be asymmetric.
        verbose : bool
            `True` if you want the print statements.
        '''
        # Synthetic Population class
        self.synthetic_population = synthetic_population
        self.verbose = verbose
        self.chunksize = chunksize
        
        if rates_file is None:
            # store rates in the same file as the synthetic population
            self.rates_file = self.synthetic_population.pop_file
        else:
            self.rates_file = rates_file
            with pd.HDFStore(self.rates_file, 'r') as f:
                if '/rates' in f.keys():
                    raise ValueError("The provided rates_file already contains rates!")
        
        # Set model to DEFAULT or provided MODEL parameters
        # Allows for partial model specification
        if MODEL is None:
            self.MODEL = DEFAULT_MODEL
        else:
            for key in MODEL:
                if key not in DEFAULT_MODEL:
                    raise ValueError(key + " is not a valid parameter name!")
            
            # write the DEFAULT_MODEL with updates parameters to self.MODEL.
            self.MODEL = DEFAULT_MODEL
            self.MODEL.update(MODEL)
    
        # TODO: ADD the GRB MODEL property check.
        #if self.MODEL['compute_GRB_properties']:
        #    columns = self.synthetic_population.head(0).columns.to_numpy()
    
    @property
    def centers_metallicity_bins(self):
        """Return the centers of the metallicity bins.

        Returns
        -------
        array double
            Returns sampled metallicities of the populattion. This correponds
            to the center of each metallicity bin.

        """
        return np.sort(self.synthetic_population.metallicities)
    
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
    def redshift_from_cosmic_time_interpolator(self):
        """Interpolator to compute the cosmological redshift given the cosmic time.

        Returns
        -------
        object
            Returns the trained interpolator object.

        """
        # astropy z_at_value method is too slow to compute z_mergers efficinty
        # we must implment interpolation
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
    
    def chi_eff(self, m_1, m_2, a_1, a_2, tilt_1, tilt_2):
        return (m_1*a_1*np.cos(tilt_1)+m_2*a_2*np.cos(tilt_2))/(m_1+m_2)

    def m_chirp(self, m_1, m_2):
        return (m_1*m_2)**(3./5)/(m_1+m_2)**(1./5)

    def mass_ratio(self, m_1, m_2):
        q = m_2/m_1
        q[q>1.] = 1./q[q>1.]
        return q
    
    def get_data(self, var, index=None):
        """Resample class elements.

        Parameters
        ----------
        var : string
            Key of self.df. This method support a list of custom definitions for
            DCO studies like q, m_tot, m_chirp, chi_eff.
        index : list of int
            List of indicies of the binaries of interest.

        Returns
        -------
        array
            Array containing the self.df[var][index].

        """
        def read_data(columns, index):
            '''read the data from the synthetic population file
            
            Parameters
            -----------
            columns : list of strings
                List of columns to read from the synthetic population file.
            index : list of int
                List of indicies of the binaries of interest.
            
            Returns
            -------
            pandas.DataFrame
                DataFrame containing the data of interest.
            '''
            if index is not None:
                
                return pd.read_hdf(self.synthetic_population.pop_file,
                                   key='synthetic',
                                   columns=columns,
                                   where='index in index')
            else:
                return pd.read_hdf(self.synthetic_population.pop_file,
                                   key='synthetic',
                                   columns=columns)
            
        if type(index) is np.ndarray:
            index = index.tolist()
            
        columns = self.synthetic_population.head(0).columns.to_numpy()
        if var not in columns:
            # compute some useful quantities that are not defined in df_synthetic
            if var == 'q':
                with pd.HDFStore(self.synthetic_population.pop_file, 'r') as store:
                    df = read_data(['S1_mass', 'S2_mass'], index)
                    values = self.mass_ratio(df["S1_mass"], df["S2_mass"])
            elif var == 'm_tot':
                df = read_data(['S1_mass', 'S2_mass'], index)
                values = df["S1_mass"] + df["S2_mass"]
            elif var == 'm_chirp':
                df = read_data(['S1_mass', 'S2_mass'], index)
                values = self.m_chirp(df["S1_mass"], df["S2_mass"])
            elif var == 'chi_eff':
                # check if tilts are provided
                if ("S1_spin_orbit_tilt" not in columns or 
                    "S2_spin_orbit_tilt" not in columns):
                    warnings.warn('Spin tilts not in dataset! Assuming spins '
                                   'are aligned to orbital angular momentum '
                                   'to compute chi_eff!')
                    tilt_BH1 = np.zeros(len(self.synthetic_population.indices))
                    tilt_BH2 = tilt_BH1
                else:
                    df = read_data(['S1_mass', 'S2_mass', 'S1_spin', 'S2_spin',
                                    'S1_spin_orbit_tilt', 'S2_spin_orbit_tilt'],
                                   index)
                    
                    tilt_BH1 = df["S1_spin_orbit_tilt"]
                    tilt_BH2 = df["S2_spin_orbit_tilt"]
                values = self.chi_eff(df["S1_mass"],
                                      df["S2_mass"],
                                      df["S1_spin"],
                                      df["S2_spin"],
                                      tilt_BH1,
                                      tilt_BH2)
            elif var == 'underlying_mass_for_met':
                df = read_data(['metallicity'], index)
                for met, underlying_mass in zip(self.synthetic_population.metallicities, self.synthetic_population.underlying_mass_per_met):
                    df[df['metallicity'] == met] = underlying_mass
                values = df['metallicity'].astype(np.float64).to_numpy()
            else:
                raise ValueError(f'Definition for {var} not available!')
        else:
            df = read_data([var], index)
            values =  df[var].to_numpy()
        return values

    def merger_rate_weight(self, z_birth, z_merger, p_det, i):
        """Compute the merger rate weight (w_ijk) in yr^-1 units.

        Parameters
        ----------
        z_birth : array doubles
            Cosmological redshift of formation of the binary systems. This MUST
            be the same for all binaries.
        z_merger : array doubles
            Cosmological redshift of merger of the binary systems.
        p_det : array doubles
            Detection probability of the binary.
        i : array integers
            Indicies correponding to the binaries you want to select out of the
            population.

        Returns
        -------
        array doubles
            Return the cosmological weights w_ijk (detection rate contibution)
            of the binary k in the metallicity bin j born at redshift birth i
            as in Eq. (B.8) of Bavera et al. (2020).

        """
        # get SFR at a given population birth redshift
        SFR_at_z_birth = star_formation_rate(self.MODEL['SFR'], z_birth)

        # get metallicity bin edges
        met_bins = self.edges_metallicity_bins

        # distribute the DCO to the corresponding bin
        binplace = np.digitize(self.get_data("metallicity", i), met_bins)

        # compute fSFR assuming log-normal distributed metallicities
        fSFR = fractional_SFR_at_given_redshift(z_birth, SFR_at_z_birth,
                                                self.MODEL['SFR'], self.MODEL['sigma_SFR'],
                                                met_bins, binplace,
                                                self.MODEL['Z_max'], self.MODEL['select_one_met'])

        # simulated mass per given metallicity corrected for the unmodeled
        # single and binary stellar mass
        M_model = self.get_data("underlying_mass_for_met",i)

        # speed of light
        c = const.c.to('Mpc/yr').value  # Mpc/yr

        # delta cosmic time bin
        deltaT = self.MODEL['delta_t'] * 10 ** 6  # yr

        # comoving distance corresponding to the DCO merger redshift
        D_c = self.get_comoving_distance_from_redshift(z_merger)

        # DCO cosmological weights B.8 in Bavera et al. (2020)
        return 4.*np.pi * c * D_c**2 * p_det * deltaT * fSFR / M_model # yr^-1
    
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
            if isinstance(self.dlogZ, float):
                bin_met[0] = met_val[0] - self.dlogZ / 2.
                bin_met[-1] = met_val[0] + self.dlogZ / 2.
            elif isinstance(self.dlogZ, list) or isinstance(self.dlogZ, np.array):
                bin_met[0] = self.dlogZ[0]
                bin_met[-1] = self.dlogZ[1]
        
        return 10**bin_met
 
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
    
    def compute_rate_weights(self, sensitivity='infinite', overwrite=False):
        
        if sensitivity != 'infinite' and self.synthetic_population.population_selection != 'DCO':
            raise ValueError("You can only compute non-infinite sensitivity for a DCO population.")
        elif sensitivity != 'infinite' and self.synthetic_population.population_selection == 'DCO' :
            self.compute_sensitivity_rate_weights(sensitivity, overwrite=overwerite)
        elif sensitivity == 'infinite':
            self.compute_sensitivity_rate_weights(sensitivity='infinite', overwrite=overwrite)        
        else:
            raise ValueError("Sensitivity not supported!")
        
    def compute_sensitivity_rate_weights(self, sensitivity='infinite', overwrite=False):
        """Compute the cosmological weights of the DCO population.

        This function will create a directory path_to_dir/DCOs/sensitivity where
        it will save the weigths, binary indicies k, z_formation, z_merger.
        This is needed for scalability as a ~100k DCO population generates ~10M
        non-zero weights assuming a delta_t=100Myr.

        Parameters
        ----------
        sensitivity : string
            GW detector sensitivity and network configuration you want to use, see arXiv:1304.0670v3
            detector sensitivities are taken from: https://dcc.ligo.org/LIGO-T2000012-v2/public
            available sensitivity keys (for Hanford, Livingston, Virgo network):
                'O3actual_H1L1V1' : aligo_O3actual_H1.txt, aligo_O3actual_L1.txt, avirgo_O3actual.txt
                'O4low_H1L1V1' : aligo_O4low.txt, aligo_O4low.txt, avirgo_O4high_NEW.txt
                'O4high_H1L1V1' : aligo_O4high.txt, aligo_O4high.txt, avirgo_O4high_NEW.txt
                'design_H1L1V1' : AplusDesign.txt, AplusDesign.txt, avirgo_O5high_NEW.txt
                'infinite': intrinsic merging DCO population, i.e. p_det = 1
        flag_pdet : bool
            This is a control variable. In order to be sure you want to run
            infinite sensitivity set `p_det=False`.
        path_to_dir : string
            Path to the workingn directory where you want to store the
            cosmological weights.

        """
        
        PATH_TO_PDET_GRID = os.path.join(PATH_TO_POSYDON_DATA, 'POSYDON_data/selection_effects/pdet_grid.hdf5')
        pop = self.synthetic_population.population_selection
        # check if the rates have been computed
        with pd.HDFStore(self.rates_file, 'r') as store:
            if f'/rates/{self.synthetic_population.population_selection}/{sensitivity}/weights' in store.keys():
                print('The rates for this sensitivity have already been computed!')
                if overwrite:
                    print('Overwriting the rates...')
                    store.remove(f'/rates/{pop}/{sensitivity}/weights')
                    store.remove(f'/rates/{pop}/{sensitivity}/z_events')
                return

        if sensitivity != 'infinite':
            self.sel_eff = selection_effects.KNNmodel(grid_path=PATH_TO_PDET_GRID,
                                                      sensitivity_key=sensitivity)
        
        # index of all events
        indices = self.synthetic_population.indices    
        
                # redshift of each time bin
        z_birth = self.centers_redshift_bins
        
        # define interpolator
        get_redshift_from_time = self.redshift_from_cosmic_time_interpolator
        
        for i in tqdm(range(0, len(indices), self.chunksize), desc='chunk loop'):
    
            selected_indices = indices[i:i+self.chunksize]
            i_selected_indices = np.arange(0, (len(selected_indices)))
            
            # create the dataframes to store the merger redshift and weights
            # -----> z_birth
            # |
            # v
            # indices
            # Creating for each chunk allows us to store the data in an HDF5 file
            # We cannot add columns to an existing HDF5 file
            
            df_z_events = pd.DataFrame(index=selected_indices, dtype=np.float64,
                                        columns=[i for i in range(len(z_birth))])
            df_weights = pd.DataFrame(index=selected_indices, dtype=np.float64,
                                        columns=[i for i in range(len(z_birth))])
            
            with pd.HDFStore(self.synthetic_population.pop_file, 'r') as store:
                event_time = store.select_column('synthetic', 'time', start=i, stop=i+self.chunksize)
                if  't_delay' in store.select('synthetic', start=0, stop=0).columns:
                    t_delay = store.select_column('synthetic', 't_delay', start=i, stop=i+self.chunksize)
                    event_time += t_delay
                    del t_delay
                    
                
            print('looping now')
            for z in tqdm(range(len(z_birth)), desc='redshift loop'):
                
                t_event = self.get_cosmic_time_at_event(z_birth[z], event_time)
                bool_event = t_event <= cosmology.age(1e-08).value*0.9999999
                
                if np.sum(bool_event) == 0:
                    continue
            
                # get the redshift of the event        
                t_event = t_event[bool_event]
                
                z_events = get_redshift_from_time(t_event)
                
                hubble_time_indices = selected_indices[bool_event]
                i_hubble_time_indices = i_selected_indices[bool_event]
                
                p_det = np.ones(len(z_events))
                
                if sensitivity != 'infinite':
                    # get detection probabilities
                    data_slice = pd.DataFrame()
                    with pd.HDFStore(self.synthetic_population.pop_file, 'r') as store:
                        
                        data_slice['m1'] = store.select_column('synthetic', 'S1_mass', start=i, stop=i+self.chunksize)[bool_event]
                        
                        data_slice['m2'] = store.select_column('synthetic', 'S2_mass', start=i, stop=i+self.chunksize)[bool_event]
                        data_slice['q'] = self.mass_ratio(data_slice['m1'], data_slice['m2'])
                        data_slice['z'] = z_events
                        data_slice['chieff'] = self.chi_eff(data_slice['m1'],
                                                            data_slice['m2'],
                                                            store.select_column('synthetic', 'S1_spin', start=i, stop=i+self.chunksize)[bool_event],
                                                            store.select_column('synthetic', 'S2_spin', start=i, stop=i+self.chunksize)[bool_event],
                                                            store.select_column('synthetic', 'S1_spin_orbit_tilt', start=i, stop=i+self.chunksize)[bool_event],
                                                            store.select_column('synthetic', 'S2_spin_orbit_tilt', start=i, stop=i+self.chunksize)[bool_event])
                    
                    data_slice.fillna(0, inplace=True)
                    # get the detection probability
                    p_det = self.sel_eff.predict_pdet(data_slice)
                    
                bool_det = p_det > 0.
                if np.sum(bool_det) == 0:
                    continue
                    # also slice on the detection probability now!
                z_events = z_events[bool_det]
                
                detectable_indices = hubble_time_indices[bool_det]
                i_detectable_indices = i_hubble_time_indices[bool_det]
                
                p_det = p_det[bool_det]
                        
                w_ijk = self.merger_rate_weight(z_birth[z] * np.ones(len(z_events)),
                                                z_events,
                                                p_det,
                                                detectable_indices)
                
                df_weights.iloc[i_detectable_indices, z] = w_ijk
                df_z_events.iloc[i_detectable_indices, z] = z_events
            
            pop = self.synthetic_population.population_selection
            # store the dataframes in the rates file
            with pd.HDFStore(self.rates_file, 'a') as store:
                store.append(f'rates/{pop}/{sensitivity}/z_events',
                                df_z_events,
                                format='table')
                store.append(f'rates/{pop}/{sensitivity}/weights',
                                df_weights,
                                format='table')

        with pd.HDFStore(self.rates_file, 'a') as store:        
            store.put(f'rates/{pop}/{sensitivity}/z_birth',pd.Series(z_birth))
            store.put(f'rates/{pop}/{sensitivity}/MODEL',  pd.Series(self.MODEL))

    def get_shell_comovig_volume(self, z_hor_i, z_hor_f, sensitivity='infinite'):
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
        
    def compute_rate_density(self, mt_channels=False, sensitivity='infinite'):
        """Compute the GRB/DCO rate density.

        Parameters
        ----------
        w_ijk : array doubles
            Cosmological weights computed with Eq. B.8 of Bavera et at. (2020).
        z_event : array doubles
            Cosmolgocial redshift of the event you are tracking.
        observable : string
            Event you are tracking, available:
            'DCOs': merger event of a DCO system
            'GRB1': gamma ray bursts of star 1
            'GRB2': gamma ray bursts of star 2
        sensitivity : string
            This takes into account the detector sensitivity, available:
            'infinite': p_det = 1
            'beamed': TODO for GRBs

        Returns
        -------
        array doubles
            Return the DCOs merger rate density (Gpc^-3 yr^-1) as a function of
            cosmolgocial redshift as in Eq. (D.1) in Bavera et al. (2022)
            arXiv:2106.15841

        """
        
            
            
        pop = self.synthetic_population.population_selection
        
        if sensitivity != 'beaming':
                    # check if the rates have been computed
            with pd.HDFStore(self.rates_file, 'r') as store:
                if f'/rates/{pop}/{sensitivity}/weights' not in store.keys():
                    print('The rates for this sensitivity have not been computed yet!')
        
            w_ijk = pd.read_hdf(self.rates_file,
                                key=f'rates/{pop}/{sensitivity}/weights')
            z_events = pd.read_hdf(self.rates_file,
                                   key=f'rates/{pop}/{sensitivity}/z_events')
        else: 
                    # check if the rates have been computed
            with pd.HDFStore(self.rates_file, 'r') as store:
                if f'/rates/{pop}/infinite/weights' not in store.keys():
                    print('The rates for this sensitivity have not been computed yet!')
        
            w_ijk = pd.read_hdf(self.rates_file,
                                key=f'rates/{pop}/infinite/weights')
            z_events = pd.read_hdf(self.rates_file,
                                   key=f'rates/{pop}/infinite/z_events')
            
            
        z_horizon = self.edges_redshift_bins
        n = len(z_horizon)
                
        # compute the rate density
        self.R = {'total' : np.zeros(n-1)}
        
        if mt_channels:
            channels = pd.read_hdf(self.synthetic_population.pop_file,
                                  key='synthetic',
                                  columns=['channel']).to_numpy()
            for ch in np.unique(channels):
                self.R[ch] = np.zeros(n-1)
        if self.verbose:
            print('Computing the rate density...')

        for i in range(1, n):
            # select events in the redshift bin
             # compute Eq. (D.1) in Bavera et al. (2022) arXiv:2106.15841
            bool_events = (z_events > z_horizon[i-1]) & (z_events <= z_horizon[i])
            if pop == 'DCO':
                # sum over all events in the redshift bin
                self.R['total'][i-1] = np.nansum(w_ijk[bool_events].to_numpy())\
                            / self.get_shell_comovig_volume(z_horizon[i-1],
                                                            z_horizon[i],
                                                            sensitivity)
            if pop == 'GRB':

                if sensitivity == 'beaming':

                    with pd.HDFStore(self.synthetic_population.pop_file, 'r') as store:
                        S1_beaming = np.nan_to_num(store.select_column('synthetic', 'S1_f_beaming'))
                        S2_beaming = np.nan_to_num(store.select_column('synthetic', 'S2_f_beaming'))
                    
                    sensitivity='infinite'

                    w_ijk = (w_ijk * S1_beaming[:,None]) + (w_ijk * S2_beaming[:,None])

                    
                self.R['total'][i-1] = np.nansum(w_ijk[bool_events].to_numpy())\
                            / self.get_shell_comovig_volume(z_horizon[i-1],
                                                            z_horizon[i],
                                                            'infinite')
                        
        if mt_channels:
            if self.verbose:
                print('Computing the rate density per formation channel...')
            # compute the rate density per formation channel
            for ch in np.unique(channels):
                if self.verbose:
                    print(f'Channel {ch}...')
                bool_ch = (channels == ch)
                channel_events = z_events[bool_ch]
                channel_weights = w_ijk[bool_ch]
                for i in range(1, n):
                    bool_events = (channel_events > z_horizon[i-1]) & (channel_events <= z_horizon[i])
                    self.R[ch][i-1] = np.nansum(channel_weights[bool_events].to_numpy())\
                                    / self.get_shell_comovig_volume(z_horizon[i-1],
                                                                    z_horizon[i],
                                                                    sensitivity)
        return self.R

    
    def resample_synthetic_population(self, output_file=None, sensitivity='infinite', export_cols=None):
        '''Function to resample the synthetic population based on their weights'''
            # index, z_formation, z_event, w_ijk, export_col, pop, reset_grb_properties
        
        if output_file is None:
            output_file = self.rates_file
            
        pop = self.synthetic_population.population_selection
            
        # check if sensitivity has been computed
        with pd.HDFStore(self.rates_file, 'r') as store:
            if f'/rates/{pop}/{sensitivity}/weights' not in store.keys():
                print('The rates for this sensitivity and population have not been computed yet!')
            
        # only select non-zero weights
        w_ijk = pd.read_hdf(self.rates_file,
                              key=f'rates/{pop}/{sensitivity}/weights')
        z_events = pd.read_hdf(self.rates_file,
                              key=f'rates/{pop}/{sensitivity}/z_events')
        
        # flatten all events and weights into a single dataframe
        w_ijk_flat = w_ijk.to_numpy(na_value=0).flatten()
        indices = w_ijk.index.to_numpy()
        indices_flat = np.repeat(indices, len(w_ijk.columns))
        z_formation_flat = np.repeat(r.centers_redshift_bins, len(indices))
        z_event_flat = z_events.to_numpy(na_value=0).flatten()
        
        # select events with non-zero weights
        bool_weights = w_ijk_flat > 0.
        df = pd.DataFrame(index=indices_flat[bool_weights],
                  data={'z_formation': z_formation_flat[bool_weights],
                        'weight': w_ijk_flat[bool_weights],
                        'z_event': z_event_flat[bool_weights]})
        
        # adding additional data increases the memory usage significantly
        save_cols = ['metallicity','time','t_delay','S1_state','S2_state',
                     'S1_mass','S2_mass','S1_spin','S2_spin',
                     'orbital_period','eccentricity', 'q', 'm_tot',
                     'm_chirp', 'chi_eff']
        
        columns = self.synthetic_population.head(0).columns
        if "channel" in columns:
            save_cols.append('channel')
        if export_cols is not None:
            for c in export_cols:
                if c not in save_cols:
                    save_cols.append(c)
        
        for c in save_cols:            
        # initialise the columns to allow pd.update to work
            if c not in df.columns:
                df[c] = pd.NA
                    
        # should be the same as the original indices
        indices = np.unique(df.index)
        if self.verbose: print("adding the following columns:")
        
        for c in save_cols:
            if self.verbose: print(c)
                
            values =  self.get_data(c, indices)
            # update repeating indices with the same value
            tmp_df = pd.DataFrame(index=indices, data={c: values})
            df.update(tmp_df[c])

        if self.verbose: print("exporting the population...")
        df = df.fillna(0.)
        
        with pd.HDFStore(output_file, 'a') as store:
            store.append(f'rates/{pop}/{sensitivity}/sampled_population',
                         df,
                         format='table',
                         data_columns=True)

        return Population(output_file)

    def plot_rate_density(self, DCO=False, GRB=False, **kwargs):
        """Plot DCO and GRB rate densities."""
        if not DCO and not GRB:
            raise ValueError('You need to choose at least one population to plot!')
        if DCO:
            #if self.dco_z_rate_density is None or self.dco_rate_density is None:
            #    raise ValueError('First you need to compute the merger rate density!')
            #else:
            z_dco = self.centers_redshift_bins
            rate_dco =  self.R
        else:
            z_dco = None
            rate_dco = None
        if GRB:
            #if self.grb_z_rate_density is None or self.grb_rate_density is None:
            #    raise ValueError('First you need to compute the GRB rate density!')
            #else:
            z_grb = self.centers_redshift_bins
            rate_grb =  self.R
        else:
            z_grb = None
            rate_grb = None
        plot_pop.plot_rate_density(z_dco, rate_dco, z_grb, rate_grb, **kwargs)
        




class Rates2(object):
    """Compute DCO rates."""
    def __init__(self, df, verbose=False, **kwargs):
        """Compute DCO rates.

        Parameters
        ----------
        df : object
            Pandas dataframe containing the synthetic binary population.
        delta_t : float
            Cosmic time bin size used to bin distribute the synthetc binary
            population.
        SFR : string
            Star-formation-rate history you want to use:
            'Madau+Dickinson14', 'Madau+Fragos17', 'Neijssel+19', 'IllustrisTNG'
        sigma_SFR : string
            Standard deviation of the trucated log-normal metallicity distribution
            assumed for the star-formation rate history in the case you select
            'Madau+Dickinson14', 'Madau+Fragos17', 'Neijssel+19',
        Z_max : float
            The absolute maximum metallicity of the star formation rate history.
            NOTE: This should be used when assuming a truncated log-normal
            distribution of metallcities to not overpredict 1, also it can be
            used with the IllustrisTNG SFR to prevent unphysically high
            metallicities.
        select_one_met : bool
            If `True` your synthetic binary population should contain only one
            descrete metallicity.
        dlogZ : float or range
            Used when `select_one_met=True` to select the metallicity range
            around which you want to integrate the SFR history. If float value
            value is provided we assume a symmetric range dlogZ/2 around the
            provided metallcity else we consider the range provided
            [Z_min, Z_max] which could also be asymmetric.
        verbose : bool
            `True` if you want the print statements.


        TODO: add the definitios of the variables used by this class.

        """

        self.df = df
        self.verbose = verbose

        if kwargs:
            for key in kwargs:
                if key not in DEFAULT_MODEL:
                    raise ValueError(key + " is not a valid parameter name!")
            for varname in DEFAULT_MODEL:
                default_value = DEFAULT_MODEL[varname]
                setattr(self, varname, kwargs.get(varname, default_value))
        else:
            for varname in DEFAULT_MODEL:
                default_value = DEFAULT_MODEL[varname]
                setattr(self, varname, default_value)
                
        if self.compute_GRB_properties:
            for key in GRB_PROPERTIES:
                if key not in self.df:
                    warnings.warn("GRB properties not in the dataset, "
                                  "computing them right now!")
                    self.df = get_GRB_properties(self.df, 
                                                 self.GRB_efficiency, 
                                                 self.GRB_beaming,
                                                 self.E_GRB_iso_min)

    ######################################################
    ###   DCO detection rate and merger rate density   ###
    ###    see arXiv:1906.12257, arXiv:2010.16333,     ###
    ###        arXiv:2106.15841, arXiv:221210924,      ###
    ###        arXiv:2204.02619, arXiv:2011.10057      ###
    ###              arXiv:2012.02274                  ###
    ######################################################

    def chi_eff(self, m_1, m_2, a_1, a_2, tilt_1, tilt_2):
        return (m_1*a_1*np.cos(tilt_1)+m_2*a_2*np.cos(tilt_2))/(m_1+m_2)

    def m_chirp(self, m_1, m_2):
        return (m_1*m_2)**(3./5)/(m_1+m_2)**(1./5)

    def mass_ratio(self, m_1, m_2):
        q = m_2/m_1
        q[q>1.] = 1./q[q>1.]
        return q

    def get_data(self, var, index=None):
        """Resample class elements.

        Parameters
        ----------
        var : string
            Key of self.df. This method support a list of custom definitions for
            DCO studies like q, m_tot, m_chirp, chi_eff.
        index : list of int
            List of indicies of the binaries of interest.

        Returns
        -------
        array
            Array containing the self.df[var][index].

        """
        if var not in self.df:
            # compute some useful quantities that are not defined in df_synthetic
            if var == 'q':
                values = self.mass_ratio(self.df["S1_mass"], self.df["S2_mass"])
            elif var == 'm_tot':
                values = self.df["S1_mass"] + self.df["S2_mass"]
            elif var == 'm_chirp':
                values = self.m_chirp(self.df["S1_mass"], self.df["S2_mass"])
            elif var == 'chi_eff':
                # check if tilts are provided
                if ("S1_spin_orbit_tilt" not in self.df or 
                    "S2_spin_orbit_tilt" not in self.df):
                    warnings.warn('Spin tilts not in dataset! Assuming spins '
                                   'are aligned to orbital angular momentum '
                                   'to compute chi_eff!')
                    tilt_BH1 = np.zeros(len(self.df["S1_spin"]))
                    tilt_BH2 = tilt_BH1
                else:
                    tilt_BH1 = self.df["S1_spin_orbit_tilt"]
                    tilt_BH2 = self.df["S2_spin_orbit_tilt"]
                values = self.chi_eff(self.df["S1_mass"],
                                      self.df["S2_mass"],
                                      self.df["S1_spin"],
                                      self.df["S2_spin"],
                                      tilt_BH1,
                                      tilt_BH2)
            else:
                raise ValueError(f'Definition for {var} not available!')
        else:
            values =  self.df[var]

        if index is not None:
            return values[index].values
        else:
            return values.values

    def get_redshift_from_cosmic_time_interpolator(self):
        """Interpolator to compute the cosmological redshift given the cosmic time.

        Returns
        -------
        object
            Returns the trained interpolator object.

        """
        # astropy z_at_value method is too slow to compute z_mergers efficinty
        # we must implment interpolation
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
        interpolator = self.get_redshift_from_cosmic_time_interpolator()
        return interpolator(t_cosm)

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

    def get_comoving_disntance_from_redshift(self, z):
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

    def get_cosmic_time_at_dco_merger(self, z_birth):
        """Get cosmic time at DCO merger.

        Parameters
        ----------
        z_birth : double
            Cosmological redshift of formation of the DCO system
            (must be the same for every binary).

        Returns
        -------
        double
            Cosmic time in Gyr at DCO merger for all binaries born at z_birth.

        """
        n = self.df.shape[0]
        t_birth = self.get_cosmic_time_from_redshift(z_birth) * np.ones(n) # Gyr
        return t_birth + (self.df["time"] + self.df["t_delay"]) * 10 ** (-3)  #Gyr

    def get_redshift_at_dco_merger(self, z_birth):
        """Get redshift of merger of DCOs.

        Parameters
        ----------
        z_birth : double
            Redshift of formation of the DCO system (must be the same for every
            binary).

        Returns
        -------
        double
            Redshift of merger of DCOs born at z_birth.

        """
        n = self.df.shape[0]
        t_merger = self.get_cosmic_time_at_dco_merger(z_birth)
        z_merger = np.ones(n) * np.nan
        bool_merger = t_merger < self.get_cosmic_time_from_redshift(0.) * np.ones(n)  # check if the binary merges
        z_merger[bool_merger] = z_at_value(cosmology.age, t_merger[bool_merger] * u.Gyr)
        return z_merger

    def get_centers_metallicity_bins(self):
        """Return the centers of the metallicity bins.

        Returns
        -------
        array double
            Returns sampled metallicities of the populattion. This correponds
            to the center of each metallicity bin.

        """
        return np.unique(self.df['metallicity'].values)

    def get_edges_metallicity_bins(self):
        """Return the edges of the metallicity bins.

        Returns
        -------
        array double
            Returns the edges of all metallicity bins. We assume metallicities
            were binned in log-space.

        """
        met_val = np.log10(self.get_centers_metallicity_bins())
        bin_met = np.zeros(len(met_val)+1)
        # if more than one metallicty bin
        if len(met_val) > 1 :
            bin_met[0] = met_val[0] - (met_val[1] - met_val[0]) / 2.
            bin_met[-1] = met_val[-1] + (met_val[-1] - met_val[-2]) / 2.
            bin_met[1:-1] = met_val[:-1] + (met_val[1:] - met_val[:-1]) / 2.
        # one metallicty bin
        elif len(met_val) == 1 :
            if isinstance(self.dlogZ, float):
                bin_met[0] = met_val[0] - self.dlogZ / 2.
                bin_met[-1] = met_val[0] + self.dlogZ / 2.
            elif isinstance(self.dlogZ, list) or isinstance(self.dlogZ, np.array):
                bin_met[0] = self.dlogZ[0]
                bin_met[-1] = self.dlogZ[1]

        return 10**bin_met

    def get_centers_redshift_bins(self):
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
        self.n_redshift_bin_centers = int(cosmology.age(0).to('Myr').value/self.delta_t)
        for i in range(self.n_redshift_bin_centers+1):
            t_birth.append(t_birth_bin[i] - self.delta_t*1e-3/2.) # Gyr
            t_birth_bin.append(t_birth_bin[i] - self.delta_t*1e-3) # Gyr
        t_birth = np.array(t_birth)
        # compute the redshift
        z_birth = []
        for i in range(self.n_redshift_bin_centers+1):
            z_birth.append(z_at_value(cosmology.age, t_birth[i] * u.Gyr))
        z_birth = np.array(z_birth)

        return z_birth

    def get_edges_redshift_bins(self):
        """Compute redshift bin edges.

        Returns
        -------
        array doubles
            We devide the cosmic time history of the Universe in equally spaced
            bins of cosmic time of self.delta_t (100 Myr default) an compute the
            redshift corresponding to edges of these bins.

        """
        # generate t_birth at the middle of each self.delta_t bin
        t_birth_bin = [cosmology.age(0.).value]
        for i in range(self.n_redshift_bin_centers+1):
            t_birth_bin.append(t_birth_bin[i] - self.delta_t*1e-3) # Gyr
        # compute the redshift
        z_birth_bin = []
        for i in range(self.n_redshift_bin_centers):
            # do not count first edge, we add z=0. later
            z_birth_bin.append(z_at_value(cosmology.age, t_birth_bin[i+1] * u.Gyr))
        # add the first and last bin edge at z=0. and z=inf.=100
        z_birth_bin = np.array([0.]+z_birth_bin+[100.])

        return z_birth_bin

    def merger_rate_weight(self, z_birth, z_merger, p_det, i):
        """Compute the merger rate weight (w_ijk) in yr^-1 units.

        Parameters
        ----------
        z_birth : array doubles
            Cosmological redshift of formation of the binary systems. This MUST
            be the same for all binaries.
        z_merger : array doubles
            Cosmological redshift of merger of the binary systems.
        p_det : array doubles
            Detection probability of the binary.
        i : array integers
            Indicies correponding to the binaries you want to select out of the
            population.

        Returns
        -------
        array doubles
            Return the cosmological weights w_ijk (detection rate contibution)
            of the binary k in the metallicity bin j born at redshift birth i
            as in Eq. (B.8) of Bavera et al. (2020).

        """
        # get SFR at a given population birth redshift
        SFR_at_z_birth = star_formation_rate(self.SFR, z_birth)

        # get metallicity bin edges
        met_bins = self.get_edges_metallicity_bins()

        # distribute the DCO to the corresponding bin
        binplace = np.digitize(self.get_data("metallicity", i), met_bins)

        # compute fSFR assuming log-normal distributed metallicities
        fSFR = fractional_SFR_at_given_redshift(z_birth, SFR_at_z_birth,
                                                self.SFR, self.sigma_SFR,
                                                met_bins, binplace,
                                                self.Z_max, self.select_one_met)

        # simulated mass per given metallicity corrected for the unmodeled
        # single and binary stellar mass
        M_model = self.get_data("underlying_mass_for_met",i)

        # speed of light
        c = const.c.to('Mpc/yr').value  # Mpc/yr

        # delta cosmic time bin
        deltaT = self.delta_t * 10 ** 6  # yr

        # comoving distance corresponding to the DCO merger redshift
        D_c = self.get_comoving_disntance_from_redshift(z_merger)

        # DCO cosmological weights B.8 in Bavera et al. (2020)
        return 4.*np.pi * c * D_c**2 * p_det * deltaT * fSFR / M_model # yr^-1


    def compute_merger_rate_weights(self, sensitivity, flag_pdet=True, path_to_dir='./', extention='npz'):
        """Compute the cosmological weights of the DCO population.

        This function will create a directory path_to_dir/DCOs/sensitivity where
        it will save the weigths, binary indicies k, z_formation, z_merger.
        This is needed for scalability as a ~100k DCO population generates ~10M
        non-zero weights assuming a delta_t=100Myr.

        Parameters
        ----------
        sensitivity : string
            GW detector sensitivity and network configuration you want to use, see arXiv:1304.0670v3
            detector sensitivities are taken from: https://dcc.ligo.org/LIGO-T2000012-v2/public
            available sensitivity keys (for Hanford, Livingston, Virgo network):
                'O3actual_H1L1V1' : aligo_O3actual_H1.txt, aligo_O3actual_L1.txt, avirgo_O3actual.txt
                'O4low_H1L1V1' : aligo_O4low.txt, aligo_O4low.txt, avirgo_O4high_NEW.txt
                'O4high_H1L1V1' : aligo_O4high.txt, aligo_O4high.txt, avirgo_O4high_NEW.txt
                'design_H1L1V1' : AplusDesign.txt, AplusDesign.txt, avirgo_O5high_NEW.txt
                'infinite': intrinsic merging DCO population, i.e. p_det = 1
        flag_pdet : bool
            This is a control variable. In order to be sure you want to run
            infinite sensitivity set `p_det=False`.
        path_to_dir : string
            Path to the workingn directory where you want to store the
            cosmological weights.

        """
        # check if the folder three exists, otherwise create it
        DCO_dir = os.path.join(path_to_dir,'DCOs')
        if 'DCOs' not in os.listdir(path_to_dir):
            os.makedirs(DCO_dir)

        sensitivity_dir = os.path.join(DCO_dir, f'{sensitivity}_sensitivity')
        if f'{sensitivity}_sensitivity' not in os.listdir(DCO_dir):
            os.makedirs(sensitivity_dir)

        # index of all DCOs
        n = self.df.shape[0]
        index = np.arange(0, n)

        # redshif of each time bin
        z_birth = self.get_centers_redshift_bins()

        # define interpolator
        get_redshift_from_time = self.get_redshift_from_cosmic_time_interpolator()

        # hashmap to store everythig
        data = {'index': {}, 'z_merger': {}, 'weights': {}}
        
        # load and store detector selection effects interpolator
        if sensitivity != 'infinite':
            self.sel_eff = selection_effects.KNNmodel(grid_path=PATH_TO_PDET_GRID,
                                                      sensitivity_key=sensitivity)

        # loop over all redshift bins
        for i in tqdm(range(len(z_birth))):

            # compute the merger time of each DCOs
            t_merger = self.get_cosmic_time_at_dco_merger(z_birth[i])

            # detectable population
            if flag_pdet == True and sensitivity != 'infinite':

                # sort out systems not merging within the Hubble time
                bool_merger = t_merger < cosmology.age(1e-08).value*0.9999999 * np.ones(n)

                # if there are no merging DCO, continue
                if len(index[bool_merger]) == 0:
                    data['index'][str(z_birth[i])] = np.array([])
                    data['z_merger'][str(z_birth[i])] = np.array([])
                    data['weights'][str(z_birth[i])] = np.array([])
                    continue

                # get some quantities to compute detection probabilities
                data_slice = pd.DataFrame()
                data_slice['m1'] = self.get_data('S1_mass', index[bool_merger])
                data_slice['q'] = self.get_data('q', index[bool_merger])
                z_m = get_redshift_from_time(t_merger[bool_merger])
                data_slice['z'] = z_m
                data_slice['chieff'] = self.get_data('chi_eff', index[bool_merger])

                # compute detection probabilities
                p_det = self.sel_eff.predict_pdet(data_slice)

                # sort out undetectable sources
                bool_detect = p_det > 0.

                # if there are no detectable DCO, continue
                if len(index[bool_merger][bool_detect]) == 0:
                    data['index'][str(z_birth[i])] = np.array([])
                    data['z_merger'][str(z_birth[i])] = np.array([])
                    data['weights'][str(z_birth[i])] = np.array([])
                    continue

                # store index and redshift of merger
                data['index'][str(z_birth[i])] = index[bool_merger][bool_detect]
                data['z_merger'][str(z_birth[i])] = z_m[bool_detect]

                # compute and store marger rate weights as in eq. B.8 in Bavera et al. (2020)
                z_b = np.ones(len(index[bool_merger][bool_detect]))*z_birth[i]
                w_ijk = self.merger_rate_weight(z_b, z_m[bool_detect],
                                                p_det[bool_detect],
                                                index[bool_merger][bool_detect])
                data['weights'][str(z_birth[i])] = w_ijk

            # intrinsic population
            elif flag_pdet == False and sensitivity == 'infinite':

                # sort out systems not merging within the Hubble time
                bool_merger = t_merger < cosmology.age(1e-08).value*0.9999999 * np.ones(n)

                # if there are no merging DCO, continue
                if len(index[bool_merger]) == 0:
                    data['index'][str(z_birth[i])] = np.array([])
                    data['z_merger'][str(z_birth[i])] = np.array([])
                    data['weights'][str(z_birth[i])] = np.array([])
                    continue

                # get merger redshifts
                z_m = get_redshift_from_time(t_merger[bool_merger])

                # set detection probabilities to 1
                p_det = np.ones(len(index[bool_merger]))

                # store index and redshift of merger
                data['index'][str(z_birth[i])] = index[bool_merger]
                data['z_merger'][str(z_birth[i])] = z_m

                # compute and store marger rate weights as in eq. B.8 in Bavera et al. (2020)
                z_b = np.ones(len(index[bool_merger]))*z_birth[i]
                w_ijk = self.merger_rate_weight(z_b, z_m, p_det, index[bool_merger])
                data['weights'][str(z_birth[i])] = w_ijk

            else:
                raise ValueError('Missmatch between sensitivity and flag_pdet!')

        # TODO: implemet the saving to h5 files
        if extention == 'npz':
            if self.verbose:
                print('Formatting the data ....')
            data_to_save = [[],[],[],[]]
            for i, dict in enumerate([data[key] for key in data.keys()]):
                for key in dict.keys():
                    data_to_save[i].extend(dict[key].tolist())
                    if i == 0:
                        data_to_save[3].extend(np.ones(len(dict[key]))*float(key))
            if self.verbose:
                print('Saving the data ....')
            for i, key in enumerate(data.keys()):
                if key == 'index':
                    fmt_str = '%i'
                else:
                    fmt_str = '%.8E'
                np.savez(os.path.join(sensitivity_dir, f"{key}.npz"),
                         key=data_to_save[i], fmt=fmt_str)

            np.savez(os.path.join(sensitivity_dir,"z_formation.npz"),
                     key=data_to_save[3], fmt='%.8E')
        else:
            raise ValueError('Extension not supported!')

    def load_merger_rate_weights(self, sensitivity, path_to_dir='./', extention='npz'):
        """Load the cosmological weights of the DCO populatio.

        Parameters
        ----------
        sensitivity : string
            GW detector sensitivity and network configuration you want to use, see arXiv:1304.0670v3
            detector sensitivities are taken from: https://dcc.ligo.org/LIGO-T2000012-v2/public
            available sensitivity keys (for Hanford, Livingston, Virgo network):
                'O3actual_H1L1V1' : aligo_O3actual_H1.txt, aligo_O3actual_L1.txt, avirgo_O3actual.txt
                'O4low_H1L1V1' : aligo_O4low.txt, aligo_O4low.txt, avirgo_O4high_NEW.txt
                'O4high_H1L1V1' : aligo_O4high.txt, aligo_O4high.txt, avirgo_O4high_NEW.txt
                'design_H1L1V1' : AplusDesign.txt, AplusDesign.txt, avirgo_O5high_NEW.txt
                'infinite': intrinsic merging DCO population, i.e. p_det = 1
        path_to_dir : string
            Path to the directory where you the cosmological weights are stored.

        Returns
        -------
        array doubles
            Return the cosmological weights, z_formation, z_merger and binary
            index k associated to each weighted binary.

        """
        dir_ = os.path.join(path_to_dir, f'DCOs/{sensitivity}_sensitivity')

        if extention == 'npz':
            if self.verbose:
                print('Loading the data ...')
            index = np.load(os.path.join(dir_, 'index.npz'), allow_pickle=True)['key']
            z_formation = np.load(os.path.join(dir_, 'z_formation.npz'), allow_pickle=True)['key']
            z_merger = np.load(os.path.join(dir_, 'z_merger.npz'), allow_pickle=True)['key']
            weights = np.load(os.path.join(dir_, 'weights.npz'), allow_pickle=True)['key']
        else:
            raise ValueError('Extension not supported!')
        return index, z_formation, z_merger, weights


    def get_shell_comovig_volume(self, z_hor_i, z_hor_f, sensitivity='infinite'):
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
                        * (self.get_comoving_disntance_from_redshift(z)
                        * 10**(-3.))**2. / E(z))
            else:
                # TODO: peanut-shaped antenna patter comoving volume calculation
                raise ValueError('Sensitivity not supported!')
        return sp.integrate.quad(f, z_hor_i, z_hor_f, args=(sensitivity))[0] # Gpc^3


    def compute_rate_density(self, w_ijk, z_event, observable='DCO', sensitivity='infinite', index=None):
        """Compute the GRB/DCO rate density.

        Parameters
        ----------
        w_ijk : array doubles
            Cosmological weights computed with Eq. B.8 of Bavera et at. (2020).
        z_event : array doubles
            Cosmolgocial redshift of the event you are tracking.
        observable : string
            Event you are tracking, available:
            'DCOs': merger event of a DCO system
            'GRB1': gamma ray bursts of star 1
            'GRB2': gamma ray bursts of star 2
        sensitivity : string
            This takes into account the detector sensitivity, available:
            'infinite': p_det = 1
            'beamed': TODO for GRBs

        Returns
        -------
        array doubles
            Return the DCOs merger rate density (Gpc^-3 yr^-1) as a function of
            cosmolgocial redshift as in Eq. (D.1) in Bavera et al. (2022)
            arXiv:2106.15841

        """

        z_hor = self.get_edges_redshift_bins()
        n = len(z_hor)

        if observable=='DCO':
            z_merger_DCO = z_event
            Rate_DCO = np.zeros(n-1)
            if sensitivity=='infinite':
                for i in range(1,n):
                    # compute Eq. (D.1) in Bavera et al. (2022) arXiv:2106.15841
                    cond_DCO = np.logical_and(z_merger_DCO>z_hor[i-1],
                                                    z_merger_DCO<=z_hor[i])
                    Rate_DCO[i-1] = (sum(w_ijk[cond_DCO])
                                     /self.get_shell_comovig_volume(z_hor[i-1],
                                                                    z_hor[i],
                                                                    sensitivity))
                return Rate_DCO # Gpc^-3 yr^-1
            else:
                raise ValueError('Unsupported sensitivity!')
        
        elif observable in ['GRB1','GRB2']:
            s = observable[-1]
            z_GRB = z_event
            if sensitivity=='beamed':
                if index is not None:
                    f_beaming = self.get_data(f"S{s}_f_beaming",index)
                    flag_GRB = self.get_data(f"GRB{s}",index)
                else:
                    raise ValueError('Missing f_beaming parameter!')
                Rate_GRB = np.zeros(n-1)
                for i in range(1,n):
                    cond_GRB = np.logical_and(np.logical_and(z_GRB>z_hor[i-1], z_GRB<=z_hor[i]),
                                              flag_GRB)
                    f_fb = f_beaming[cond_GRB]
                    Rate_GRB[i-1] = (sum(w_ijk[cond_GRB]*f_fb)
                                      /self.get_shell_comovig_volume(z_hor[i-1],
                                                                     z_hor[i],
                                                                     sensitivity='infinite'))
                return Rate_GRB

            elif sensitivity=='infinite':
                if index is not None:
                    flag_GRB = self.get_data(f"GRB{s}",index)
                else:
                    raise ValueError('Missing f_beaming parameter!')
                Rate_GRB = np.zeros(n-1)
                for i in range(1,n):
                    if index is None:
                        raise ValueError('Provide index comlumn to identify GRB systems.')
                    cond_GRB = np.logical_and(np.logical_and(z_GRB>z_hor[i-1], z_GRB<=z_hor[i]),
                                              flag_GRB)
                    Rate_GRB[i-1] = (sum(w_ijk[cond_GRB])
                                      /self.get_shell_comovig_volume(z_hor[i-1],
                                                                     z_hor[i],
                                                                     sensitivity))
                return Rate_GRB
            else:
                raise ValueError('Unknown sensitivity!')
        else:
            raise ValueError('Unknown observable!')

    ##############################
    #####  GRB class methods #####
    ##############################
    
    def get_time_GRB(self, z_birth, event=None):
        """Get the time of the GRB.

        Parameters
        ----------
        z_brith : double
            Redshif of birth.
        event : string
            Event you are tracking, either first or second core collpase:
            'CC1', 'CC2'.

        Returns
        -------
        t_BRB : double
            Cosmic time of the GRB event in Gyr.

        """
        if event not in ['CC1', 'CC2']:
            raise ValueError(f'Unknown event {event}!')
        n = self.df.shape[0]
        t_birth = self.get_cosmic_time_from_redshift(z_birth) * np.ones(n)  # Gyr
        t_GRB = t_birth + self.df[f"time_{event}"] * 10 ** (-3) # Gyr
        return t_GRB

            
    def compute_GRB_rate_weights(self, sensitivity='infinity', path_to_dir='./', extention='npz'):
        """Compute the cosmological weights of the transient events associated to the population.

        This function will create a directory path_to_dir/GRBs/sensitivity where
        it will save the weigths, binary indicies k, z_formation, z_grb.
        This is needed for scalability as a ~100k DCO population generates ~10M
        non-zero weights assuming a delta_t=100Myr.

        Parameters
        ----------
        sensitivity : string
            Assume there are no selection effects. Available:
            'infinite': whole GRB population, i.e. p_det = 1
        path_to_dir : string
            Path to the workingn directory where you want to store the
            cosmological weights.

        """
        
        # check if the folder three exists, otherwise create it
        GRB_dir = os.path.join(path_to_dir,'GRBs')
        if 'GRBs' not in os.listdir(path_to_dir):
            os.makedirs(GRB_dir)

        sensitivity_dir = os.path.join(GRB_dir, f'{sensitivity}_sensitivity')
        if f'{sensitivity}_sensitivity' not in os.listdir(GRB_dir):
            os.makedirs(sensitivity_dir)

        # index of all DCOs
        n = self.df.shape[0]
        index = np.arange(0, n)

        # redshif of each time bin
        z_birth = self.get_centers_redshift_bins()

        # define interpolator
        get_redshift_from_time = self.get_redshift_from_cosmic_time_interpolator()

        # hashmap to store everythig
        data = {'index_1': {}, 'z_grb_1': {}, 'weights_1': {},
                'index_2': {}, 'z_grb_2': {}, 'weights_2': {}}

        # loop over all redshift bins
        for i in tqdm(range(len(z_birth))):
            
            for s, event in enumerate(['CC1', 'CC2']):
                # compute the CC time of each compact object
                t_CC = self.get_time_GRB(z_birth[i], event)

                # intrinsic population
                if sensitivity == 'infinite':
                    
                    # sort out system not emitting GRBs
                    bool_GRB = np.logical_and(t_CC < cosmology.age(1e-08).value*0.9999999 , self.df[f"GRB{s+1}"])
                    
                    # if there are no system emitting any GRB, continue
                    if len(index[bool_GRB]) == 0:
                        data[f'index_{s+1}'][str(z_birth[i])] = np.array([])
                        data[f'z_grb_{s+1}'][str(z_birth[i])] = np.array([])
                        data[f'weights_{s+1}'][str(z_birth[i])] = np.array([])
                        continue
                
                    # get GRB redshifts
                    z_GRB = get_redshift_from_time(t_CC[bool_GRB])

                    # set detection probabilities to 1
                    p_det = np.ones(len(index[bool_GRB]))

                    # store index and redshift of merger
                    data[f'index_{s+1}'][str(z_birth[i])] = index[bool_GRB]
                    data[f'z_grb_{s+1}'][str(z_birth[i])] = z_GRB

                    # compute and store cosmological rate weights as in eq. B.8 in Bavera et al. (2020)
                    z_b = np.ones(len(index[bool_GRB]))*z_birth[i]
                    w_ijk = self.merger_rate_weight(z_b, z_GRB, p_det, index[bool_GRB])
                    data[f'weights_{s+1}'][str(z_birth[i])] = w_ijk

                else:
                    raise ValueError('Unknown sensitivity!')

        # TODO: implemet the saving to h5 files
        if extention == 'npz':
            if self.verbose:
                print('Formatting the data ....')
            for s in [1,2]:
                data_to_save = [[],[],[],[]]
                for i, dict in enumerate([data[key] for key in data.keys() if f'{s}' in key]):
                    for key in dict.keys():
                        data_to_save[i].extend(dict[key].tolist())
                        if i == 0:
                            data_to_save[-1].extend(np.ones(len(dict[key]))*float(key))
                if self.verbose:
                    print('Saving the data ....')
                for i, key in enumerate([key for key in data.keys() if f'{s}' in key]):
                    if key == 'index':
                        fmt_str = '%i'
                    else:
                        fmt_str = '%.8E'
                    np.savez(os.path.join(sensitivity_dir, f"{key}.npz"),
                            key=data_to_save[i], fmt=fmt_str)

                np.savez(os.path.join(sensitivity_dir,f"z_formation_{s}.npz"),
                        key=data_to_save[-1], fmt='%.8E')
        else:
            raise ValueError('Extension not supported!')


    def load_grb_rate_weights(self, sensitivity, path_to_dir='./', extention='npz'):
        """Load the cosmological weights of the transient events associated to the population.

        Parameters
        ----------
        sensitivity : string
            Assume there are no selection effects. Available:
            'infinite': whole GRB population, i.e. p_det = 1
        path_to_dir : string
            Path to the directory where you the cosmological weights are stored.

        Returns
        -------
        array doubles
            Return the cosmological weights, z_formation, z_GRB and binary
            index k associated to each weighted binary.

        """
        dir_ = os.path.join(path_to_dir, f'GRBs/{sensitivity}_sensitivity')

        if extention == 'npz':
            if self.verbose:
                print('Loading the data ...')
            index_1 = np.load(os.path.join(dir_, 'index_1.npz'), allow_pickle=True)['key']
            z_formation_1 = np.load(os.path.join(dir_, 'z_formation_1.npz'), allow_pickle=True)['key']
            z_grb_1 = np.load(os.path.join(dir_, 'z_grb_1.npz'), allow_pickle=True)['key']
            weights_1 = np.load(os.path.join(dir_, 'weights_1.npz'), allow_pickle=True)['key']
            index_2 = np.load(os.path.join(dir_, 'index_2.npz'), allow_pickle=True)['key']
            z_formation_2 = np.load(os.path.join(dir_, 'z_formation_2.npz'), allow_pickle=True)['key']
            z_grb_2 = np.load(os.path.join(dir_, 'z_grb_2.npz'), allow_pickle=True)['key']
            weights_2 = np.load(os.path.join(dir_, 'weights_2.npz'), allow_pickle=True)['key']
        else:
            raise ValueError('Extension not supported!')
        return index_1, z_formation_1, z_grb_1, weights_1, index_2, z_formation_2, z_grb_2, weights_2