import numpy as np
import pandas as pd
import posydon.popsyn.selection_effects as selection_effects
from posydon.config import PATH_TO_POSYDON_DATA
import os
from tqdm import tqdm

PATH_TO_PDET_GRID = os.path.join(PATH_TO_POSYDON_DATA, 'POSYDON_data/selection_effects/pdet_grid.hdf5')

# RAW FUNCTIONS! NEED TO BE EDITED TO WORKS WITH THE NEW STRUCTURE!!!!

def create_GRB_population(output_file=None, GRB_properties={}, additional_oneline_cols=None):
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


def chi_eff(m_1, m_2, a_1, a_2, tilt_1, tilt_2):
    '''Calculate the effective spin of two masses.'''
    return (m_1*a_1*np.cos(tilt_1)+m_2*a_2*np.cos(tilt_2))/(m_1+m_2)

def m_chirp(m_1, m_2):
    '''Calculate the chirp mass of two masses.'''
    return (m_1*m_2)**(3./5)/(m_1+m_2)**(1./5)

def mass_ratio(m_1, m_2):
    '''Calculate the mass ratio of two masses.'''
    q = m_2/m_1
    q[q>1.] = 1./q[q>1.]
    return q

def BBH_selection_function(history_chunk, oneline_chunk, formation_pathway_chunks):
    '''A BBH selection function to create a transient population of BBHs mergers.
    
    This is an example function for selecting BBH mergers and some of their properties.
    Additional properties from the history or oneline can be added to the transient population.
    The output dataframe contains the following columns:
    - time : the time of the event
    - metallicity : the metallicity of the event
    All other columns are optional.
    
    Parameters
    ----------
    history_chunk : pd.DataFrame
        The history chunk of the DCO population.
    oneline_chunk : pd.DataFrame
        The oneline chunk of the DCO population.
    formation_pathway_chunks : pd.DataFrame (can be made optional)
        The formation pathway chunk of the DCO population.
        
    Returns
    -------
    df_synthetic : pd.DataFrame
        A DataFrame containing the transient population of BBHs.
        This DataFrame contains the following columns:
        - time : the time of the event
        - metallicity : the metallicity of the event
    
    '''
    
    indices = oneline_chunk.index.to_numpy()
    df_synthetic = pd.DataFrame(index = indices)
    
    df_synthetic['time'] = history_chunk[history_chunk['event'] == 'END']['time'] * 1e-6 #Myr
    mask = (history_chunk['S1_state'] == 'BH') & (history_chunk['S2_state'] == 'BH') & (history_chunk['step_names'] == 'step_SN') & (history_chunk['state'] == 'detached')

    df_synthetic['t_inspiral'] = df_synthetic['time'] - history_chunk[mask]['time']*1e-6
    
    df_synthetic['metallicity'] = oneline_chunk['metallicity']/0.0142
    df_synthetic['S1_state']  = history_chunk[mask]['S1_state']
    df_synthetic['S2_state']  = history_chunk[mask]['S2_state']
    df_synthetic['S1_mass'] = history_chunk[mask]['S1_mass']
    df_synthetic['S2_mass'] = history_chunk[mask]['S2_mass']
    df_synthetic['S1_spin'] = history_chunk[mask]['S1_spin']
    df_synthetic['S2_spin'] = history_chunk[mask]['S2_spin']
    df_synthetic['S1_spin_orbit_tilt'] = oneline_chunk['S1_spin_orbit_tilt']
    df_synthetic['S2_spin_orbit_tilt'] = oneline_chunk['S2_spin_orbit_tilt']
    df_synthetic['orbital_period'] = history_chunk[mask]['orbital_period']
    df_synthetic['chirp_mass'] = m_chirp(history_chunk[mask]['S1_mass'], history_chunk[mask]['S2_mass'])
    df_synthetic['mass_ratio'] = mass_ratio(history_chunk[mask]['S1_mass'], history_chunk[mask]['S2_mass'])
    df_synthetic['chi_eff'] = chi_eff(history_chunk[mask]['S1_mass'], history_chunk[mask]['S2_mass'], history_chunk[mask]['S1_spin'], history_chunk[mask]['S2_spin'], oneline_chunk['S1_spin_orbit_tilt'], oneline_chunk['S2_spin_orbit_tilt'])
    df_synthetic['eccentricity'] = history_chunk[mask]['eccentricity']

    df_synthetic = pd.concat([df_synthetic, formation_pathway_chunks[['channel']]], axis=1)
    
    return df_synthetic
    

def DCO_detactability(sensitivity, transient_pop_chunk, z_events_chunk, z_weights_chunk, verbose=False):
    '''Calculate the observability of a DCO population.
    
    Parameters
    ----------
    sensitivity : float
        The sensitivity of the detector.
        Available sensitivities are (for Hanford, Livingston, Virgo network):
        - 'O3actual_H1L1V1' : aligo_O3actual_H1.txt, aligo_O3actual_L1.txt, avirgo_O3actual.txt
        - 'O4low_H1L1V1' : aligo_O4low.txt, aligo_O4low.txt, avirgo_O4high_NEW.txt
        - 'O4high_H1L1V1' : aligo_O4high.txt, aligo_O4high.txt, avirgo_O4high_NEW.txt
        - 'design_H1L1V1' : AplusDesign.txt, AplusDesign.txt, avirgo_O5high_NEW.txt
        
        GW detector sensitivity and network configuration you want to use, see arXiv:1304.0670v3
        detector sensitivities are taken from: https://dcc.ligo.org/LIGO-T2000012-v2/public
    
    '''
    available_sensitiveies = ['O3actual_H1L1V1', 'O4low_H1L1V1', 'O4high_H1L1V1', 'design_H1L1V1']
    if sensitivity not in available_sensitiveies:
        raise ValueError(f'Unknown sensitivity {sensitivity}. Available sensitivities are {available_sensitiveies}')
    else:
        sel_eff = selection_effects.KNNmodel(grid_path=PATH_TO_PDET_GRID,
                                                      sensitivity_key=sensitivity)
    
    data_slice = pd.DataFrame()
    data_slice['m1'] = transient_pop_chunk['S1_mass']
    
    if 'q' in transient_pop_chunk.columns:
        data_slice['q'] = transient_pop_chunk['q']
    else:
        data_slice['q'] = mass_ratio(transient_pop_chunk['S1_mass'].to_numpy(),  transient_pop_chunk['S2_mass'].to_numpy())
    
    if 'chi_eff' in transient_pop_chunk.columns:
        data_slice['chieff'] = transient_pop_chunk['chi_eff']
        
    else:
        data_slice['chieff'] = chi_eff(transient_pop_chunk['S1_mass'],
                                        transient_pop_chunk['S2_mass'],
                                        transient_pop_chunk['S1_spin'],
                                        transient_pop_chunk['S2_spin'],
                                        transient_pop_chunk['S1_spin_orbit_tilt'],
                                        transient_pop_chunk['S2_spin_orbit_tilt'])
    
    detectable_weights = z_weights_chunk.to_numpy()
    for i in tqdm(range(z_events_chunk.shape[1]), total=z_events_chunk.shape[1], disable= not verbose):
        data_slice['z'] = z_events_chunk.iloc[:,i]
        mask = ~np.isnan(data_slice['z']).to_numpy()
        detectable_weights[mask, i] = detectable_weights[mask, i] * sel_eff.predict_pdet(data_slice[mask])
        
    return pd.DataFrame(detectable_weights, index=z_events_chunk.index, columns=z_events_chunk.columns)
    
    