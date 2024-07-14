import numpy as np
import pandas as pd
import posydon.popsyn.selection_effects as selection_effects
from posydon.config import PATH_TO_POSYDON_DATA
import os
from tqdm import tqdm


import warnings
# This is to suppress the performance warnings from pandas
# These warnings are not important for the user
# We should alter the code to remove these warnings, but for now we suppress them
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

PATH_TO_PDET_GRID = os.path.join(PATH_TO_POSYDON_DATA, 'POSYDON_data/selection_effects/pdet_grid.hdf5')

def GRB_selection(history_chunk, oneline_chunk, formation_channels_chunk=None, S1_S2='S1'):
    """A GRB selection function to create a transient population of LGRBs.
    
    This function requires a wrapper function to be used to create a transient population, because
    of the extra parameters that are needed to be passed to the function (S1_S2).
    
    >>> def GRB_selection_wrapper(history_chunk, oneline_chunk, formation_channels_chunk=None):
    ...     return GRB_selection(history_chunk, oneline_chunk, formation_channels_chunk, S1_S2='S1')
        
    This is an example function for selecting LGRBs and some of their properties.
    Additional properties from the history or oneline can be added to the transient population.
    The output dataframe contains the following columns:
    - time : the time of the event
    - metallicity : the metallicity of the event
    - S1_m_disk_radiated : the mass of the disk radiated by the first star
    - S2_m_disk_radiated : the mass of the disk radiated by the second star
    
    Columns with pre and post SN data:
    - S1_mass : the mass of the first star
    - S2_mass : the mass of the second star
    - S1_spin : the spin of the first star
    - S2_spin : the spin of the second star
    - orbital_period : the orbital period of the binary
    - eccentricity : the eccentricity of the binary

    Parameters
    ----------
    history_chunk : pd.DataFrame
        The history chunk of the population.
    oneline_chunk : pd.DataFrame
        The oneline chunk of the population.
    formation_channels_chunk : pd.DataFrame (can be made optional)
        The formation pathway chunk of the population.
    S1_S2 : str
        The star that will be the progenitor of the GRB. Can be either 'S1' or 'S2'.
    
    Returns
    -------
    pd.DataFrame
        A DataFrame containing the transient population of LGRBs with the columns mentioned above.
    
    Example
    -------
    
    >>> GRB_selection_wrapper(pop.history[0], pop.oneline[0], pop.formation_channels.loc[0])
    >>> pop.create_transient_population(GRB_selection_wrapper, name='GRB_S1')
    """
    columns_pre_post  = ['orbital_period', 'eccentricity', 'S1_spin', 'S2_spin']
    columns = []
    oneline_columns = ['metallicity', 'S1_m_disk_radiated', 'S2_m_disk_radiated']
   
    if S1_S2 == 'S1':
        indices_selection = oneline_chunk.index[oneline_chunk['S1_m_disk_radiated'] > 0.0].to_numpy()
        oneline_chunk = oneline_chunk.drop(columns=['S2_m_disk_radiated'])
    elif S1_S2 == 'S2':
        indices_selection = oneline_chunk.index[oneline_chunk['S2_m_disk_radiated'] > 0.0].to_numpy()
        oneline_chunk = oneline_chunk.drop(columns=['S1_m_disk_radiated'])
    else:
        raise ValueError('S1_S2 must be either S1 or S2')
    # no events in this chunk
    if len(indices_selection) == 0:
        return pd.DataFrame()
    # filter out the events that are not relevant for the LGRB formation
    selection = history_chunk.loc[indices_selection]    
    if S1_S2 == 'S1':
        S_mask = (selection['S1_state'] == 'BH') & (selection['S1_state'] != 'BH').shift(1) & (selection['step_names'] == 'step_SN')
    elif S1_S2 == 'S2':
        S_mask = (selection['S2_state'] == 'BH') & (selection['S2_state'] != 'BH').shift(1) & (selection['step_names'] == 'step_SN')
    
    GRB_df_synthetic = pd.DataFrame(index=indices_selection)
    # Pre and post SN data
    post_SN_hist = selection[S_mask]
    pre_mask = S_mask.shift(-1, fill_value=False)
    pre_SN_hist = selection[pre_mask]
    
    if S1_S2 == 'S1':
        columns_pre_post.append('S1_mass')
        columns.append('S2_mass')
    elif S1_S2 == 'S2':
        columns_pre_post.append('S2_mass')
        columns.append('S1_mass')
        
    
    for col in columns_pre_post:
        GRB_df_synthetic[col+'_preSN'] = pre_SN_hist[col].values
        GRB_df_synthetic[col+'_postSN'] = post_SN_hist[col].values
        
    # unchanged columns
    for col in columns:
        GRB_df_synthetic[col] = post_SN_hist[col].values
        
    # oneline data
    for col in oneline_chunk.columns:
        GRB_df_synthetic[col] = oneline_chunk.loc[indices_selection][col].values
    
    if any(formation_channels_chunk != None):
        formation_channels_chunk = formation_channels_chunk.loc[indices_selection]
        if S1_S2 == 'S1':
            GRB_df_synthetic['channel'] = formation_channels_chunk['channel'].str.split('_CC1').str[0].apply(lambda x: x+'_CC1')
        elif S1_S2 == 'S2':
            GRB_df_synthetic['channel'] = formation_channels_chunk['channel'].str.split('_CC2').str[0].apply(lambda x: x+'_CC2')
    
    # calculate the time!
    GRB_df_synthetic['time'] = post_SN_hist['time'].values * 1e-6 # convert to Myr
        
    return GRB_df_synthetic


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

def BBH_selection_function(history_chunk, oneline_chunk, formation_channels_chunk=None):
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
    formation_channels_chunk : pd.DataFrame (can be made optional)
        The formation pathway chunk of the DCO population.
        
    Returns
    -------
    df_transients : pd.DataFrame
        A DataFrame containing the transient population of BBHs.
        This DataFrame contains the following columns:
        - time : the time of the event
        - metallicity : the metallicity of the event
    
    '''
    
    indices = oneline_chunk.index.to_numpy()
    df_transients = pd.DataFrame(index = indices)
    
    df_transients['time'] = history_chunk[history_chunk['event'] == 'END']['time'] * 1e-6 #Myr
    mask = (history_chunk['S1_state'] == 'BH') & (history_chunk['S2_state'] == 'BH') & (history_chunk['step_names'] == 'step_SN') & (history_chunk['state'] == 'detached')

    df_transients['t_inspiral'] = df_transients['time'] - history_chunk[mask]['time']*1e-6
    df_transients['metallicity'] = oneline_chunk['metallicity']
    df_transients['S1_state']  = history_chunk[mask]['S1_state']
    df_transients['S2_state']  = history_chunk[mask]['S2_state']
    df_transients['S1_mass'] = history_chunk[mask]['S1_mass']
    df_transients['S2_mass'] = history_chunk[mask]['S2_mass']
    df_transients['S1_spin'] = history_chunk[mask]['S1_spin']
    df_transients['S2_spin'] = history_chunk[mask]['S2_spin']
    df_transients['S1_spin_orbit_tilt'] = oneline_chunk['S1_spin_orbit_tilt']
    df_transients['S2_spin_orbit_tilt'] = oneline_chunk['S2_spin_orbit_tilt']
    df_transients['orbital_period'] = history_chunk[mask]['orbital_period']
    df_transients['chirp_mass'] = m_chirp(history_chunk[mask]['S1_mass'], history_chunk[mask]['S2_mass'])
    df_transients['mass_ratio'] = mass_ratio(history_chunk[mask]['S1_mass'], history_chunk[mask]['S2_mass'])
    df_transients['chi_eff'] = chi_eff(history_chunk[mask]['S1_mass'], history_chunk[mask]['S2_mass'], history_chunk[mask]['S1_spin'], history_chunk[mask]['S2_spin'], oneline_chunk['S1_spin_orbit_tilt'], oneline_chunk['S2_spin_orbit_tilt'])
    df_transients['eccentricity'] = history_chunk[mask]['eccentricity']

    if formation_channels_chunk is not None:
        df_transients = pd.concat([df_transients, formation_channels_chunk[['channel']]], axis=1)
    
    return df_transients
    

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
        if np.sum(mask) == 0:
            detectable_weights[mask, i] = 0.0
        else:
            detectable_weights[mask, i] = detectable_weights[mask, i] * sel_eff.predict_pdet(data_slice[mask])
        
    return pd.DataFrame(detectable_weights, index=z_events_chunk.index, columns=z_events_chunk.columns)
    
    
