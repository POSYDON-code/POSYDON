__authors__ = [
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Max Briel <max.briel@unige.ch>",
    ]

import os
import numpy as np
import matplotlib as mpl
import pandas as pd
from posydon.config import PATH_TO_POSYDON, PATH_TO_POSYDON_DATA
from posydon.visualization.plot_defaults import DEFAULT_LABELS
from posydon.utils.constants import Zsun
from posydon.grids.psygrid import PSyGrid
from posydon.utils.common_functions import convert_metallicity_to_string
from posydon.utils.constants import Zsun
import matplotlib.pyplot as plt
from posydon.visualization.plot_defaults import DEFAULT_LABELS

plt.style.use(os.path.join(PATH_TO_POSYDON, "posydon/visualization/posydon.mplstyle"))

cm = mpl.colormaps.get_cmap('tab20')
COLORS = [cm.colors[i] for i in range(len(cm.colors)) if i%2==0] + [cm.colors[i] for i in range(len(cm.colors)) if i%2==1]


def plot_perley16_rate_density(ax=None):
    '''Plot the GRB rate density from Perley et al. (2016)
    
    plot the long gamma-ray burst rate density from Perley et al. (2016)
    if ax is None, it plotted on the latest axis, otherwise it is added to the given axis
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes (optional)
        Axes object to plot the data on
    '''
    # Perley et al. (2016)
    
    z_P16 = np.array([0.3,0.75,1.35,1.75,2.3,3.1,4,5.2,7.])
    rate_P16 = np.array([2.07e-10,1.46e-9,1.816e-9,3.93e-9,8.27e-9,6.90e-9,4.20e-9,3.74e-9,1.11e-9])*1e9
    rate_P16_lower = rate_P16-np.array([1.12e-10,1.08e-9,1.30e-9,3.02e-9,6.78e-9,5.29e-9,2.82e-9,2.397e-9,5.96e-10])*1e9
    rate_P16_upper = np.array([5.40e-10,2.03e-9,2.60e-9,5.29e-9,1.03e-8,9.13e-9,6.78e-9,6.78e-9,4.27e-9])*1e9-rate_P16
    rate_P16_error_y =np.array([rate_P16_lower.tolist(),rate_P16_upper.tolist()])
    rate_P16_left = z_P16-np.array([0.1,0.5,1,1.5,2.,2.6,3.5,4.5,6])
    rate_P16_right = np.array([0.5,1,1.5,2.,2.6,3.5,4.5,6,8])-z_P16
    rate_P16_error_x =np.array([rate_P16_left.tolist(),rate_P16_right.tolist()])
    if ax is None:
        plt.errorbar(z_P16,rate_P16,xerr=rate_P16_error_x,yerr=rate_P16_error_y, fmt='.', color='black',
                    label=r'$\mathcal{R}_\mathrm{LGRB}(E_\mathrm{iso}^{45-450\,\mathrm{keV}} > 10^{51} \, \mathrm{erg})$ SHOALS survey')
    else:
        ax.errorbar(z_P16,rate_P16,xerr=rate_P16_error_x,yerr=rate_P16_error_y, fmt='.', color='black',
                    label=r'$\mathcal{R}_\mathrm{LGRB}(E_\mathrm{iso}^{45-450\,\mathrm{keV}} > 10^{51} \, \mathrm{erg})$ SHOALS survey')

def plot_rate_density(intrinsic_rates, channels=False, **kwargs):
    plt.figure()

    plt.plot(intrinsic_rates.index, intrinsic_rates['total'], label='total', color='black')

    if channels:
        for i, ch in enumerate([key for key in intrinsic_rates.keys() if key != 'total']):
            if ch != 'total':
                plt.plot(intrinsic_rates.index, intrinsic_rates[ch], label=ch, color=COLORS[i])
    plt.yscale('log')
    plt.ylabel(r'$\mathcal{R} \,[\mathrm{Gpc}^{-3}\,\mathrm{yr}^{-1}]$')
    plt.xlabel(r'$z$')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    if 'ylim' in kwargs:
        plt.ylim(kwargs['ylim'])
    if 'xlim' in kwargs:
        plt.xlim(kwargs['xlim'])
    if 'path' in kwargs and isinstance(kwargs['path'],str):
        plt.savefig(kwargs['path'])
    if 'show' in kwargs and kwargs['show']:
        plt.show()

def plot_merger_efficiency(met, merger_efficiency, show=True, path=None, channels=False):
    title = r'Event efficiency'
    plt.figure()
    plt.title(title)
    plt.plot(met/Zsun, merger_efficiency['total'], label='total', color='black')
    if channels:
        unique_channels = [key for key in merger_efficiency.keys() if key != 'total']
        if len(unique_channels) > len(cm.colors):
            raise ValueError('Too many channels to plot!')
        for i, ch in enumerate(unique_channels):
            if ch != 'total':
                plt.plot(met/Zsun, merger_efficiency[ch], label=ch, color=COLORS[i])
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel(r'$Z/Z_\odot$')
    plt.ylabel('\#events [$M_\odot^{-1}$]')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    if path:
        plt.savefig(path)
    if show:
        plt.show()


def plot_hist_properties(df, ax=None, df_intrinsic=None, df_observable=None,
                        channel=None,
                        show=True, path=None, alpha=0.5,
                        range=None, bins=20, normalise=False,color=COLORS[0], label='', **kwargs):
    if ax is None:
        fig, ax = plt.subplots(1,1, figsize=(5,5))

    df_columns = df.columns
    if 'intrinsic' in df_columns and 'observable' in df_columns:
        title = r'Intrinsic vs. observable (dashed) population'
    elif 'intrinsic' in df_columns:
        title = r'Intrinsic population'
    elif 'observable' in df_columns:
        title = r'Observable population'
    else:
        raise ValueError('You should provide either an intrinsic or a '
                         'detectable population.')

    if 'intrinsic' in df_columns:
        values =  df['intrinsic']
        if normalise:
            values /= sum(values)

        ax.hist(df['property'],
                 weights=values,
                 color=color,
                 alpha=alpha,
                 range=range,
                 bins=bins,
                 label=label+' intrinsic')

    if 'observable' in df_columns:
        values = df['observable']
        if normalise:
            values /= sum(values)

        ax.hist(df['property'],
                 weights=values,
                 color=color,
                 alpha=alpha,
                 range=range,
                 histtype=u'step',
                 linestyle='--',
                 bins=bins,
                 label=label+' observable')

    ax.set_title(title)

    if 'xlabel' in kwargs:
        ax.set_xlabel(kwargs['xlabel'])
    if normalise:
        ax.set_ylabel(r'PDF')
    else:
        ax.set_ylabel(r'\#events in bin')

    if 'yscale' in kwargs:
        ax.set_yscale(kwargs['yscale'])
    if 'xscale' in kwargs:
        ax.set_xscale(kwargs['xscale'])

    if path:
        plt.savefig(path)
    if show:
        plt.show()


def plot_popsyn_over_grid_slice(pop, grid_type, met_Zsun,
                                termination_flag='combined_TF12',
                                slices=None, channel=None,
                                plot_dir='./', prop=None, prop_range=None,
                                log_prop=False, alpha=0.3, s=5.,
                                show_fig=True, save_fig=True, close_fig=True,
                                plot_extension='png', verbose=False):
    '''Plot the population synthesis data over a grid slice
    
    Outputs one or multiple plots of the population synthesis data over a grid slice.
    Either stores the plots in a directory or shows them.
    
    Parameters
    ----------
    pop : Population
        Population object
    grid_type : str
        Grid type to get data for.
        Options are 'HMS-HMS', 'CO-HMS_RLO', 'CO-HeMS', 'CO-HeMS_RLO'.
    met_Zsun : float
        Metallicity of the population in solar units
    termination_flag : str
        Termination flag for the grid
    slices : list (optional)
        List of slices to plot
    channel : str (optional)
        Formation channel to be plotted
    plot_dir : str (optional)
        Directory to save the plots
    prop : str (optional)
        Property to plot on top of the grid slice
    prop_range : list (optional)
        Range of the property to plot
    log_prop : bool (optional)
        Logarithmic scale for the property
    alpha : float (optional)
        Transparency of the data points
    s : float (optional)
        Size of the data points
    show_fig : bool (optional)
        Show the figure
    save_fig : bool (optional)
        Save the figure
    close_fig : bool (optional)
        Close the figure
    plot_extension : str (optional)
        File extension for saving the plot
    verbose : bool (optional)
        Print information
    '''

    # Check if step_names in pop.history data
    if 'step_names' not in pop.history.columns:
        raise ValueError('Formation channel information not available in popsynth data.')
    # check if formation channel information is avaialbe

    if channel is not None and 'channel' not in pop.population.columns:
        raise ValueError('Formation channel information not available in popsynth data.')

    # get population data from the popsyn for the given metallicity
    data = get_population_data(pop=pop, metallicity=met_Zsun,
                               grid_type=grid_type, channel=channel, prop=prop)
    # Setup the grid 
    met = convert_metallicity_to_string(met_Zsun)
    grid_path = os.path.join(PATH_TO_POSYDON_DATA, f'{grid_type}/{met}_Zsun.h5')
    grid = PSyGrid(verbose=False)
    grid.load(grid_path)
    bin_centers,bin_edges, PLOT_PROPERTIES, tmp_fname, slice_3D_var_str = \
                            setup_grid_slice_plotting(grid, grid_type, plot_extension)

    for i, bin_center in enumerate(bin_centers):
        # skip slices not in the list
        if slices is not None and round(bin_center,2) not in np.around(slices,2):
            if verbose:
                print(f'Skipping {round(bin_center,2)}')
            continue
        
        slice_3D_var_range = (bin_edges[i], bin_edges[i+1])
        
        # add additional information about the plot_slice to the title
        if slice_3D_var_str == 'mass_ratio':
            PLOT_PROPERTIES['title'] = f'q = {bin_center:.2f}'
        elif slice_3D_var_str == 'star_2_mass':
            PLOT_PROPERTIES['title'] = '$M_\mathrm{CO}'+f' = {bin_center:.2f} M_{{\odot}}$'
        # Check for channel
        if channel is not None:
            PLOT_PROPERTIES['title'] += f'\n {channel}'
        
        # change the file name for the plot
        PLOT_PROPERTIES['fname'] = tmp_fname % bin_center
        
        # change size of the figure for properties
        if prop is not None:
            PLOT_PROPERTIES['figsize'] = (4,5.5)
        
        # plot the grid slice
        plot_grid_slice(grid,
                        slice_3D_var_str,
                        slice_3D_var_range,
                        termination_flag=termination_flag, 
                        PLOT_PROPERTIES=PLOT_PROPERTIES)
        
        # plot data on top of the grid slice
        plot_population_data(data, slice_3D_var_str, slice_3D_var_range, 
                             prop, prop_range, log_prop, alpha, s)

        if save_fig:
            plt.savefig(os.path.join(plot_dir, PLOT_PROPERTIES['fname']), bbox_inches='tight')
        if show_fig:
            plt.show()
        if close_fig:
            plt.close()

def get_population_data(pop, metallicity, grid_type, channel=None, prop=None):
    '''get the population data of a given metallicity for plotting
    
    Parameters
    ----------
    pop : Population
        Population object
    metallicity : float
        Metallicity of the population in solar units
    grid_type : str
        Grid type to get data for.
        Options are 'HMS-HMS', 'CO-HMS_RLO', 'CO-HeMS', 'CO-HeMS_RLO'.
    channel : str
        Formation channel to be plotted
        
    Returns
    -------
    data : DataFrame
        Population data of the given metallicity and channel and grid type
    '''
    
    # 1. mask channel
    if channel is not None:
        channel_mask = pop.population['channel'] == channel
    else:
        channel_mask = True
    
    # 2. mask metallicity
    metallicity_mask = pop.population['metallicity'] == metallicity
    
    # 3. Get the history data for the selected population based on the grid/step names
    selected_indices = pop.population.index[channel_mask & metallicity_mask].tolist()
    
    where_string = "(index in "+str(selected_indices)+")"
    
    data = pop.history.select(where=where_string,
                              columns=['S1_mass',
                                       'S2_mass',
                                       'orbital_period',
                                       'event',
                                       'step_names']).copy()
    
    # select the right data based on the grid type
    if grid_type == 'HMS-HMS':
        data = data[data['event'] == 'ZAMS']
    elif grid_type == 'CO-HMS_RLO':
        data = data[data['event'] == 'oRLO2']
        # swap the masses since data has CO as the primary
        tmp_S1 = data['S2_mass'].copy()
        data['S2_mass'] = data['S1_mass'].copy()
        data['S1_mass'] = tmp_S1
        
    elif grid_type == 'CO-HeMS':
        data = data[data['event'] == 'step_CE']
        # swap the masses since data has CO as the primary
        tmp_S1 = data['S2_mass'].copy()
        data['S2_mass'] = data['S1_mass'].copy()
        data['S1_mass'] = tmp_S1
    else:
        print('grid type not recognized')

    data['q'] = data['S2_mass']/data['S1_mass']

    # only relevant for Compact Object (CO) grids
    data.loc[data['q']>1, 'q'] = 1./data['q'][data['q']>1]
    
    # add prop of DCO merger to the data
    if prop is not None:
        if prop not in pop.population.columns:
            raise ValueError(f'Property {prop} not available in pop.population.')
        data[prop] = pop.population[prop][channel_mask & metallicity_mask].values
    
    # remove unnecessary columns
    data.drop(columns=['event', 'step_names'], inplace=True)
    return data

def setup_grid_slice_plotting(grid, grid_type, plot_extension):
    '''Setup the values for plotting the grid slice
    
    Parameters
    ----------
    grid : PSyGrid
        PSyGrid object
    grid_type : str
        Grid type to get data for.
        Options are 'HMS-HMS', 'CO-HMS_RLO', 'CO-HeMS', 'CO-HeMS_RLO'.
    plot_extension : str
        File extension for saving the plot
    
    Returns
    -------
    bin_centers : list
        List of bin centers
    bin_edges : list
        List of bin edges
    PLOT_PROPERTIES : dict
        Dictionary of plot properties
    slice_3D_var_str : str
        String of the 3D variable to slice the grid
    tmp_fname : str
        Template file name for saving the plot
    '''
    
    if grid_type == "HMS-HMS":
        grid_q_unique = np.unique(np.around(grid.initial_values['star_2_mass']/grid.initial_values['star_1_mass'],2))
        delta_q = (grid_q_unique[1:]-grid_q_unique[:-1])*0.5
        q_edges = (grid_q_unique[:-1]+delta_q).tolist()
        
        # set output variables
        bin_edges = [0.]+q_edges+[1.]
        bin_centers = grid_q_unique.tolist()
        tmp_fname = 'grid_q_%1.2f.' + plot_extension
        slice_3D_var_str='mass_ratio'
        
    elif 'CO' in grid_type:
        m_COs = np.unique(np.around(grid.initial_values['star_2_mass'],1))
        m_COs_edges = 10**((np.log10(np.array(m_COs)[1:])+np.log10(np.array(m_COs)[:-1]))*0.5)
        m2 = [0.]+m_COs_edges.tolist()+[2*m_COs_edges[-1]]
        
        bin_edges = m2
        bin_centers = m_COs
        tmp_fname = 'grid_m_%1.2f.' + plot_extension
        slice_3D_var_str='star_2_mass'
    else:
        print('grid type not recognized')
        
    PLOT_PROPERTIES = {
        'show_fig' : False,
        'close_fig' : False,
        'log10_x' : True,
        'log10_y' : True,
        'legend2D' : {'bbox_to_anchor' : (1.03, 0.5)},
    }

    return bin_centers, bin_edges, PLOT_PROPERTIES, tmp_fname, slice_3D_var_str

def plot_population_data(data,
                         slice_3D_var_str,
                         slice_3D_var_range,
                         prop=None,
                         prop_range=None,
                         log_prop=False,
                         alpha=0.3,
                         s=5):
    '''Plot the population data based on the grid slice parameters
    
    Parameters
    ----------
    data : DataFrame
        Population data
    slice_3D_var_str : str
        String of the 3D variable to slice the grid
    slice_3D_var_range : list
        Range of the 3D variable to slice the grid
    prop : str (optional)
        Property to plot on top of the grid slice
    prop_range : list (optional)
        Range of the property to plot
    log_prop : bool (optional)
        Logarithmic scale for the property
    alpha : float (optional)
        Transparency of the data points
    s : float (optional)
        Size of the data points
    '''

    # get only slice data
    if slice_3D_var_str == 'mass_ratio':
        var = 'q'
    elif slice_3D_var_str == 'star_2_mass':
        var = 'S2_mass'
    else:
        raise ValueError(f"Unknown slice_3D_var_str: {slice_3D_var_str}")
    
    # make it exclusive for the upper limit
    mask = (data[var] >= slice_3D_var_range[0]) & (data[var] < slice_3D_var_range[1])
    
    if prop is not None:
        if prop not in data.columns:
            raise ValueError(f'Property {prop} not available in popsynth data.')
        if prop_range is None:
            prop_range = [data[prop].min(), data[prop].max()]
        
        vmin = prop_range[0]
        vmax = prop_range[1]
        prop_values = data[prop][mask].values
        
        if log_prop:
            prop_values = np.log10(prop_values)
            vmin = np.log10(vmin)
            vmax = np.log10(vmax)
            
        plt.scatter(np.log10(data['S1_mass'][mask]),
                    np.log10(data['orbital_period'][mask]),
                    c=prop_values,
                    s=s,
                    vmin=vmin,
                    vmax=vmax,
                    marker='v',
                    cmap='viridis',
                    alpha=alpha,
                    zorder=0.5)
        if prop in DEFAULT_LABELS:
            if log_prop:
                label = DEFAULT_LABELS[prop][1]
            else:
                label = DEFAULT_LABELS[prop][0]
        else:
            label = prop
            if log_prop:
                label = 'log10_'+label
        plt.colorbar(label=label, orientation='horizontal')
    else:
        plt.scatter(np.log10(data['S1_mass'][mask]),
                    np.log10(data['orbital_period'][mask]),
                    s=s,
                    marker='v',
                    color='black',
                    alpha=alpha,
                    zorder=0.5)
    
def plot_grid_slice(grid, slice_3D_var_str, slice_3D_var_range, termination_flag='combined_TF12', PLOT_PROPERTIES=None):
    grid.plot2D('star_1_mass', 'period_days', None,
                termination_flag=termination_flag,
                grid_3D=True, slice_3D_var_str=slice_3D_var_str,
                slice_3D_var_range=slice_3D_var_range,
                verbose=False, **PLOT_PROPERTIES)