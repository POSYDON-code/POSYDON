__authors__ = [
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Max Briel <max.briel@unige.ch>",
    ]

import os
import numpy as np
import matplotlib as mpl
import pandas as pd
from posydon.utils.common_functions import PATH_TO_POSYDON
from posydon.visualization.plot_defaults import DEFAULT_LABELS
from posydon.utils.constants import Zsun
from posydon.grids.psygrid import PSyGrid
from posydon.utils.common_functions import convert_metallicity_to_string
from posydon.utils.constants import Zsun
import matplotlib.pyplot as plt
from posydon.visualization.plot_defaults import DEFAULT_LABELS

PATH_TO_POSYDON_DATA = os.environ.get("PATH_TO_POSYDON_DATA",'./')

plt.style.use(os.path.join(PATH_TO_POSYDON, "posydon/visualization/posydon.mplstyle"))

cm = mpl.colormaps.get_cmap('tab20')
COLORS = [cm.colors[i] for i in range(len(cm.colors)) if i%2==0] + [cm.colors[i] for i in range(len(cm.colors)) if i%2==1]

def plot_merger_rate_density(z, rate_density, zmax=10., channels=False,
                             GWTC3=False, label='DCO', **kwargs):

    plt.plot(z[z<zmax], rate_density['total'][z<zmax], label=f'{label} total', color='black')
    if channels:
        for i, ch in enumerate([key for key in rate_density.keys() if key != 'total']):
            if ch != 'total':
                plt.plot(z[z<zmax], rate_density[ch][z<zmax], label=ch, color=COLORS[i])

    if GWTC3:
        plt.errorbar(0.2,17.9, yerr=[[0],[26.1]], fmt='',color='black', label='$\mathcal{R}_\mathrm{BBH}$ GWTC-3')

def plot_grb_rate_density(z, rate_density, zmax=10., channels=False,
                          grb_components=False, Perley16=False, **kwargs):

    plt.plot(z[z<zmax], rate_density['total'][z<zmax], label='GRB total', color='black', linestyle='--')
    if grb_components:
        plt.plot(z[z<zmax], rate_density['total_GRB1'][z<zmax], label='total_GRB1', color='black', linestyle=':')
        plt.plot(z[z<zmax], rate_density['total_GRB2'][z<zmax], label='total_GRB2', color='black', linestyle='-.')
    if channels:
        for i, ch in enumerate([key for key in rate_density.keys() if (key != 'total' and 'GRB' not in key)]):
            if 'total' not in ch and 'GRB' not in ch:
                if sum(rate_density[ch][z<zmax]) == 0.:
                    continue
                if 'label_channels' in kwargs and kwargs['label_channels']:
                    label = ch
                else:
                    label = None
                plt.plot(z[z<zmax], rate_density[ch][z<zmax], label=label, color=COLORS[i], linestyle='--')
                if grb_components:
                    if 'label_channels' in kwargs and kwargs['label_channels']:
                        label1 = ch+'_GRB1'
                        label2 = ch+'_GRB2'
                    else:
                        label1 = None
                        label2 = None
                    plt.plot(z[z<zmax], rate_density[ch+'_GRB1'][z<zmax], label=label1, color=COLORS[i], linestyle=':')
                    plt.plot(z[z<zmax], rate_density[ch+'_GRB2'][z<zmax], label=label2, color=COLORS[i], linestyle='-.')

    # Perley et al. (2016)
    if Perley16:
        z_P16 = np.array([0.3,0.75,1.35,1.75,2.3,3.1,4,5.2,7.])
        rate_P16 = np.array([2.07e-10,1.46e-9,1.816e-9,3.93e-9,8.27e-9,6.90e-9,4.20e-9,3.74e-9,1.11e-9])*1e9
        rate_P16_lower = rate_P16-np.array([1.12e-10,1.08e-9,1.30e-9,3.02e-9,6.78e-9,5.29e-9,2.82e-9,2.397e-9,5.96e-10])*1e9
        rate_P16_upper = np.array([5.40e-10,2.03e-9,2.60e-9,5.29e-9,1.03e-8,9.13e-9,6.78e-9,6.78e-9,4.27e-9])*1e9-rate_P16
        rate_P16_error_y =np.array([rate_P16_lower.tolist(),rate_P16_upper.tolist()])
        rate_P16_left = z_P16-np.array([0.1,0.5,1,1.5,2.,2.6,3.5,4.5,6])
        rate_P16_right = np.array([0.5,1,1.5,2.,2.6,3.5,4.5,6,8])-z_P16
        rate_P16_error_x =np.array([rate_P16_left.tolist(),rate_P16_right.tolist()])
        plt.errorbar(z_P16,rate_P16,xerr=rate_P16_error_x,yerr=rate_P16_error_y, fmt='.', color='black',
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
    title = r'Merger efficiency'
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




def plot_popsyn_over_grid_slice(pop, grid_type, met_Zsun, slices=None, channel=None,
                                plot_dir='./', prop=None, prop_range=None,
                                log_prop=False, alpha=0.3, s=5.,
                                show_fig=True, save_fig=True, close_fig=True,
                                plot_extension='png', verbose=False):

    # Check if step_names in pop.history data
    if 'step_names' not in pop.history.columns:
        raise ValueError('Formation channel information not available in popsynth data.')

    # load grid
    met = convert_metallicity_to_string(met_Zsun)
    grid_path = os.path.join(PATH_TO_POSYDON_DATA, f'POSYDON_data/{grid_type}/{met}_Zsun.h5')
    grid = PSyGrid(verbose=False)
    grid.load(grid_path)

    # check if formation channel information is avaialbe
    if channel is not None and 'channel' not in pop.columns:
        raise ValueError('Formation channel information not available in popsynth data.')

    if 'CO' in grid_path:
        # compact object mass slices
        # TODO: THIS SELECTION DOES NOT WORK!
        m_COs = np.unique(np.around(grid.initial_values['star_2_mass'],1))
        m_COs_edges = 10**((np.log10(np.array(m_COs)[1:])+np.log10(np.array(m_COs)[:-1]))*0.5)
        m2 = [0.]+m_COs_edges.tolist()+[2*m_COs_edges[-1]]
        vars = m_COs
        fname = 'grid_m_%1.2f.' + plot_extension
        title = '$m_\mathrm{CO}=%1.2f\,M_\odot$'
        slice_3D_var_str='star_2_mass'
    elif 'HMS-HMS' in grid_path:
        # mass ratio slices
        qs = np.unique(np.around(grid.initial_values['star_2_mass']/grid.initial_values['star_1_mass'],2))
        dq = (qs[1:]-qs[:-1])*0.5
        dq_edges = (qs[:-1]+dq).tolist()
        dq_edges = [0.]+dq_edges+[1.]
        vars = qs.tolist()
        fname = 'grid_q_%1.2f.' + plot_extension
        title = '$q=%1.2f$'
        slice_3D_var_str='mass_ratio'
    else:
        raise ValueError('Grid type not supported!')


    if channel is not None:
        channel_sel = 'channel == '+str(channel)
    else:
        channel_sel = ''

    for i, var in enumerate(vars):
        if slices is not None and round(var,2) not in np.around(slices,2):
            if verbose:
                print(f'Skipping {round(var,2)}')
            continue

        met_indices = np.where(pop.select(where=channel_sel, columns=['metallicity']) == met_Zsun)[0].tolist()
        if 'HMS-HMS' in grid_path:
            slice_3D_var_range = (dq_edges[i],dq_edges[i+1])
            if len(met_indices) == 0:
                data = pd.DataFrame(columns=['S1_mass','S2_mass', 'orbital_period'])
            else:
                sel = 'index in '+str(met_indices)+' & event == "ZAMS"'
                data = pop.history.select(where=sel, columns=['S1_mass','S2_mass', 'orbital_period'])

            q = data['S2_mass'].values/data['S1_mass'].values
            q[q>1] = 1./q[q>1]
            mask = (q>=slice_3D_var_range[0]) & (q<=slice_3D_var_range[1])
        elif 'CO' in grid_path:
            if i == len(m2):
                continue
            slice_3D_var_range = (m2[i],m2[i+1])
            # select popsynth binaries in the given compact object mass
            # TODO: implement the case of reversal mass ratio
            if 'CO-HMS_RLO' in grid_path:
                if len(met_indices) == 0:
                    data = pd.DataFrame(columns=['S1_mass','S2_mass', 'orbital_period'])
                else:
                    sel = 'index in '+str(met_indices)+' & event == "oRLO2"'
                    data = pop.history.select(where=sel, columns=['S1_mass','S2_mass', 'orbital_period'])
                
                m_CO = data['S1_mass'].values
                mask = (m_CO>=slice_3D_var_range[0]) & (m_CO<=slice_3D_var_range[1])

            elif 'CO-HeMS' in grid_path:
                if len(met_indices) == 0:
                    data = pd.DataFrame(columns=['S1_mass','S2_mass', 'orbital_period'])
                else:
                    sel = 'index in '+str(met_indices)+' & step_names == "step_CE"'
                    data = pop.history.select(where=sel, columns=['S1_mass','S2_mass', 'orbital_period'])
                    
                m_CO = data['S1_mass'].values
                mask = (m_CO>=slice_3D_var_range[0]) & (m_CO<=slice_3D_var_range[1])
            else:
                raise ValueError('Grid type not supported!')
        # TODO: skip plotting slice if there are no data
        try:
            PLOT_PROPERTIES = {
                'figsize' : (4,3.5) if prop is None else (4,5.5),
                'path_to_file' : plot_dir,
                'show_fig' : False,
                'close_fig' : False,
                #'fname' : fname%var,
                'title' : title%var if channel is None else title%var + '\n' + channel,
                'log10_x' : True,
                'log10_y' : True,
                'legend2D' : {'bbox_to_anchor' : (1.03, 0.5)},
            }
            # grid
            grid.plot2D('star_1_mass', 'period_days', None,
                            termination_flag='combined_TF12',
                            grid_3D=True, slice_3D_var_str=slice_3D_var_str,
                            slice_3D_var_range=slice_3D_var_range,
                            verbose=False, **PLOT_PROPERTIES)

            log10_m1 = np.log10(data.loc[mask,'S1_mass'].values)
            log10_p = np.log10(data.loc[mask,'orbital_period'].values)
            # plot color map of a given DCO variable
            if prop is not None:
                if prop not in pop.columns:
                    raise ValueError(f'Property {prop} not available in popsynth data.')

                # basically only for DCO systems
                met_indices = np.array(met_indices)[mask].tolist()
                prop_sel = 'index in '+str(met_indices) + ' & event == "CO_contact"'
                prop_data = pop.history.select(where=prop_sel, columns=[prop])
                prop_values = prop_data[prop].values
                vmin = prop_range[0]
                vmax = prop_range[1]
                if log_prop:
                    prop_values = np.log10(prop_values)
                    vmin = np.log10(vmin)
                    vmax = np.log10(vmax)
                plt.scatter(log10_m1, log10_p, c=prop_values, s=s, vmin=vmin, vmax=vmax, marker='v', cmap='viridis', alpha=alpha, zorder=0.5)
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
                plt.scatter(log10_m1, log10_p, s=s, marker='v', color='black', alpha=alpha, zorder=0.5)

            if save_fig:
                plt.savefig(os.path.join(plot_dir, fname%var), bbox_inches='tight')
            if show_fig:
                plt.show()
            if close_fig:
                plt.close()
        except:
            raise ValueError(f'Failed to plot slice %1.2f'%var)
