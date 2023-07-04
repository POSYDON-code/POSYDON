__authors__ = ["Simone Bavera <Simone.Bavera@unige.ch>"]
    
import os
import numpy as np
import matplotlib.pyplot as plt
from posydon.utils.common_functions import PATH_TO_POSYDON
from posydon.visualization.plot_defaults import DEFAULT_LABELS
from posydon.utils.constants import Zsun

plt.style.use(os.path.join(PATH_TO_POSYDON, "posydon/visualization/posydon.mplstyle"))

COLORS = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
          'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']
COLORS *= 10

def plot_merger_rate_density(z, rate_density, ylim=(1e-1,4e2), zmax=10., show=True, 
                             path=None, channels=False, GWTC3=False, **kwargs):

    plt.plot(z[z<zmax], rate_density['total'][z<zmax], label='BBH total', color='black')
    if channels:
        for ch in [key for key in rate_density.keys() if key != 'total']:
            if ch != 'total':
                plt.plot(z[z<zmax], rate_density[ch][z<zmax], label=ch)
    
    if GWTC3:
        plt.errorbar(0.2,17.9, yerr=[[0],[26.1]], fmt='',color='black', label='$\mathcal{R}_\mathrm{BBH}$ GWTC-3')
        
def plot_grb_rate_density(z, rate_density, ylim=(1e-1,4e2), zmax=10., show=True, path=None, channels=False,
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
        
def plot_rate_density(z_dco=None, R_dco=None, z_grb=None, R_grb=None, **kwargs):
    plt.figure()   
    title1 = ''
    title2 = '' 
    if z_dco is not None and R_dco is not None:
        title1 += 'DCO merger'
        plot_merger_rate_density(z_dco, R_dco, **kwargs)
        # do not display twice the channel labels
        kwargs['label_channels'] = False 
    if z_grb is not None and R_grb is not None:
        title2 += 'GRB beamed'
        plot_grb_rate_density(z_grb, R_grb, **kwargs)
    if title1 and title2:
        title = title1 + ' \& ' + title2 + ' rate densities'
    else:
        title = title1 + title2 + ' rate density'
    plt.title(title)
    plt.yscale('log')
    plt.ylabel(r'$\mathcal{R} \,[\mathrm{Gpc}^{-3}\,\mathrm{yr}^{-1}]$')
    plt.xlabel(r'$z$')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    if 'ylim' in kwargs:
        plt.ylim(kwargs['ylim'])
    if 'path' in kwargs and isinstance(kwargs['path'],str):
        plt.savefig(path)
    if 'show' in kwargs and kwargs['show']:
        plt.show()

def plot_merger_efficiency(met, merger_efficiency, show=True, path=None, channels=False):
    title = r'Merger efficiency'
    plt.figure()
    plt.title(title)
    plt.plot(met/Zsun, merger_efficiency['total'], label='total', color='black')
    if channels:
        for ch in merger_efficiency.keys():
            if ch != 'total':
                plt.plot(met/Zsun, merger_efficiency[ch], label=ch)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel(r'$Z/Z_\odot$')
    plt.ylabel(r'\#DCOs [$M_\odot^{-1}$]')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    if path:
        plt.savefig(path)
    if show:
        plt.show()

        
def plot_hist_properties(x, df_intrinsic=None, df_observable=None,
                        channel=None,
                        show=True, path=None, alpha=0.5,
                        range=None, bins=20, xlog=False, ylog=False,
                        **kwargs):
    if df_intrinsic is not None and df_observable is not None:
        title = r'Intrinsic vs. observable (dashed) population'
    elif df_intrinsic is not None:
        title = r'Intrinsic population'
    elif df_observable is not None:
        title = r'Observable population'
    else:
        raise ValueError('You should provide either an intrinsic or a '
                         'detectable population.')

    plt.figure()
    if df_intrinsic is not None:
        # reweight to have GRB1/GRB2 ratio
        df_intrinsic['weight'] /= sum(df_intrinsic['weight'])
        if channel is not None:
            sel = (df_intrinsic['weight'] > 0) & (df_intrinsic['channel'] == channel)
            title += f'\n{channel}'
        else:
            sel = df_intrinsic['weight'] > 0
        if isinstance(x, str):
            sel_tmp = sel & ~df_intrinsic[x].isna()
            values, weights = df_intrinsic.loc[sel_tmp,[x,'weight']].values.T
            if xlog:
                mask = values > 0
                values = np.log10(values[mask])
                weights = weights[mask]
            plt.hist(values, weights=weights, color=COLORS[0],
                     alpha=alpha, range=range, bins=bins)
        elif isinstance(x, list):
            for i, x_i in enumerate(x):
                if "S1" in x_i:
                    label = 'S1'
                elif "S2" in x_i:
                    label = 'S2'
                else:
                    label = x_i
                sel_tmp = sel & ~df_intrinsic[x_i].isna()
                values, weights = df_intrinsic.loc[sel_tmp,[x_i,'weight']].values.T
                if xlog:
                    mask = values > 0
                    values = np.log10(values[mask])
                    weights = weights[mask]
                plt.hist(values, weights=weights,
                         color=COLORS[i], label=label,
                         alpha=alpha, range=range, bins=bins)
    if df_observable is not None:
        # reweight to have GRB1/GRB2 ratio
        df_observable['weight'] /= sum(df_observable['weight'])
        if channel is not None:
            sel = (df_observable['weight'] > 0) & (df_observable['channel'] == channel)
            if channel not in title:
                title += f'\n{channel}'
        else:
            sel = df_observable['weight'] > 0
        if isinstance(x, str):
            sel_tmp = sel & ~df_observable[x].isna()
            values, weights = df_observable.loc[sel_tmp,[x,'weight']].values.T
            if xlog:
                mask = values > 0
                values = np.log10(values[mask])
                weights = weights[mask]
            plt.hist(values, weights=weights,
                     color=COLORS[0],
                     histtype=u'step', linestyle='--',
                     alpha=alpha, range=range, bins=bins)
        elif isinstance(x, list):
            for i, x_i in enumerate(x):
                if "S1" in x_i:
                    label = 'S1'
                elif "S2" in x_i:
                    label = 'S2'
                else:
                    label = x_i
                sel_tmp = sel & ~df_observable[x_i].isna()
                values, weights = df_observable.loc[sel_tmp,[x_i,'weight']].values.T
                if xlog:
                    mask = values > 0
                    values = np.log10(values[mask])
                    weights = weights[mask]
                plt.hist(values, weights=weights,
                         color=COLORS[i], label=label,
                         histtype=u'step', linestyle='--',
                         alpha=alpha, range=range, bins=bins)
    plt.title(title)
    if ylog:
        plt.yscale('log')
    plt.ylabel(r'PDF')
    try:
        if isinstance(x, str):
            if not xlog:
                plt.xlabel(DEFAULT_LABELS[x][0])
            else:
                plt.xlabel(DEFAULT_LABELS[x][1])
        else:
            if not xlog:
                plt.xlabel(DEFAULT_LABELS[x[0]][0])
            else:
                plt.xlabel(DEFAULT_LABELS[x[0]][1])
    except:
        if isinstance(x, str):
            if not xlog:
                plt.xlabel(x)
            else:
                plt.xlabel('log10_'+x)
        else:
            if not xlog:
                plt.xlabel(x[0])
            else:
                plt.xlabel('log10_'+x[0])
    if isinstance(x, list):
        plt.legend(loc=1)
    if path:
        plt.savefig(path)
    if show:
        plt.show()