"""
Simple utility for generating detection weights

Uses grid of detection probabilities to estimate detection probabilities

Anticipates data as Pandas dataframe with series ['m1', 'q', 'z', 'chieff']
"""

__authors__ = ["Michael Zevin <michael.zevin@ligo.org>",
               "Simone Bavera <Simone.Bavera@unige.ch>",]

import time

import numpy as np
import pandas as pd
from sklearn.neighbors import KNeighborsRegressor


class KNNmodel():
    """
    K-nearest neighbor model that instantiates based on detection probability grid

    When instantiating, must supply path to the grid, and key that represents
    GW network and sensitivity.
    """

    def __init__(self, grid_path, sensitivity_key, verbose=False):
        """
        Instantiates KNNmodel class and trains the KNN.

        grid_path : string
            Path to grid of detection probabilities.

        sensitivity_key : string
            GW detector sensitivity and network configuration you want to use,
                see arXiv:1304.0670v3
            detector sensitivities are taken from: https://dcc.ligo.org/LIGO-T2000012-v2/public
                available sensitivity keys (for Hanford, Livingston, Virgo network):
                    'O3actual_H1L1V1' : aligo_O3actual_H1.txt, aligo_O3actual_L1.txt, avirgo_O3actual.txt
                    'O4low_H1L1V1' : aligo_O4low.txt, aligo_O4low.txt, avirgo_O4high_NEW.txt
                    'O4high_H1L1V1' : aligo_O4high.txt, aligo_O4high.txt, avirgo_O4high_NEW.txt
                    'design_H1L1V1' : AplusDesign.txt, AplusDesign.txt, avirgo_O5high_NEW.txt
                detection probabilities are calculated using the IMRPhenomXHM approximant with a network SNR threshold of 10

        verbose : boolean
            Adds verbosity.
        """
        start = time.time()
        if verbose:
            print('\ntraining nearest neighbor algorithm...')

        # read in grid
        grid = pd.read_hdf(grid_path, key=sensitivity_key)

        # get values from grid for training
        pdets = np.asarray(grid['pdet'])
        m1_grid = np.asarray(grid['m1'])
        q_grid = np.asarray(grid['q'])
        z_grid = np.asarray(grid['z'])
        chieff_grid = np.asarray(grid['chieff'])

        # get bounds based on grid
        m1_bounds = (np.round(m1_grid.min(), 5), np.round(m1_grid.max(), 5))
        q_bounds = (np.round(q_grid.min(), 5), np.round(q_grid.max(), 5))
        z_bounds = (np.round(z_grid.min(), 5), np.round(z_grid.max(), 5))
        chieff_bounds = (np.round(chieff_grid.min(), 5), np.round(chieff_grid.max(), 5))
        self.m1_bounds, self.q_bounds, self.z_bounds, self.chieff_bounds = m1_bounds, \
            q_bounds, z_bounds, chieff_bounds

        # normalize to unit cube
        logm1_grid_norm = self.normalize(np.log10(m1_grid), np.log10(m1_bounds[0]), np.log10(m1_bounds[1]))
        q_grid_norm = self.normalize(q_grid, q_bounds[0], q_bounds[1])
        logz_grid_norm = self.normalize(np.log10(z_grid), np.log10(z_bounds[0]), np.log10(z_bounds[1]))
        chieff_grid_norm = self.normalize(chieff_grid, chieff_bounds[0], chieff_bounds[1])

        # train nearest neighbor algorithm
        X = np.transpose(np.vstack([logm1_grid_norm, q_grid_norm, logz_grid_norm, chieff_grid_norm]))
        y = np.transpose(np.atleast_2d(pdets))
        nbrs = KNeighborsRegressor(n_neighbors=10, weights='distance', algorithm='ball_tree', leaf_size=30, p=2, metric='minkowski')
        nbrs.fit(X, y)

        self.model = nbrs
        if verbose:
            print('   finished! It took {:0.2f}s to train the model on {:d} systems\n'.format(time.time()-start, len(pdets)))


    def predict_pdet(self, data, verbose=False):
        """
        Gives relative weight to each system in `data` based on its proximity to the points on the grid.
        Each system in `data` should have a primary mass `m1`, mass ratio `q`, redshift `z`, and effective spin `chieff`
        This function will determine detection probabilities using nearest neighbor algorithm in [log(m1), q, log(z), chieff] space
        Need to specify bounds (based on the trained grid) so that the grid and data get normalized properly

        data : Pandas dataframe
            Data you wish to predict detection probabilities for.
            Required series in the dataframe:
                'm1' : primary source-frame mass
                'q' : mass ratio (secondary mass/primary mass)
                'z' : redshift of merger
                'chieff' : effective inspiral spin

        verbose : boolean
            Adds verbosity.
        """
        start = time.time()
        if verbose:
            print('determining detection probabilities for data...')

        # get values from dataset and normalize
        m1_data = np.asarray(data['m1'])
        q_data = np.asarray(data['q'])
        z_data = np.asarray(data['z'])
        chieff_data = np.asarray(data['chieff'])

        logm1_data_norm = self.normalize(np.log10(m1_data), np.log10(self.m1_bounds[0]), np.log10(self.m1_bounds[1]))
        q_data_norm = self.normalize(q_data, self.q_bounds[0], self.q_bounds[1])
        logz_data_norm = self.normalize(np.log10(z_data), np.log10(self.z_bounds[0]), np.log10(self.z_bounds[1]))
        chieff_data_norm = self.normalize(chieff_data, self.chieff_bounds[0], self.chieff_bounds[1])

        # get pdets for the testing data
        X_fit = np.transpose(np.vstack([logm1_data_norm, q_data_norm,
                                        logz_data_norm, chieff_data_norm]))
        pdets = self.model.predict(X_fit).flatten()
        assert all([((p<=1) & (p>=0)) for p in pdets]), 'pdet is not between 0 and 1'

        if verbose:
            print('   finished! It took {:0.2f}s to determine detection probabilities for {:d} systems'.format(time.time()-start, len(X_fit)))
        return pdets


    @staticmethod
    def normalize(x, xmin, xmax, a=0, b=1):
        """
        normalizes data on range [a,b]
        """
        data_norm = (b-a)*(x-xmin) / (xmax-xmin) + a
        return data_norm
