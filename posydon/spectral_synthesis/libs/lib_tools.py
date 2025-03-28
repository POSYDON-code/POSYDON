"""
Tools and functions to be used for the Spec Lib functionalities

"""

__authors__ = [
    "Eirini Kasdagli <kasdaglie@ufl.edu>"
]
import itertools
import os
import h5py
import numpy as np
import pymsg
from scipy.spatial import KDTree


MSG_DIR = os.environ['MSG_DIR']
GRID_DIR = os.path.join(MSG_DIR, 'data', 'grids')
grid_cache = {}


def spec_grid(name):
    """Place docstring here."""
    specgrid_filepath = os.path.join(GRID_DIR, name)
    return pymsg.SpecGrid(specgrid_filepath)

def get_grid(name):
    grid = spec_grid(name)
    grid_path = os.path.join(GRID_DIR, name)
    grid_val = {}
    axis_to_ignore = ['v_lin_seq']
    for label in grid.axis_labels:
        min_ax = grid.axis_x_min[label]
        max_ax = grid.axis_x_max[label]
        with h5py.File(grid_path,"r") as f:
            keys = f['vgrid'].keys()
            for key in keys:
                if key not in axis_to_ignore:
                    data = list(f['vgrid'][key]['x'])
                    min_data = min(data)
                    max_data = max(data)
                    if (min_data == min_ax ) & (max_data == max_ax):
                        grid_val[label] = data
                        axis_to_ignore.append(key)
    
    return grid_val


def get_grid_points(name):
    grid = spec_grid(name)
    grid_val = get_grid(name)
    combinations = itertools.product(*grid_val.values())
    labels = grid_val.keys()
    combination_dicts = [dict(zip(labels, values)) for values in combinations]
    grid_points = []
    lam = np.linspace(3000, 7000, 1000)
    lam_c = 0.5*(lam[1:] + lam[:-1])
    for x in combination_dicts:
        try:
            F_test = np.asarray(grid.flux(x, lam))
            grid_points.append(x)
        except Exception as e:
            pass
            
    return grid_points


def get_nearest_neighbor(name,x):
    if name not in grid_cache:
        grid_cache[name] = get_grid_points(name)  # Cache grid points
        
    grid_points = grid_cache[name]
    labels = list(grid_points[0].keys())
    points = np.array([[point[l] for l in labels] for point in grid_points])
    #The values we need to find the NN 
    nn_point = np.array([x[l] for l in labels])
    tree = KDTree(points)
    _, index = tree.query(nn_point)
    nn_x = grid_points[index]
    return nn_x