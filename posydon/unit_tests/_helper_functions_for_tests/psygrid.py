"""Helper function(s) for tests requiring a PSyGrid

used in:
    - posydon/unit_tests/girds/test_psygrid.py
    - posydon/unit_tests/girds/test_post_processing.py
"""

__authors__ = [
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>"
]

import os
import numpy as np
import h5py
import json
from posydon.grids.psygrid import PSyGrid, H5_UNICODE_DTYPE

# helper functions
def get_PSyGrid(dir_path, idx, binary_history, star_history, profile,\
                n_runs=6, initial_values_dtypes=None,\
                final_values_dtypes=None):
    """Create a PSyGrid file with some runs

    Parameters
    ----------
    dir_path : str
        The path to the directory where to create the PSyGrid.
    idx : int
        Grid index.
    binary_history : ndarray
        Binary history data.
    star_history : ndarray
        Star history data.
    profile : ndarray
        Final profile data.
    n_runs : int (default: 6)
        Number of runs with data. run0 has no data, thus n_runs+1 entries in
        the grid.
    initial_values_dtypes : dtype specifcation or None (default: None)
        If specified infers the initial dtypes from this list.
    final_values_dtypes : dtype specifcation or None (default: None)
        If specified infers the final dtypes from this list.

    Returns
    -------
    str
        Path to the PSyGrid file.

    """
    path = os.path.join(dir_path, f"grid{idx}.h5")
    Grid = PSyGrid()
    Grid.filepath = path
    Grid.generate_config()
    MESA_paths = []
    H5_STRING = h5py.string_dtype()
    with h5py.File(path, "w") as hdf5_file:
        for r in range(n_runs+1):
            if r == 0:
                MESA_paths.append(os.path.join(dir_path, "m1_1.0_m2_1.0"\
                 +"_initial_period_in_days_0.0_initial_z_0.01_idx_0"))
                hdf5_file.create_group("/grid/run0/")
                continue
            MESA_paths.append(os.path.join(dir_path, "m1_1.0_m2_1.0"\
             +"_initial_period_in_days_{}_initial_z_0.01_idx_{}".format(\
             binary_history['period_days'][0]+r-1, r)))
            hdf5_file.create_dataset(f"/grid/run{r}/binary_history",\
                                     data=binary_history)
            if r == 4:
                new_star_history = np.copy(star_history)
                new_star_history['center_he4'][-1] = 0.1
                hdf5_file.create_dataset(f"/grid/run{r}/history1",\
                                         data=new_star_history)
                hdf5_file.create_dataset(f"/grid/run{r}/history2",\
                                         data=new_star_history)
            else:
                hdf5_file.create_dataset(f"/grid/run{r}/history1",\
                                         data=star_history)
                hdf5_file.create_dataset(f"/grid/run{r}/history2",\
                                         data=star_history)
            hdf5_file.create_dataset(f"/grid/run{r}/final_profile1",\
                                     data=profile)
            hdf5_file.create_dataset(f"/grid/run{r}/final_profile2",\
                                     data=profile)
        BH_types = binary_history.dtype.names
        SH_types = star_history.dtype.names
        if initial_values_dtypes is None:
            ini_dtypes = [('Z', '<f8'), ('period_days', '<f8')]
        else:
            ini_dtypes = initial_values_dtypes
        ini_vals = [tuple([np.nan]*len(ini_dtypes))]
        for i in range(n_runs):
            ini_vals +=[tuple([binary_history[k][0]+i\
                               if k in BH_types else 1.0\
                               for (k, t) in ini_dtypes])]
        ini_val = np.array(ini_vals, dtype=ini_dtypes)
        hdf5_file.create_dataset("/grid/initial_values", data=ini_val)
        # either required by post_process_grid or BinaryStar.from_run
        if final_values_dtypes is None:
            fin_dtypes = [('star_1_mass', '<f8'), ('star_2_mass', '<f8'),\
                          ('period_days', '<f8'), ('age', '<f8'),\
                          ('lg_mstar_dot_1', '<f8'),\
                          ('lg_mstar_dot_2', '<f8'),\
                          ('lg_system_mdot_1', '<f8'),\
                          ('lg_system_mdot_2', '<f8'),\
                          ('lg_wind_mdot_1', '<f8'),\
                          ('lg_wind_mdot_2', '<f8'), ('S1_log_R', '<f8'),\
                          ('S2_log_R', '<f8'), ('S1_center_h1', '<f8'),\
                          ('S2_center_h1', '<f8'), ('S1_center_he4', '<f8'),\
                          ('S2_center_he4', '<f8'), ('S1_center_c12', '<f8'),\
                          ('S2_center_c12', '<f8'), ('S1_center_n14', '<f8'),\
                          ('S2_center_n14', '<f8'), ('S1_center_o16', '<f8'),\
                          ('S2_center_o16', '<f8'), ('S1_surface_h1', '<f8'),\
                          ('S2_surface_h1', '<f8'), ('S1_surface_he4', '<f8'),\
                          ('S2_surface_he4', '<f8'),\
                          ('S1_surface_c12', '<f8'),\
                          ('S2_surface_c12', '<f8'),\
                          ('S1_surface_n14', '<f8'),\
                          ('S2_surface_n14', '<f8'),\
                          ('S1_surface_o16', '<f8'),\
                          ('S2_surface_o16', '<f8'), ('S1_log_LH', '<f8'),\
                          ('S2_log_LH', '<f8'), ('S1_log_LHe', '<f8'),\
                          ('S2_log_LHe', '<f8'), ('S1_log_Lnuc', '<f8'),\
                          ('S2_log_Lnuc', '<f8'), ('S1_center_gamma', '<f8'),\
                          ('S2_center_gamma', '<f8'),\
                          ('S1_he_core_mass', '<f8'),\
                          ('S2_he_core_mass', '<f8'),\
                          ('S1_co_core_mass', '<f8'),\
                          ('S2_co_core_mass', '<f8'),\
                          ('S1_avg_c_in_c_core_at_He_depletion', '<f8'),\
                          ('S1_co_core_mass_at_He_depletion', '<f8'),\
                          ('termination_flag_1', H5_STRING),\
                          ('termination_flag_2', H5_STRING),\
                          ('interpolation_class', H5_STRING)]
        else:
            fin_dtypes = final_values_dtypes
        # run0: initial MT
        fin_vals_run = []
        for (k, t) in fin_dtypes:
            if k in ['S1_avg_c_in_c_core_at_He_depletion',\
                     'S1_co_core_mass_at_He_depletion']:
                fin_vals_run += [None]
            elif k == 'termination_flag_1':
                fin_vals_run += ["Terminate because of overflowing initial "\
                                 +"model"]
            elif k == 'termination_flag_2':
                fin_vals_run += ["initial_RLOF"]
            elif k == 'interpolation_class':
                fin_vals_run += ["initial_MT"]
            else:
                fin_vals_run += [np.nan]
        fin_vals = [tuple(fin_vals_run)]
        if n_runs >= 1:
            # run1: no MT, WD
            fin_vals_run = []
            for (k, t) in fin_dtypes:
                if k in BH_types:
                    fin_vals_run += [binary_history[k][-1]]
                elif k[3:] in SH_types:
                    fin_vals_run += [star_history[k[3:]][-1]]
                elif k in ['S1_avg_c_in_c_core_at_He_depletion',\
                           'S1_co_core_mass_at_He_depletion']:
                    fin_vals_run += [None]
                elif k == 'termination_flag_1':
                    fin_vals_run += ["gamma_center_limit"]
                elif k == 'termination_flag_2':
                    fin_vals_run += ["no_RLOF"]
                elif k == 'interpolation_class':
                    fin_vals_run += ["no_MT"]
                else:
                    fin_vals_run += [np.nan]
            fin_vals += [tuple(fin_vals_run)]
        if n_runs >= 2:
            # run2: stable MT, WD error
            fin_vals_run = []
            for (k, t) in fin_dtypes:
                if k in ['S1_avg_c_in_c_core_at_He_depletion',\
                         'S1_co_core_mass_at_He_depletion']:
                    fin_vals_run += [None]
                elif k == 'termination_flag_1':
                    fin_vals_run += ["gamma_center_limit"]
                elif k == 'termination_flag_2':
                    fin_vals_run += ["case_A1"]
                elif k == 'interpolation_class':
                    fin_vals_run += ["stable_MT"]
                else:
                    fin_vals_run += [np.nan]
            fin_vals += [tuple(fin_vals_run)]
        if n_runs >= 3:
            # run3: stable reverse MT, star 2 becomes WD
            fin_vals_run = []
            for (k, t) in fin_dtypes:
                if k == 'S1_center_gamma':
                    fin_vals_run += [0.1]
                elif k in BH_types:
                    fin_vals_run += [binary_history[k][-1]]
                elif k[3:] in SH_types:
                    fin_vals_run += [star_history[k[3:]][-1]]
                elif k in ['S1_avg_c_in_c_core_at_He_depletion',\
                           'S1_co_core_mass_at_He_depletion']:
                    fin_vals_run += [None]
                elif k == 'termination_flag_1':
                    fin_vals_run += ["gamma_center_limit"]
                elif k == 'termination_flag_2':
                    fin_vals_run += ["case_A1/A2"]
                elif k == 'interpolation_class':
                    fin_vals_run += ["stable_reverse_MT"]
                else:
                    fin_vals_run += [np.nan]
            fin_vals += [tuple(fin_vals_run)]
        if n_runs >= 4:
            # run4: stable MT, max age
            fin_vals_run = []
            for (k, t) in fin_dtypes:
                if k in BH_types:
                    fin_vals_run += [binary_history[k][-1]]
                elif k[3:] in SH_types:
                    fin_vals_run += [star_history[k[3:]][-1]]
                elif k in ['S1_avg_c_in_c_core_at_He_depletion',\
                           'S1_co_core_mass_at_He_depletion']:
                    fin_vals_run += [None]
                elif k == 'termination_flag_1':
                    fin_vals_run += ["max_age"]
                elif k == 'termination_flag_2':
                    fin_vals_run += ["case_B1"]
                elif k == 'interpolation_class':
                    fin_vals_run += ["stable_MT"]
                else:
                    fin_vals_run += [np.nan]
            fin_vals += [tuple(fin_vals_run)]
        if n_runs >= 5:
            # run5: stable reverse MT, CC of star 2 with values at He depletion
            fin_vals_run = []
            for (k, t) in fin_dtypes:
                if k in BH_types:
                    fin_vals_run += [binary_history[k][-1]]
                elif k[3:] in SH_types:
                    fin_vals_run += [star_history[k[3:]][-1]]
                elif k == 'S1_avg_c_in_c_core_at_He_depletion':
                    fin_vals_run += [0.4]
                elif k == 'S1_co_core_mass_at_He_depletion':
                    fin_vals_run += [1.0]
                elif k == 'termination_flag_1':
                    fin_vals_run += ["Secondary enters pulsational "\
                                     +"pair-instability regime"]
                elif k == 'termination_flag_2':
                    fin_vals_run += ["case_B1/B2"]
                elif k == 'interpolation_class':
                    fin_vals_run += ["stable_reverse_MT"]
                else:
                    fin_vals_run += [np.nan]
            fin_vals += [tuple(fin_vals_run)]
        if n_runs >= 6:
            # run6: unstable MT
            fin_vals_run = []
            for (k, t) in fin_dtypes:
                if k in BH_types:
                    fin_vals_run += [binary_history[k][-1]]
                elif k[3:] in SH_types:
                    fin_vals_run += [star_history[k[3:]][-1]]
                elif k in ['S1_avg_c_in_c_core_at_He_depletion',\
                           'S1_co_core_mass_at_He_depletion']:
                    fin_vals_run += [None]
                elif k == 'termination_flag_1':
                    fin_vals_run += ["Reached maximum mass transfer rate: "\
                                     +"1d-1"]
                elif k == 'termination_flag_2':
                    fin_vals_run += ["case_C1"]
                elif k == 'interpolation_class':
                    fin_vals_run += ["unstable_MT"]
                else:
                    fin_vals_run += [np.nan]
            fin_vals += [tuple(fin_vals_run)]
        for i in range(7, n_runs+1):
            # ith run: no MT, WD
            fin_vals_run = []
            for (k, t) in fin_dtypes:
                if k in BH_types:
                    fin_vals_run += [binary_history[k][-1]]
                elif k[3:] in SH_types:
                    fin_vals_run += [star_history[k[3:]][-1]]
                elif k in ['S1_avg_c_in_c_core_at_He_depletion',\
                           'S1_co_core_mass_at_He_depletion']:
                    fin_vals_run += [None]
                elif k == 'termination_flag_1':
                    fin_vals_run += ["gamma_center_limit"]
                elif k == 'termination_flag_2':
                    fin_vals_run += ["no_RLOF"]
                elif k == 'interpolation_class':
                    fin_vals_run += ["no_MT"]
                else:
                    fin_vals_run += [np.nan]
            fin_vals += [tuple(fin_vals_run)]
        fin_val = np.array(fin_vals, dtype=fin_dtypes)
        hdf5_file.create_dataset("/grid/final_values", data=fin_val)
        hdf5_file.attrs["config"] = json.dumps(str(dict(Grid.config)))
        rel_paths = np.array(MESA_paths, dtype=H5_STRING)
        hdf5_file.create_dataset("relative_file_paths", data=rel_paths)
    return path


def get_simple_PSyGrid(dir_path, idx, binary_history, star_history, profile):
    """Create a PSyGrid file with two runs (+ the empty run0) and reduced
    columns in the initial and final values.

    Parameters
    ----------
    dir_path : str
        The path to the directory where to create the PSyGrid.
    idx : int
        Grid index.
    binary_history : ndarray
        Binary history data.
    star_history : ndarray
        Star history data.
    profile : ndarray
        Final profile data.

    Returns
    -------
    str
        Path to the PSyGrid file.

    """
    return get_PSyGrid(dir_path, idx, binary_history, star_history, profile,\
                       n_runs=2, initial_values_dtypes=[('period_days',\
                                                         '<f8')],\
                       final_values_dtypes=[('period_days', '<f8'),\
                                            ('termination_flag_1',\
                                             h5py.string_dtype())])
