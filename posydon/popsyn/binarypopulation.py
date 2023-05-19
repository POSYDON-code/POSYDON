"""Create, evolve and save a binary star population.

Large populations are RAM limited when holding an arbitrary
number of BinaryStar instances. Therefore, by default the BinaryPopulation
will generate binaries, evolve, and try to save them to disk one at a time.

Create a BinaryPopulation instance from an inifile:
I. CREATING A POPULATION
------------------------
a) One-liner for creating a BinaryPopulation from an inifile:

> BinaryPopulation.from_ini('<PATH_TO_POSYDON>' \
            '/posydon/popsyn/population_params_default.ini')
"""


__authors__ = [
    "Kyle Akira Rocha <kylerocha2024@u.northwestern.edu>",
    "Jeffrey Andrews <jeffrey.andrews@northwestern.edu>",
    "Konstantinos Kovlakas <Konstantinos.Kovlakas@unige.ch>",
    "Devina Misra <devina.misra@unige.ch>",
    "Simone Bavera <Simone.Bavera@unige.ch>",
]


__credits__ = ["Nam Tran <tranhn03@gmail.com>"]


import pandas as pd
import numpy as np
import warnings
import traceback
import atexit
import os
from tqdm import tqdm
import psutil
import random

from posydon.binary_evol.binarystar import BinaryStar
from posydon.binary_evol.singlestar import SingleStar
from posydon.binary_evol.simulationproperties import SimulationProperties
from posydon.popsyn.star_formation_history import get_formation_times

from posydon.popsyn.independent_sample import generate_independent_samples
from posydon.utils.common_functions import (orbital_period_from_separation,
                                            orbital_separation_from_period)
from posydon.popsyn.defaults import default_kwargs
from posydon.popsyn.io import binarypop_kwargs_from_ini
from posydon.utils.constants import Zsun


# 'event' usually 10 but 'detached (Integration failure)' can occur
HISTORY_MIN_ITEMSIZE = {'state': 30, 'event': 10, 'step_names': 15,
                        'S1_state': 31, 'S2_state': 31,
                        'mass_transfer_case': 7,
                        'S1_SN_type': 5, 'S2_SN_type': 5}
ONELINE_MIN_ITEMSIZE = {'state_i': 30, 'state_f': 30,
                        'event_i': 10, 'event_f': 10,
                        'step_names_i': 15, 'step_names_f': 15,
                        'S1_state_i': 31, 'S1_state_f': 31,
                        'S2_state_i': 31, 'S2_state_f': 31,
                        'mass_transfer_case_i': 7, 'mass_transfer_case_f': 7,
                        'S1_SN_type': 5, 'S2_SN_type': 5}

# BinaryPopulation will enforce a constant metallicity accross all steps that
# load stellar or binary models by checked this list of steps.
STEP_NAMES_LOADING_GRIDS = [
    'step_HMS_HMS', 'step_CO_HeMS', 'step_CO_HMS_RLO', 'step_detached','step_isolated','step_disrupted'
]

class BinaryPopulation:
    """Handle a binary star population."""

    def __init__(self, **kwargs):
        """Initialize the binary population object.

        Parameters
        ----------
        number_of_binaries : int
            Size of the population
        population_properties : SimulationProperties
            Instance of simulationproperties holding steps.
        """
        # Set population kwargs - first set defaults, then add updates
        self.kwargs = default_kwargs.copy()
        for key, arg in kwargs.items():
            self.kwargs[key] = arg
        # Have a binary fraction change the number_of binaries.
        #Eirini's change
        self.binary_fraction = self.kwargs.get('binary_fraction')
        self.number_of_binaries = self.kwargs.get('number_of_binaries')

        self.population_properties = self.kwargs.get('population_properties',
                                                     SimulationProperties())
        atexit.register(lambda: BinaryPopulation.close(self))

        self.population_properties.max_simulation_time = self.kwargs.get(
            'max_simulation_time')  # years

        entropy = self.kwargs.get('entropy', None)
        seq = np.random.SeedSequence(entropy=entropy)

        self.comm = self.kwargs.pop('comm', None)
        if self.comm is not None:
            self.rank = self.comm.Get_rank()
            self.size = self.comm.Get_size()
            seed_seq = [i for i in seq.spawn(self.size)][self.rank]
        else:
            seed_seq = seq

        self.RNG = np.random.default_rng(seed=seed_seq)
        self.kwargs['RNG'] = self.RNG
        self.entropy = self.RNG.bit_generator._seed_seq.entropy
        self.kwargs['entropy'] = self.entropy

        self.manager = PopulationManager(**self.kwargs)

        # use manager methods
        self.to_df = self.manager.to_df
        self.to_oneline_df = self.manager.to_oneline_df
        self.find_failed = self.manager.find_failed

    @classmethod
    def from_ini(cls, path, verbose=False):
        """Create a BinaryPopulation instance from an inifile.

        Parameters
        ----------
        path : str
            Path to an inifile to load in.

        verbose : bool
            Print useful info.

        Returns
        -------
        BinaryPopulation
            A new instance of a BinaryPopulation.
        """
        kwargs = binarypop_kwargs_from_ini(path, verbose=verbose)
        return cls(**kwargs)

    def evolve(self, **kwargs):
        """Evolve a binary population.

        Parameters
        ----------
        indices : list, optional
            Custom binary indices to use. Default is range(number_of_binaries).
            If running with MPI, indices are split between processes if given.
        breakdown_to_df : bool, True
            Breakdown a binary after evolution, converting to dataframe and
            removing the binary instance from memory.
        tqdm : bool, False
            Show tqdm progress bar during evolution.

        Returns
        -------
        None
        """
        # combine kw defined at init and any passed here
        kw = {**self.kwargs, **kwargs}
        tqdm_bool = kw.get('tqdm', False)
        breakdown_to_df_bool = kw.get('breakdown_to_df', True)
        from_hdf_bool = kw.get('from_hdf', False)

        if self.comm is None:   # do regular evolution
            indices = kw.get('indices',
                             list(range(self.number_of_binaries)))
            params = {'indices':indices,
                      'tqdm':tqdm_bool,
                      'breakdown_to_df':breakdown_to_df_bool,
                      'from_hdf':from_hdf_bool}
            self.kwargs.update(params)

            self._safe_evolve(**self.kwargs)
        else:
            # do MPI evolution
            indices = kw.get('indices',
                            list(range(self.number_of_binaries)))
            indices_split = np.array_split(indices, self.size)
            batch_indices = indices_split[self.rank]
            mpi_tqdm_bool = True if (tqdm_bool and self.rank == 0) else False

            params = {'indices':batch_indices,
                      'tqdm':mpi_tqdm_bool,
                      'breakdown_to_df':breakdown_to_df_bool,
                      'from_hdf':from_hdf_bool}
            self.kwargs.update(params)

            self._safe_evolve(**self.kwargs)

    def _safe_evolve(self, **kwargs):
        """Evolve binaries in a population, catching warnings/exceptions."""
        if not self.population_properties.steps_loaded:
            # Enforce the same metallicity for all grid steps
            for step_name, tup in self.population_properties.kwargs.items():

                if step_name in STEP_NAMES_LOADING_GRIDS:
                    step_function, step_kwargs = tup # unpack params
                    step_kwargs['metallicity'] = self.kwargs.get('metallicity', 1)

                    # update the step kwargs, override metallicity
                    modified_tup = (step_function, step_kwargs)
                    self.population_properties.kwargs[step_name] = modified_tup

            self.population_properties.load_steps()

        indices = kwargs.get('indices', list(range(self.number_of_binaries)))

        indices_for_iter = (tqdm(indices) if kwargs.get('tqdm', False)
                            else indices)
        breakdown_to_df = kwargs.get('breakdown_to_df', True)
        optimize_ram = kwargs.get("optimize_ram", True)
        self.kwargs['optimize_ram'] = optimize_ram

        ram_per_cpu = kwargs.get("ram_per_cpu", None)

        if optimize_ram:
            dump_rate = kwargs.get(
                "dump_rate", int(len(indices_for_iter) / 10))
            dump_rate = np.max([dump_rate, 1])

        # Set temporary directory for population batches
        temp_directory = os.path.abspath(
            kwargs.get("temp_directory", "batches"))
        self.kwargs['temp_directory'] = temp_directory

        kw = self.kwargs.copy()

        # Create temporary directory if it doesn't exist
        # Built to handle MPI
        if self.comm is None:
            if not os.path.exists(temp_directory):
                os.makedirs(temp_directory)
        else:
            if self.rank == 0:
                if not os.path.exists(temp_directory):
                    os.makedirs(temp_directory)

        filenames = []

        for j, index in enumerate(indices_for_iter):

            if kwargs.get('from_hdf', False):
                #generator
                binary = self.manager.from_hdf(index, restore=True).pop()
            else:
                binary = self.manager.generate(index=index, **self.kwargs)
            binary.properties = self.population_properties

            with warnings.catch_warnings(record=True) as w:
                try:
                    binary.evolve()
                except Exception:
                    #print(traceback.print_exc())
                    binary.event = 'FAILED'
                    binary.traceback = traceback.format_exc()
                if len(w) > 0:
                    warnings.simplefilter("always")
                    binary.warning_message = [x.message for x in w]

            if breakdown_to_df:
                self.manager.breakdown_to_df(binary, **self.kwargs)

            # Save data at some frequency
            if (optimize_ram
                    and j % dump_rate == 0 and j != 0 and ram_per_cpu is None):

                # Create filenames for each batch
                if(self.comm is None):
                    path = os.path.join(temp_directory, f"{j}_evolution.batch")
                else:
                    path = os.path.join(temp_directory,
                                        f"{j}_evolution.batch.{self.rank}")

                # save binaries to disk for RAM optimization
                self.manager.save(path, mode="w", **kw)
                filenames.append(path)

                if(breakdown_to_df):
                    self.manager.clear_dfs()
                else:
                    self.manager.remove(self.manager.binaries.copy())

            # Check to see if used memory is greater than 99% of allowed memory
            # rss gives memory usage in bytes, so divide by 2^30 to get GBs
            elif (optimize_ram and ram_per_cpu is not None
                  and psutil.Process().memory_info().rss / (1024**3)
                  >= 0.9 * ram_per_cpu):

                if(self.comm is None):
                    path = os.path.join(temp_directory, f"{j}_evolution.batch")
                else:
                    path = os.path.join(temp_directory,
                                        f"{j}_evolution.batch.{self.rank}")

                # save binaries to disk for RAM optimization
                self.manager.save(path, mode="w", **kw)
                filenames.append(path)

                if(breakdown_to_df):
                    self.manager.clear_dfs()
                else:
                    self.manager.remove(self.manager.binaries.copy())


        # handling case if dump rate is not multiple of population size
        if ((len(self.manager.binaries) != 0
                or len(self.manager.history_dfs) != 0)
                and optimize_ram):

            if(self.comm is None):
                path = os.path.join(temp_directory, "leftover_evolution.batch")
            else:
                path = os.path.join(temp_directory,
                                    f"leftover_evolution.batch.{self.rank}")

            # save binaries to disk for RAM optimization
            self.manager.save(path, mode="w", **kw)
            filenames.append(path)

            if(breakdown_to_df):
                self.manager.clear_dfs()
            else:
                self.manager.remove(self.manager.binaries.copy())

        if optimize_ram:
            # combining files
            if self.comm is None:
                self.combine_saved_files(os.path.join(temp_directory,
                                                      "evolution.combined"),
                                         filenames, mode = "w")
            else:
                self.combine_saved_files(
                    os.path.join(temp_directory,
                                 f"evolution.combined.{self.rank}"),
                    filenames, mode = "w")

        else:
            if self.comm is None:
                self.manager.save(os.path.join(temp_directory,
                                               "evolution.combined"),
                                  mode='w',
                                  **kwargs)
            else:
                self.manager.save(
                    os.path.join(temp_directory,
                                 f"evolution.combined.{self.rank}"),
                    mode='w', **kwargs)

    def save(self, save_path, mode='a', **kwargs):
        """Save BinaryPopulation to hdf file."""
        optimize_ram = self.kwargs['optimize_ram']
        temp_directory = self.kwargs['temp_directory']

        if self.comm is None:
            if optimize_ram:
                os.rename(os.path.join(temp_directory, "evolution.combined"),
                          save_path)
            else:
                self.manager.save(save_path, mode=mode, **kwargs)
        else:
            absolute_filepath = os.path.abspath(save_path)
            dir_name = os.path.dirname(absolute_filepath)
            file_name = os.path.basename(absolute_filepath)

            if os.path.isdir(absolute_filepath):
                file_name = 'backup_save_pop_data.h5'
                file_path = os.path.join(dir_name, file_name)
                warnings.warn('The provided path is a directory - saving '
                              'to {0} instead.'.format(file_path), Warning)

            self.comm.Barrier()

            if self.rank == 0:

                file_name = os.path.basename(absolute_filepath)
                tmp_files = [os.path.join(
                    self.kwargs["temp_directory"], f"evolution.combined.{i}")
                             for i in range(self.size)]


                self.combine_saved_files(absolute_filepath, tmp_files, mode = mode)

            else:
                return

    def make_temp_fname(self):
        """Get a valid filename for the temporary file."""
        temp_directory = self.kwargs['temp_directory']
        return os.path.join(temp_directory, f"evolution.combined.{self.rank}")
        # return os.path.join(dir_name, '.tmp{}_'.format(rank) + file_name)

    def combine_saved_files(self, absolute_filepath, file_names, mode = "a"):
        """Combine various temporary files in a given folder."""
        dir_name = os.path.dirname(absolute_filepath)

        history_cols = pd.read_hdf(file_names[0], key='history').columns
        oneline_cols = pd.read_hdf(file_names[0], key='oneline').columns

        history_tmp = pd.read_hdf(file_names[0], key='history')

        history_min_itemsize = {key: val for key, val in
                                HISTORY_MIN_ITEMSIZE.items()
                                if key in history_cols}
        oneline_min_itemsize = {key: val for key, val in
                                ONELINE_MIN_ITEMSIZE.items()
                                if key in oneline_cols}

        with pd.HDFStore(absolute_filepath, mode = mode) as store:
            for f in file_names:
                # strings itemsize set by first append max value,
                # which may not be largest string
                try:
                    store.append('history', pd.read_hdf(f, key='history'),
                                 min_itemsize=history_min_itemsize)
                    store.append('oneline', pd.read_hdf(f, key='oneline'),
                                 min_itemsize=oneline_min_itemsize)
                    os.remove(f)
                except Exception:
                    print(traceback.format_exc(), flush=True)

    def close(self):
        """Close loaded h5 files from SimulationProperties."""
        self.population_properties.close()
        self.manager.close()

    def __getstate__(self):
        """Prepare the BinaryPopulation to be 'pickled'."""
        # In order to be generally picklable, we need to discard the
        # communicator object before trying.
        d = self.__dict__
        d["comm"] = None
        prop = d['population_properties']
        if prop.steps_loaded:
            prop.close()
        return d

    def __iter__(self):
        """Iterate the binaries."""
        return iter(self.manager)

    def __getitem__(self, key):
        """Get the k-th binary."""
        return self.manager[key]

    def __len__(self):
        """Get the number of binaries in the population."""
        return len(self.manager)

    def __repr__(self):
        """Report key properties of the object."""
        s = "<{}.{} at {}>\n".format(
            self.__class__.__module__, self.__class__.__name__, hex(id(self))
        )
        for key, arg in self.kwargs.items():
            s += "{}: {}\n".format(key, arg)
        return s


class PopulationManager:
    """Manage a population of binaries."""

    def __init__(self, file_name=None, **kwargs):
        """Initialize a PopulationManager instance."""
        self.kwargs = kwargs.copy()
        self.binaries = []
        self.indices = []
        self.history_dfs = []
        self.oneline_dfs = []

        self.binary_generator = BinaryGenerator(**kwargs)
        self.entropy = self.binary_generator.entropy

        if file_name:
            self.store = pd.HDFStore(file_name, mode='r',)
            atexit.register(lambda: PopulationManager.close(self))

    def close(self):
        """Close the HDF5 file."""
        if hasattr(self, 'store'):
            self.store.close()

    def append(self, binary):
        """Add a binary instance internaly."""
        if isinstance(binary, (list, np.ndarray)):
            self.indices.append([b.index for b in binary])
            self.binaries.extend(list(binary))
        elif isinstance(binary, BinaryStar):
            self.indices.append(binary.index)
            self.binaries.append(binary)
        else:
            raise ValueError('Must be BinaryStar or list of BinaryStars')

    def remove(self, binary):
        """Remove a binary instance."""
        if isinstance(binary, (list, np.ndarray)):
            for b in binary:
                self.binaries.remove(b)
                self.indices.remove(b.index)

        elif isinstance(binary, BinaryStar):
            self.binaries.remove(binary)
            self.indices.remove(binary.index)
        else:
            raise ValueError('Must be BinaryStar or list of BinaryStars')

    def clear_dfs(self):
        """Remove all dfs."""
        self.history_dfs = []
        self.oneline_dfs = []

    def breakdown_to_df(self, binary, **kwargs):
        """Breakdown to a pandas DataFrame.

        Breakdown a binary into more convenient data type, store it, and
        remove the BinaryStar instance from self.

        """
        try:
            history = binary.to_df(**kwargs)
            self.history_dfs.append(history)
            oneline = binary.to_oneline_df(**kwargs)
            self.oneline_dfs.append(oneline)
            self.remove(binary)
        except Exception as err:
            print("Error during breakdown of {0}:\n{1}".
                  format(str(binary), err))

    def to_df(self, selection_function=None, **kwargs):
        """Convert all binaries to dataframe."""
        if len(self.binaries) == 0 and len(self.history_dfs) == 0:
            return
        is_callable = callable(selection_function)
        holder = []
        if len(self.binaries) > 0:
            for binary in self.binaries:
                if not is_callable or (is_callable
                                       and selection_function(binary)):
                    holder.append(binary.to_df(**kwargs))
        elif len(self.history_dfs) > 0:
            holder.extend(self.history_dfs)

        if len(holder) > 0:
            return pd.concat(holder, axis=0, ignore_index=False)

    def to_oneline_df(self, selection_function=None, **kwargs):
        """Convert all binaries to oneline dataframe."""
        if len(self.binaries) == 0 and len(self.oneline_dfs) == 0:
            return
        is_callable = callable(selection_function)
        holder = []
        if len(self.binaries) > 0:
            for binary in self.binaries:
                if not is_callable or (is_callable
                                       and selection_function(binary)):
                    holder.append(binary.to_oneline_df(**kwargs))
        elif len(self.oneline_dfs) > 0:
            holder.extend(self.oneline_dfs)

        if len(holder) > 0:
            return pd.concat(holder, axis=0, ignore_index=False)

    def find_failed(self,):
        """Find any failed binaries in the population."""
        if len(self) > 0:
            return [b for b in self if b.event == 'FAILED']
        elif len(self.history_dfs) > 0:
            failed_dfs = [f for f in self.history_dfs
                          if f['event'].iloc[-1] == 'FAILED']
            if bool(failed_dfs):
                return pd.concat(failed_dfs, axis=0, ignore_index=False)
            else:
                return failed_dfs

    def generate(self, **kwargs):
        """Generate a binary by drawing from the binary_generator.

        This can be a callable or a generator.
        """
        binary = self.binary_generator.draw_initial_binary(**kwargs)
        self.append(binary)
        return binary

    def from_hdf(self, indices, where=None, restore=False):
        """Load a BinaryStar instance from an hdf file of a saved population.

        Parameters
        ----------
        indices : int, list
            Selects the binaries to load.
        where : str
            Query performed on disk to select history and oneline DataFrames.
        restore : bool
            Restore binaries back to initial conditions.

        """
        if where is None:
            query_str = 'index==indices'
        else:
            query_str = str(where)

        hist = self.store.select(key='history', where=query_str)
        oneline = self.store.select(key='oneline', where=query_str)

        binary_holder = []
        for i in np.unique(hist.index):
            binary = BinaryStar.from_df(
                hist.loc[i],
                extra_columns=self.kwargs.get('extra_columns', []))

            # if the binary has failed
            if bool(oneline.loc[[i]]['FAILED'].values[-1]):
                setattr(binary, 'event', 'FAILED')

            bin_scalars = self.kwargs.get('scalar_names', [])
            s1_scalars = self.kwargs.get(
                'S1_kwargs', {}).get('scalar_names', [])
            s2_scalars = self.kwargs.get(
                'S2_kwargs', {}).get('scalar_names', [])
            if any([bool(bin_scalars), bool(s1_scalars), bool(s2_scalars)]):
                oline_bin = BinaryStar.from_oneline_df(
                        oneline.loc[[i]],
                        extra_columns=self.kwargs.get('extra_columns', [])
                )
                for name in bin_scalars:
                    val = getattr(oline_bin, name)
                    setattr(binary, name, val)
                for name in s1_scalars:
                    val = getattr(oline_bin.star_1, name)
                    setattr(binary.star_1, name, val)
                for name in s2_scalars:
                    val = getattr(oline_bin.star_2, name)
                    setattr(binary.star_2, name, val)

                del oline_bin

            binary_holder.append(binary)
            self.append(binary)

        if restore:
            [b.restore() for b in binary_holder]

        return binary_holder

    def save(self, fname, mode='a', **kwargs):
        """Save binaries to an hdf file using pandas HDFStore.

        Any object dtype columns not parsed by infer_objects() is converted to
        a string.

        Parameters
        ----------
        fname : str
            Name of hdf file saved.
        **kwargs

        Returns
        -------
        None
        """
        def set_dtypes_oneline(data):
            """Change the dtypes for consistency in saving."""
            for col in data.columns:
                if 'natal_kick_array' in col:
                    # First convert None's to NaN's
                    # data[col][data[col] == 'None'] = np.nan
                    data.loc[data[col] == 'None', col] = np.nan

                    # Next, convert dtype
                    data = data.astype({col: np.float64})

            return data

        with pd.HDFStore(fname, mode=mode) as store:

            history_df = self.to_df(**kwargs)

            # Set dtypes for saving
            history_df = history_df.infer_objects()

            object_to_str = {name: 'str' for i, name
                             in enumerate(history_df.columns)
                             if history_df.iloc[:, i].dtype == 'object'}
            history_df = history_df.astype(object_to_str)

            store.append('history', history_df, data_columns=True)

            online_df = self.to_oneline_df(**kwargs)
            online_df = online_df.infer_objects()
            object_to_str = {name: 'str' for i, name
                             in enumerate(online_df.columns)
                             if online_df.iloc[:, i].dtype == 'object'}
            online_df = online_df.astype(object_to_str)
            online_df = set_dtypes_oneline(online_df)

            store.append('oneline', online_df, data_columns=True)

#         store = pd.HDFStore(fname, mode=mode)

#         history_df = self.to_df(**kwargs)

#         # Set dtypes for saving
#         history_df = history_df.infer_objects()

#         # TODO: move these data type conversions into to_df, to_oneline_df
#         object_to_str = {name:'str'
#                          for i, name in enumerate(history_df.columns)
#                         if history_df.iloc[:,i].dtype == 'object'}
#         history_df = history_df.astype(object_to_str)

#         store.append('history', history_df, data_columns=True)

#         online_df = self.to_oneline_df(**kwargs)
#         online_df = online_df.infer_objects()
#         object_to_str = {name:'str'
#                          for i, name in enumerate(online_df.columns)
#                         if online_df.iloc[:,i].dtype == 'object'}
#         online_df = online_df.astype(object_to_str)
#         online_df = set_dtypes_oneline(online_df)

#         store.append('oneline', online_df, data_columns=True)

#         store.close()

    def __getitem__(self, key):
        """Return the key-th binary."""
        return self.binaries[key]

    def __iter__(self):
        """Iterate the binaries in the population."""
        return iter(self.binaries)

    def __len__(self):
        """Return the number of binaries in the population."""
        return len(self.binaries)

    def __bool__(self):
        """Evaluate as True if binaries have been appended."""
        return len(self) > 0


class BinaryGenerator:
    """Generate binaries to be included in a BinaryPopulation."""

    def __init__(self, sampler=generate_independent_samples,
                 RNG=None, **kwargs):
        """Initialize the BinaryGenerator instance."""
        self._num_gen = 0
        if RNG is None:
            self.RNG = np.random.default_rng()
            self.entropy = self.RNG.bit_generator._seed_seq.entropy
        else:
            assert isinstance(RNG, np.random.Generator)
            self.RNG = RNG
            self.entropy = self.RNG.bit_generator._seed_seq.entropy

        self.kwargs = kwargs.copy()
        self.sampler = sampler
        self.star_formation = kwargs.get('star_formation', 'burst')
        self.binary_fraction =  kwargs.get('binary_fraction', 1)
    def reset_rng(self):
        """Reset the RNG with the stored entropy."""
        self._num_gen = 0
        self.RNG = self.get_original_rng()

    def get_original_rng(self):
        """Return the random-number generator."""
        seq = np.random.SeedSequence(entropy=self.entropy)
        RNG = np.random.default_rng(seed=seq)
        return RNG

    def get_binary_by_iter(self, n=1, **kwargs):
        """Get the nth binary as if n calls were made to draw_intial_binary."""
        RNG = self.get_original_rng()
        original_num_gen = self._num_gen
        if n != 0:  # draw n-1 samples with new original RNG
            self.sampler(
                number_of_binaries=int(n-1), RNG=RNG, **kwargs)
        binary = self.draw_initial_binary(index=n, RNG=RNG, **kwargs)
        self._num_gen = original_num_gen
        return binary

    def draw_initial_samples(self, orbital_scheme='separation', **kwargs):
        """Generate all random varibles."""
        binary_fraction = self.kwargs.get('binary_fraction', 1)
        if not ('RNG' in kwargs.keys()):
            kwargs['RNG'] = self.RNG
        # a, e, M_1, M_2, P
        sampler_output = self.sampler(orbital_scheme, **kwargs)
        if orbital_scheme == 'separation':
            separation, eccentricity, m1, m2 = sampler_output
            orbital_period = orbital_period_from_separation(separation, m1, m2)
        elif orbital_scheme == 'period':
            orbital_period, eccentricity, m1, m2 = sampler_output
            separation = orbital_separation_from_period(orbital_period, m1, m2)
        else:
            raise ValueError("Allowed orbital schemes: separation or period.")

        # formation times
        N_binaries = len(orbital_period)
        formation_times = get_formation_times(N_binaries, **kwargs)

        # indices
        indices = np.arange(self._num_gen, self._num_gen+N_binaries, 1)

        output_dict = {
            'binary_index': indices,
            'binary_fraction':binary_fraction,
            'time': formation_times,
            'separation': separation,
            'eccentricity': eccentricity,
            'orbital_period': orbital_period,
            'S1_mass': m1,
            'S2_mass': m2,
        }
        self._num_gen += N_binaries
        return output_dict

    def draw_initial_binary(self, **kwargs):
        """Return a binary for evolution in a population.

        Parameters
        ----------
        index : int
            Sets binary index. Defaults to number of generated binaries.

        Returns
        -------
        binary : BinaryStar

        """
        sampler_kwargs = kwargs.copy()
        sampler_kwargs['number_of_binaries'] = 1
        sampler_kwargs['RNG'] = kwargs.get('RNG', self.RNG)
        output = self.draw_initial_samples(**sampler_kwargs)

        default_index = output['binary_index'].item()
        #Eirini's comments:
        # Randomly generated variables


        if self.RNG.uniform() < self.binary_fraction:
            formation_time = output['time'].item()
            separation = output['separation'].item()
            orbital_period = output['orbital_period'].item()
            eccentricity = output['eccentricity'].item()
            m1 = output['S1_mass'].item()
            m2 = output['S2_mass'].item()
            Z_div_Zsun = kwargs.get('metallicity', 1.)
            zams_table = {1.: 2.703e-01,
                          0.1: 2.511e-01,
                          0.01: 2.492e-01,
                          0.001: 2.49e-01,
                          0.0001: 2.49e-01}
            Y = zams_table[Z_div_Zsun]
            Z = Z_div_Zsun*Zsun
            X = 1. - Z - Y

            binary_params = dict(
                index=kwargs.get('index', default_index),
                time=formation_time,
                state="detached",
                event="ZAMS",
                separation=separation,
                orbital_period=orbital_period,
                eccentricity=eccentricity,
            )
            star1_params = dict(
                mass=m1,
                state="H-rich_Core_H_burning",
                metallicity=Z,
                center_h1=X,
                center_he4=Y,
            )
            star2_params = dict(
                mass=m2,
                state="H-rich_Core_H_burning",
                metallicity=Z,
                center_h1=X,
                center_he4=Y,
            )
        #If binary_fraction not default a initially single star binary is created.
        else:
            formation_time = output['time'].item()
            separation = output['separation'].item()
            orbital_period = output['orbital_period'].item()
            eccentricity = output['eccentricity'].item()
            m1 = output['S1_mass'].item()
            m2 = output['S2_mass'].item()
            Z_div_Zsun = kwargs.get('metallicity', 1.)
            zams_table = {1.: 2.703e-01,
                          0.1: 2.511e-01,
                          0.01: 2.492e-01,
                          0.001: 2.49e-01,
                          0.0001: 2.49e-01}
            Y = zams_table[Z_div_Zsun]
            Z = Z_div_Zsun*Zsun
            X = 1. - Z - Y

            binary_params = dict(
                index=kwargs.get('index', default_index),
                time=formation_time,
                state="initially_single_star",
                event="ZAMS",
                separation=separation,
                orbital_period=orbital_period,
                eccentricity=eccentricity,
            )
            star1_params = dict(
                mass=m1,
                state="H-rich_Core_H_burning",
                metallicity=Z,
                center_h1=X,
                center_he4=Y,
            )
            star2_params = dict(
                mass=m2,
                state="BH",
                metallicity=Z,
                center_h1=X,
                center_he4=Y,
            )
        #do all of the above but with state = "initial_single star"
        #in this case the second star can be a compact object.


        binary = BinaryStar(**binary_params,
                            star_1=SingleStar(**star1_params),
                            star_2=SingleStar(**star2_params))
        return binary

    def __repr__(self,):
        """Report key properties of the BinaryGenerator instance."""
        s = "<{}.{} at {}>\n".format(
            self.__class__.__module__, self.__class__.__name__, hex(id(self))
        )
        s += 'RNG: {}\n'.format(repr(self.RNG))
        s += 'entropy: {}\n'.format(self.entropy)
        for key, arg in self.kwargs.items():
            s += "{}: {}\n".format(key, arg)
        return s
