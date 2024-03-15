"""Processing of population files.

This module contains classes and functions to process population files.
Population files are HDF5 files containing the history and oneline dataframes
of a population of binary systems. The history dataframe contains the detailed
evolution of each binary system, while the oneline dataframe contains the
final state of each binary system.

Classes
-------
PopulationRunner
    A class to handle the evolution of binary populations.

Population
    A class to handle population files.
History
    A class to handle the history dataframe of a population file.
Oneline
    A class to handle the oneline dataframe of a population file.

TransientPopulation
    A class to handle transient populations.

Rates
    A class to handle the cosmic rates of in a population file.
"""

__authors__ = [
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Kyle Akira Rocha <kylerocha2024@u.northwestern.edu>",
    "Monica Gallegos-Garcia <monicagallegosgarcia2024@u.northwestern.edu>",
    "Max Briel <max.briel@unige.ch>",
]

import warnings
import numpy as np
import pandas as pd
from tqdm import tqdm
import os
from matplotlib import pyplot as plt

from posydon.utils.constants import Zsun
from posydon.popsyn.io import binarypop_kwargs_from_ini
from posydon.popsyn.normalized_pop_mass import initial_total_underlying_mass
import posydon.visualization.plot_pop as plot_pop
from posydon.utils.common_functions import convert_metallicity_to_string

from astropy.cosmology import Planck15 as cosmology
from astropy import constants as const

from posydon.popsyn.rate_calculation import (
    get_shell_comoving_volume,
    get_comoving_distance_from_redshift,
    get_cosmic_time_from_redshift,
    redshift_from_cosmic_time_interpolator,
    DEFAULT_MODEL,
    get_redshift_bin_edges,
    get_redshift_bin_centers,
)

from posydon.popsyn.star_formation_history import (
    star_formation_rate,
    SFR_Z_fraction_at_given_redshift,
)

from posydon.popsyn.binarypopulation import (
    BinaryPopulation,
    HISTORY_MIN_ITEMSIZE,
    ONELINE_MIN_ITEMSIZE,
)


###############################################################################

parameter_array = [
    "number_of_binaries",
    "binary_fraction_scheme",
    "binary_fraction_const",
    "star_formation",
    "max_simulation_time",
    "primary_mass_scheme",
    "primary_mass_min",
    "primary_mass_max",
    "secondary_mass_scheme",
    "secondary_mass_min",
    "secondary_mass_max",
    "orbital_scheme",
    "orbital_period_scheme",
    "orbital_period_min",
    "orbital_period_max",
    "eccentricity_scheme",
]


class PopulationRunner:
    """A class to handle the evolution of binary populations.

    Attributes
    ----------
    pop_params : dict
        The parameters of the population population.
    solar_metallicities : list of float
        The metallicities of the populations in solar units.
    binary_populations : list of BinaryPopulation
        The binary populations.
    verbose : bool
        If `True`, print additional information.

    Methods
    -------
    evolve()
        Evolve the binary populations.
    merge_parallel_runs(pop)
        Merge the parallel runs of the population.
    """

    def __init__(self, path_to_ini, verbose=False):
        """Initialize the binary populations from an ini file.

        Parameters
        ----------
        path_to_ini : str
            The path to the ini file.
        verbose : bool, optional
            If True, print additional information. Default is False.

        Raises
        ------
        ValueError
            If the provided `path_to_ini` does not have a '.ini' extension.

        Notes
        -----
        This method initializes the `PopulationRunner` object by reading the binary population parameters from an ini file specified by `path_to_ini`.
        It also allows for optional verbosity, which, when set to True, enables additional information to be printed in operation of the class.

        The `path_to_ini` should be a valid path to an ini file. If the provided `path_to_ini` does not have a '.ini' extension, a `ValueError` is raised.

        Examples
        --------
        >>> sp = PopulationRunner('/path/to/ini/file.ini')
        >>> sp = PopulationRunner('/path/to/ini/file.ini', verbose=True)

        """
        if ".ini" not in path_to_ini:
            raise ValueError("You did not provide a valid path_to_ini!")
        else:
            self.pop_params = binarypop_kwargs_from_ini(path_to_ini)
            self.solar_metallicities = self.pop_params["metallicity"]
            self.verbose = verbose
            if not isinstance(self.solar_metallicities, list):
                self.solar_metallicities = [self.solar_metallicities]

            self.binary_populations = []
            for MET in self.solar_metallicities:
                ini_kw = binarypop_kwargs_from_ini(path_to_ini)
                ini_kw["metallicity"] = MET
                ini_kw["temp_directory"] = (
                    convert_metallicity_to_string(MET)
                    + "_Zsun_"
                    + ini_kw["temp_directory"]
                )
                self.binary_populations.append(BinaryPopulation(**ini_kw))

    def evolve(self):
        """Evolve the binary populations.

        This method is responsible for evolving the binary populations. It iterates over each population
        in the `binary_populations` list and calls the `evolve` method on each population. After evolving
        each population, it merges them using the `merge_parallel_runs` method.

        Note:
            The `merge_parallel_runs` method is called only if the `comm` attribute of the population is `None`.

        """
        for pop in self.binary_populations:
            pop.evolve()
            if pop.comm is None:
                self.merge_parallel_runs(pop)

    def merge_parallel_runs(self, pop):
        """Merge the parallel runs of the population.

        Parameters
        ----------
        pop : BinaryPopulation
            The binary population whose files have to be merged.

        """
        if os.path.exists(convert_metallicity_to_string(pop.metallicity) + "_Zsun_population.h5"):
            raise FileExistsError(
                f"{convert_metallicity_to_string(pop.metallicity)}_Zsun_population.h5 already exists!\n"
                +"Files were not merged. You can use PopulationRunner.merge_parallel_runs() to merge the files manually."
            )
                
        path_to_batch = pop.kwargs["temp_directory"]
        tmp_files = [
            os.path.join(path_to_batch, f)
            for f in os.listdir(path_to_batch)
            if os.path.isfile(os.path.join(path_to_batch, f))
        ]
        pop.combine_saved_files(
            convert_metallicity_to_string(pop.metallicity) + "_Zsun_population.h5",
            tmp_files,
        )
        # remove files
        if len(os.listdir(path_to_batch)) == 0:
            os.rmdir(path_to_batch)


# TODO: I don't know if the merge_populations function is still needed
def merge_populations(populations, filename, verbose=False):
    """merges multiple populations into a single population file

    Parameters
    ----------
    populations : list of Population
        The populations to merge together
    filename : str
        The name of the population file to create
    verbose : bool
        If `True`, print additional information.

    Notes
    -----
    This function merges multiple populations into a single population file. It iterates over each population
    in the `populations` list and calls the `merge_parallel_runs` method on each population.
    It writes the history and oneline dataframes from the input files to the given output file.

    """

    history_cols = populations[0].history.columns
    oneline_cols = populations[0].oneline.columns

    history_min_itemsize = {
        key: val for key, val in HISTORY_MIN_ITEMSIZE.items() if key in history_cols
    }
    oneline_min_itemsize = {
        key: val for key, val in ONELINE_MIN_ITEMSIZE.items() if key in oneline_cols
    }

    # open new file
    with pd.HDFStore(filename, mode="a") as store:
        prev = 0
        for pop in populations:
            # write history
            for i in tqdm(
                range(0, len(pop), 10000), total=len(pop) // 10000, disable=not verbose
            ):
                tmp_df = pop.history[i : i + 10000]
                tmp_df.index = tmp_df.index + prev
                store.append(
                    "history",
                    tmp_df,
                    format="table",
                    data_columns=True,
                    min_itemsize=history_min_itemsize,
                )

            if verbose:
                print(f"History for {pop.filename} written to population file!")

            # write oneline
            for i in tqdm(
                range(0, len(pop), 10000), total=len(pop) // 10000, disable=not verbose
            ):
                tmp_df = pop.oneline[i : i + 10000]
                tmp_df.index = tmp_df.index + prev
                store.append(
                    "oneline",
                    tmp_df,
                    format="table",
                    data_columns=True,
                    min_itemsize=oneline_min_itemsize,
                )

            if verbose:
                print(f"Oneline for {pop.filename} written to population file!")

            prev += len(pop)

        # merge mass_per_metallicity
        tmp_df = pd.concat([pop.mass_per_metallicity for pop in populations])
        mass_per_metallicity = tmp_df.groupby(tmp_df.index).sum()

        store.put("mass_per_metallicity", mass_per_metallicity)

    # add ini parameters from the first population
    populations[0]._save_ini_params(filename)

    print(f"Populations merged into {filename}!")


##################
# Helper classes #
##################


class History:
    """A class to handle the history dataframe of a population file.

    This class provides methods to handle the history dataframe of a population file.
    It allows accessing and manipulating the history table based on various keys and conditions.

    Attributes
    ----------
    filename : str
        The path to the population file.
    verbose : bool
        If `True`, print additional information.
    chunksize : int
        The chunksize to use when reading the history file.
    lengths : pd.DataFrame
        The number of rows of each binary in the history dataframe.
    number_of_systems : int
        The number of systems in the history dataframe.
    columns : list of str
        The columns of the history dataframe.
    indices : np.ndarray
        The binary indices of the history dataframe.
    """

    def __init__(self, filename, verbose=False, chunksize=10000):
        """Initialise the history dataframe.

        This class is used to handle the history dataframe of a population file.
        On initialisation, the history_lengths are calculated and stored in the population file,
        if not present in the file. The history dataframe is not loaded into memory.

        Parameters
        ----------
        filename : str
            The path to the population file.
        verbose : bool
            If `True`, print additional information.
        chunksize : int (default=10000)
            The chunksize to use when reading the history file.
        """
        self.filename = filename
        self.verbose = verbose
        self.chunksize = chunksize
        self.lengths = None
        self.number_of_systems = None
        self.columns = None
        if not os.path.exists(filename):
            raise FileNotFoundError(f"{filename} does not exist!")

        # add history_lengths
        with pd.HDFStore(filename, mode="a") as store:
            # get the history lengths from the file
            if "/history_lengths" in store.keys():
                self.lengths = store["history_lengths"]
            else:
                if self.verbose:
                    print(
                        "history_lengths not found in population file. Calculating history lengths..."
                    )
                history_events = store.select_column("history", "index")
                self.lengths = pd.DataFrame(
                    history_events.groupby(history_events).count()
                )
                if self.verbose:
                    print("Storing history lengths in population file!")
                store.put("history_lengths", pd.DataFrame(self.lengths), format="table")

                del history_events

            self.columns = store.select("history", start=0, stop=0).columns

        self.indices = self.lengths.index.to_numpy()
        self.number_of_systems = len(self.lengths)

    def __getitem__(self, key):
        """Return the history table based on the provided key.


        Parameters
        ----------
        key : int, list of int, np.ndarray, or str
            The key to use for indexing the history dataframe.

        Returns
        -------
        pd.DataFrame
            The history table based on the provided key.

        Raises
        ------
        ValueError
            If the key type is invalid or if the column name(s) are not valid.

        Examples
        --------
        # Get a single row by index
        >>> population[0]
        Returns the history table for the row with index 0.

        # Get multiple rows by index
        >>> population[[0, 1, 2]]
        Returns the history table for the rows with indices 0, 1, and 2.

        # Get rows based on a boolean mask
        >>> mask = population['age'] > 30
        >>> population[mask]
        Returns the history table for the rows where the 'age' column is greater than 30.

        # Get a specific column
        >>> population['time']
        Returns the 'age' column from the history table.

        # Get multiple columns
        >>> population[['S1_mass', 'S2_mass']]
        Returns the 'S1_mass' and 'S2_mass' columns from the history table.
        """
        if isinstance(key, slice):
            if key.start is None:
                pre = 0
            else:
                pre = key.start
            if key.stop is None:
                chunk = self.number_of_systems
            else:
                chunk = key.stop - pre
            indices = list(range(pre, pre + chunk))
            return pd.read_hdf(self.filename, key="history", where='index in indices')       
        
        elif isinstance(key, int):
            return pd.read_hdf(self.filename, where="index == key", key="history")
        elif isinstance(key, list) and all(isinstance(x, int) for x in key):
            with pd.HDFStore(self.filename, mode="r") as store:
                return store.select("history", where="index in key")
        elif isinstance(key, np.ndarray) and (key.dtype == int):
            indices = key.tolist()
            with pd.HDFStore(self.filename, mode="r") as store:
                return store.select("history", where="index in indices")
        elif (isinstance(key, np.ndarray) and key.dtype == bool) or (isinstance(key, pd.DataFrame) and all(key.dtypes == bool)):
            out_df = pd.DataFrame()
            for i in range(0, len(key), self.chunksize):
                tmp_df = pd.read_hdf(
                    self.filename, key="history", start=i, stop=i + self.chunksize
                )[key[i : i + self.chunksize]]
                out_df = pd.concat([out_df, tmp_df])
            return out_df
        elif isinstance(key, str):
            if key in self.columns:
                return pd.read_hdf(self.filename, key="history", columns=[key])
            else:
                raise ValueError(f"{key} is not a valid column name!")
        elif isinstance(key, list) and all(isinstance(x, str) for x in key):
            if all(x in self.columns for x in key):
                return pd.read_hdf(self.filename, key="history", columns=key)
            else:
                raise ValueError(f"Not all columns in {key} are valid column names!")
        else:
            raise ValueError("Invalid key type!")

    def __len__(self):
        """Return the number of rows in the history table

        Returns
        -------
        int
            The number of rows in the history table.
        """
        return np.sum(self.lengths.values)

    def head(self, n=10):
        """Return the first n rows of the history table

        Parameters
        ----------
        n : int, optional
            The number of rows to return. Default is 10.

        Returns
        -------
        pandas.DataFrame
            The first n rows of the history table.

        """
        return pd.read_hdf(self.filename, key="history", start=0, stop=n)

    def tail(self, n=10):
        """Return the last n rows of the history table.

        Parameters:
        -----------
        n : int, optional
            Number of rows to return. Default is 10.

        Returns:
        --------
        pandas.DataFrame
            The last n rows of the history table.

        """
        return pd.read_hdf(self.filename, key="history", start=-n)

    def __repr__(self):
        """Return a string representation of the object.

        Returns:
            str: A string representation of the object.
        """
        return pd.read_hdf(self.filename, key="history").__repr__()

    def _repr_html_(self):
        """Return the HTML representation of the history dataframe.

        This method reads the history data from an HDF file and returns
        the HTML representation of the data using the `_repr_html_` method of the
        pandas DataFrame.

        Returns
        -------
        str
            The HTML representation of the history dataframe.
        """
        return pd.read_hdf(self.filename, key="history")._repr_html_()


    def select(self, where=None, start=None, stop=None, columns=None):
        """Select a subset of the history table based on the given conditions.

        This method allows you to query and retrieve a subset of data from the history table
        stored in an HDFStore file. You can specify conditions using the `where` parameter,
        which is a string representing the query condition to apply to the data. You can also
        specify the starting and ending indices of the data to select using the `start` and
        `stop` parameters. Additionally, you can choose to select specific columns by providing
        a list of column names using the `columns` parameter.

        Parameters
        ----------
        where : str, optional
            A string representing the query condition to apply to the data.
        start : int, optional
            The starting index of the data to select.
        stop : int, optional
            The ending index of the data to select.
        columns : list, optional
            A list of column names to select.

        Returns
        -------
        pandas.DataFrame
            The selected data as a DataFrame.
        """
        with pd.HDFStore(self.filename, mode="r") as store:
            return store.select(
                "history", where=where, start=start, stop=stop, columns=columns
            )


class Oneline:
    """A class to handle the oneline dataframe of a population file.

    The `Oneline` class provides methods to manipulate and retrieve data from the oneline dataframe of a population file.

    Attributes
    ----------
    filename : str
        The path to the population file.
    verbose : bool
        If `True`, print additional information.
    chunksize : int
        The chunksize to use when reading the oneline file.
    number_of_systems : int
        The number of systems in the oneline dataframe.
    indices : np.ndarray
        The binary indices of the oneline dataframe.
    columns : list of str

    Methods
    -------
    head(self, n=10)
        Get the first n rows of the oneline table.
    tail(self, n=10)
        Get the last n rows of the oneline table.
    select(self, where=None, start=None, stop=None, columns=None)
        Select a subset of the oneline table based on the given conditions.
    __init__(self, filename, verbose=False, chunksize=10000)
        Initialize the Oneline object.
    __getitem__(self, key)
        Return the oneline table based on the provided key.
    __len__(self)
        Return the number of systems in the oneline table.
    __repr__(self)
        Get a string representation of the oneline table.
    _repr_html_(self)
        Get an HTML representation of the oneline table.
    iterator(self)
        Get an iterator over the oneline table.
    iterbinaries(self)
        Get an iterator over the binary systems in the oneline table.

    """

    def __init__(self, filename, verbose=False, chunksize=10000):
        """Initialize a Oneline class instance.

        Parameters
        ----------
        filename : str
            The path to the HDFStore file containing the Oneline population data.
        verbose : bool, optional
            If True, print additional information during initialization. Default is False.
        chunksize : int, optional
            The number of rows to read from the HDFStore file at a time. Default is 10000.
        """
        self.filename = filename
        self.verbose = verbose
        self.chunksize = chunksize
        self.number_of_systems = None
        self.indices = None
        if not os.path.exists(filename):
            raise FileNotFoundError(f"{filename} does not exist!")

        with pd.HDFStore(filename, mode="r") as store:
            self.indices = store.select_column("oneline", "index").to_numpy()
            self.columns = store.select("oneline", start=0, stop=0).columns

        self.number_of_systems = len(self.indices)

    def __getitem__(self, key):
        """Get a subset of the oneline table based on the given key.

        Parameters
        ----------
        key : slice, int, list, np.ndarray, str
            The key to select the subset of the oneline table.

        Returns
        -------
        pd.DataFrame
            The subset of the oneline table.

        Raises
        ------
        ValueError
            If the key is of invalid type or contains invalid values.

        Examples
        --------
        # Get a slice of the oneline table
        >>> subset = population[10:20]

        # Get a single row from the oneline table
        >>> row = population[5]

        # Get multiple rows from the oneline table using a list of indices
        >>> rows = population[[1, 3, 5]]

        # Get rows from the oneline table using a boolean array
        >>> mask = population['age'] > 30
        >>> filtered_rows = population[mask]

        # Get a specific column from the oneline table
        >>> column = population['age']

        # Get multiple columns from the oneline table using a list of column names
        >>> columns = population[['age', 'gender']]
        """
        if isinstance(key, slice):
            if key.start is None:
                pre = 0
            else:
                pre = key.start
            if key.stop is None:
                chunk = self.number_of_systems
            else:
                chunk = key.stop - pre
            return pd.read_hdf(
                self.filename, key="oneline", start=pre, stop=pre + chunk
            )            
        elif isinstance(key, int):
            return pd.read_hdf(self.filename, where="index == key", key="oneline")
        elif isinstance(key, list) and all(isinstance(x, int) for x in key):
            return pd.read_hdf(self.filename, where="index in key", key="oneline")
        elif isinstance(key, np.ndarray) and (key.dtype == int):
            indices = key.tolist()
            return pd.read_hdf(self.filename, where="index in indices", key="oneline")
        elif isinstance(key, list) and all(isinstance(x, float) for x in key):
            raise ValueError("elements in list are not integers! Try casting to int.")
        elif isinstance(key, pd.DataFrame) and all(key.dtypes == bool):
            indices = self.indices[key.to_numpy().flatten()].tolist()
            return pd.read_hdf(self.filename, where="index in indices", key="oneline")
        elif isinstance(key, np.ndarray) and key.dtype == bool:
            indices = self.indices[key].tolist()
            return pd.read_hdf(self.filename, where="index in indices", key="oneline")
        elif isinstance(key, str):
            if key in self.columns:
                return pd.read_hdf(self.filename, key="oneline", columns=[key])
            else:
                raise ValueError(f"{key} is not a valid column!")
        elif isinstance(key, list) and all(isinstance(x, str) for x in key):
            if all(x in self.columns for x in key):
                return pd.read_hdf(self.filename, key="oneline", columns=key)
            else:
                raise ValueError(f"Not all columns in {key} are valid column names!")
        else:
            raise ValueError("Invalid key type!")

    def __len__(self):
        """
        Get the number of systems in the oneline table.

        Returns
        -------
        int
            The number of systems in the oneline table.
        """
        return self.number_of_systems

    def head(self, n=10):
        """Get the first n rows of the oneline table.

        Parameters
        ----------
        n : int, optional
            The number of rows to return. Default is 10.

        Returns
        -------
        pd.DataFrame
            The first n rows of the oneline table.
        """
        return pd.read_hdf(self.filename, key="oneline", start=0, stop=n)

    def tail(self, n=10):
        """
        Get the last n rows of the oneline table.

        Parameters
        ----------
        n : int, optional
            The number of rows to return. Default is 10.

        Returns
        -------
        pd.DataFrame
            The last n rows of the oneline table.
        """
        return pd.read_hdf(self.filename, key="oneline", start=-n)

    def __repr__(self):
        """
        Get a string representation of the oneline table.

        Returns
        -------
        str
            The string representation of the oneline table.
        """
        return pd.read_hdf(self.filename, key="oneline").__repr__()

    def _repr_html_(self):
        """
        Get an HTML representation of the oneline table.

        Returns
        -------
        str
            The HTML representation of the oneline table.
        """
        return pd.read_hdf(self.filename, key="oneline")._repr_html_()


    def select(self, where=None, start=None, stop=None, columns=None):
        """Select a subset of the oneline table based on the given conditions.

        This method allows you to filter and extract a subset of rows from the oneline table stored in an HDF file.
        You can specify conditions to filter the rows, define the range of rows to select, and choose specific columns to include in the subset.

        Parameters
        ----------
        where : str, optional
            A condition to filter the rows of the oneline table. Default is None.
        start : int, optional
            The starting index of the subset. Default is None.
        stop : int, optional
            The ending index of the subset. Default is None.
        columns : list, optional
            The column names to include in the subset. Default is None.

        Returns
        -------
        pd.DataFrame
            The selected subset of the oneline table.

        Examples
        --------
        # Select rows based on a condition
        >>> df = Oneline.select(where="S2_mass_i > 30")

        # Select rows from index 10 to 20
        >>> df = Oneline.select(start=10, stop=20)

        # Select specific columns
        >>> df = Oneline.select(columns=['S1_mass_i', 'S1_mass_f'])
        """
        with pd.HDFStore(self.filename, mode="r") as store:
            return store.select(
                "oneline", where=where, start=start, stop=stop, columns=columns
            )


class PopulationIO:
    """A class to handle the input/output of population files.

    This class provides methods to load and save population files in HDF5 format.
    It also includes methods to load and save metadata and ini parameters.

    Attributes
    ----------
    mass_per_metallicity : pandas.DataFrame
        A DataFrame containing mass per metallicity data.
    ini_params : dict
        A dictionary containing some ini parameters, described in parameter_array.

    """

    def __init__(self):
        self.verbose = False        


    def _load_metadata(self, filename):
        """Load the metadata from the file.

        Parameters
        ----------
        filename : str
            The name of the file to load the metadata from.

        Raises
        ------
        ValueError
            If the filename does not contain '.h5' extension.

        """
        if ".h5" not in filename:
            raise ValueError(
                f"{filename} does not contain .h5 in the se.\n Is this a valid population file?"
            )

        self._load_ini_params(filename)
        self._load_mass_per_metallicity(filename)

    def _save_mass_per_metallicity(self, filename):
        """Save the mass per metallicity data to the file.

        Parameters
        ----------
        filename : str
            The name of the file to save the mass per metallicity data to.

        """
        with pd.HDFStore(filename, mode="a") as store:
            store.put("mass_per_metallicity", self.mass_per_metallicity)
            if self.verbose:
                print("mass_per_metallicity table written to population file!")

    def _load_mass_per_metallicity(self, filename):
        """Load the mass per metallicity data from the file.

        Parameters
        ----------
        filename : str
            The name of the file to load the mass per metallicity data from.

        """
        with pd.HDFStore(filename, mode="r") as store:
            self.mass_per_metallicity = store["mass_per_metallicity"]
            if self.verbose:
                print("mass_per_metallicity table read from population file!")

    def _save_ini_params(self, filename):
        """Save the ini parameters to the file.

        Parameters
        ----------
        filename : str
            The name of the file to save the ini parameters to.

        """
        with pd.HDFStore(filename, mode="a") as store:
            # write ini parameters to file
            tmp_df = pd.DataFrame()
            for c in parameter_array:
                tmp_df[c] = [self.ini_params[c]]
            store.put("ini_parameters", tmp_df)

    def _load_ini_params(self, filename):
        """Load the ini parameters from the file.

        The values loaded from file are stored in the parameter_array.

        Parameters
        ----------
        filename : str
            The name of the file to load the ini parameters from.

        """
        # load ini parameters
        with pd.HDFStore(filename,mode="r",) as store:
            tmp_df = store["ini_parameters"]
            self.ini_params = {}
            for c in parameter_array:
                self.ini_params[c] = tmp_df[c][0]


##########################
# Main interface classes #
##########################


class Population(PopulationIO):
    """A class to handle population files.

    This class provides methods to handle population files. It includes methods to read and write population files,
    as well as methods to access and manipulate the history and oneline dataframes.

    Attributes
    ----------
    history : History
        The history dataframe of the population.
    oneline : Oneline
        The oneline dataframe of the population.
    formation_channels : pd.DataFrame
        The formation channels dataframe of the population.
    ini_params : dict
        The parameters from the ini file used to create the population.
    mass_per_metallicity : pd.DataFrame
        The mass per metallicity dataframe of the population.
    solar_metallicities : np.ndarray
        The solar metallicities of the population.
    metallicities : np.ndarray
        The metallicities of the population.
    indices : np.ndarray
        The indices of the binaries in the population.
    verbose : bool
        If `True`, print additional information.
    chunksize : int
        The chunksize to use when reading the population file.
    filename : str
        The path to the population file.
    number_of_systems : int
        The number of systems in the population.
    history_lengths : pd.DataFrame
        The number of rows of each binary in the history dataframe. The index is the binary index.

    Methods
    -------
    export_selection(selection, filename, chunksize)
        Export a selection of the population to a new file.
    calculate_formation_channels(mt_history=False)
        Calculate the formation channels of the population.
    __len__()
        Get the number of systems in the population.
    columns()
        Get the columns of the population.
    create_transient_population(func, transient_name, oneline_cols=None, hist_cols=None)
        Create a transient population using a given function.
    """

    def __init__(
        self, filename, metallicity=None, ini_file=None, verbose=False, chunksize=1000
    ):
        """Initialize the Population object.

        The Population object is initialised by creating History and Oneline objects,
        which refer back to the population file (filename). The formation channels are also
        linked to the population file, if present. The mass per metallicity data is loaded
        from the file, if present, or calculated and saved to the file if not present.

        If the mass per metallicity data is not present in the file, you can provide a metallicity
        and the ini file (used to create the population) to calculate and save the mass per metallicity data.
        You only need to do this once for a given population file. However, all systems in the population
        will be given the same metallicity.

        Parameters
        -----------
        filename : str
            The path to the population file
        metallicity : float, optional
            The metallicity of the population in solar units.
        ini_file : str, optional
            The path to the ini file used to create the population.
        verbose : bool, optional
            If `True`, print additional information.
        chunksize : int, optional
            The chunksize to use when reading the population file.

        Raises
        ------
        ValueError
            If the provided filename does not contain '.h5' extension.
            If the population file does not contain a history table.
            If the population file does not contain an oneline table.
            If the population file does not contain an ini_parameters table.
            If the population file does not contain a mass_per_metallicity table and no metallicity for the file was given.

        Examples
        --------
        # When the population file contains a mass_per_metallicity table
        >>> pop = Population('/path/to/population_file.h5')

        # When the population file does not contain a mass_per_metallicity table
        >>> pop = Population('/path/to/population_file.h5', metallicity=0.02, ini_file='/path/to/ini_file.ini')
        """

        self.filename = filename
        self.verbose = verbose
        self.chunksize = chunksize

        self.mass_per_metallicity = None
        self.number_of_systems = None
        self.history_lengths = None

        # if the user provided a single string instead of a list of strings
        if not (".h5" in filename):
            raise ValueError(
                f"{filename} does not contain .h5 in the se.\n Is this a valid population file?"
            )

        # read the population file
        with pd.HDFStore(filename, mode="r") as store:
            keys = store.keys()

        # check if pop contains history
        if "/history" not in keys:
            raise ValueError(f"{filename} does not contain a history table!")
        else:
            self.history = History(filename, self.verbose, self.chunksize)

        # check if pop contains oneline
        if "/oneline" not in keys:
            raise ValueError(f"{filename} does not contain an oneline table!")
        else:
            self.oneline = Oneline(filename, self.verbose, self.chunksize)

        # check if formation channels are present
        if "/formation_channels" not in keys:
            if self.verbose:
                warnings.warn(f"{filename} does not contain formation channels!")
            self._formation_channels = None
        else:
            self._formation_channels = pd.read_hdf(
                self.filename, key="formation_channels"
            )

        # if an ini file is given, read the parameters from the ini file
        if ini_file is not None:
            self.ini_params = binarypop_kwargs_from_ini(ini_file)
            self._save_ini_params(filename)
            self._load_ini_params(filename)
        else:
            if "/ini_parameters" not in keys:
                raise ValueError(
                    f"{filename} does not contain an ini_parameters table!"
                )
            else:
                self._load_ini_params(filename)

        # check if pop contains mass_per_metallicity table
        if "/mass_per_metallicity" in keys and metallicity is None:
            self._load_mass_per_metallicity(filename)
            self.solar_metallicities = self.mass_per_metallicity.index.to_numpy()
            self.metallicities = self.solar_metallicities * Zsun
        elif metallicity is None:
            raise ValueError(
                f"{filename} does not contain a mass_per_metallicity table and no metallicity for the file was given!"
            )

        # calculate the metallicity information. This assumes the metallicity is for the whole file!
        if metallicity is not None and ini_file is not None:
            if "/mass_per_metallicity" in keys:
                warnings.warn(
                    f"{filename} already contains a mass_per_metallicity table. Overwriting the table!"
                )

            simulated_mass = np.sum(self.oneline[["S1_mass_i", "S2_mass_i"]].to_numpy())
            underlying_mass = initial_total_underlying_mass(
                df=simulated_mass, **self.ini_params
            )[0]
            self.mass_per_metallicity = pd.DataFrame(
                index=[metallicity],
                data={
                    "simulated_mass": simulated_mass,
                    "underlying_mass": underlying_mass,
                    "number_of_systems": len(self.oneline),
                },
            )

            self._save_mass_per_metallicity(filename)
            self.solar_metallicities = self.mass_per_metallicity.index.to_numpy()
            self.metallicities = self.solar_metallicities * Zsun

        elif metallicity is not None and ini_file is None:
            raise ValueError(
                f"{filename} does not contain a mass_per_metallicity table and no ini file was given!"
            )

        # add number of systems
        self.history_lengths = self.history.lengths
        self.number_of_systems = self.oneline.number_of_systems
        self.indices = self.history.indices

    def export_selection(self, selection, filename, overwrite=True, history_chunksize=1):
        """Export a selection of the population to a new file

        This method exports a selection of systems from the population to a new file.
        The selected systems are specified by their indices in the population.

        By default the selected systems will overwrite the file if it already exists.
        If overwrite is set to False, the selected systems will be appended to 
        the existing file if it already exists and the indices will be shifted 
        based on the current length of data in the file that is being appended to.


        Parameters
        ----------
        selection : list of int
            The indices of the systems to export.
        filename : str
            The name of the export file to create or append to.
        chunksize : int, optional
            The number of systems to export at a time. Default is 1.

        Raises
        ------
        ValueError
            If the filename does not contain ".h5" extension.

        Warnings
        --------
        UserWarning
            If there is no "metallicity" column in the oneline dataframe and the population file contains multiple metallicities.

        Notes
        -----
        - The exported systems will be appended to the existing file if it already exists.
        - The indices of the exported systems will be shifted based on the current length of data in the file.
        - The "oneline" and "history" dataframes of the selected systems will be written to the file.
        - If available, the "formation_channels" dataframe and "history_lengths" dataframe of the selected systems will also be written to the file.
        - The "metallicity" column of the oneline dataframe will be added if it is not present, using the metallicity of the population file.
        - The "mass_per_metallicity" dataframe will be updated with the number of selected systems.

        Examples
        --------
        # Export systems with indices [0, 1, 2] to a new file named "selected.h5"
        >>> population.export_selection([0, 1, 2], "selected.h5")

        # Export systems with indices [3, 4, 5] to an existing file named "existing.h5"
        >>> population.export_selection([3, 4, 5], "existing.h5")

        # Export systems with indices [6, 7, 8] to a new file named "selected.h5" in chunks of 2
        >>> population.export_selection([6, 7, 8], "selected.h5", history_chunksize=2)
        """

        if not (".h5" in filename):
            raise ValueError(
                f"{filename} does not contain .h5 in the se.\n Is this a valid population file?"
            )

        history_cols = self.history.columns
        oneline_cols = self.oneline.columns

        history_min_itemsize = {
            key: val for key, val in HISTORY_MIN_ITEMSIZE.items() if key in history_cols
        }
        oneline_min_itemsize = {
            key: val for key, val in ONELINE_MIN_ITEMSIZE.items() if key in oneline_cols
        }
        if overwrite:
            mode = 'w'
        else:
            mode = 'a'
        with pd.HDFStore(filename, mode=mode) as store:
            # shift all new indices by the current length of data in the file
            last_index_in_file = 0

            if "/oneline" in store.keys():
                last_index_in_file = np.sort(store["oneline"].index)[-1]
            elif "/history" in store.keys():
                last_index_in_file = np.sort(store["history"].index)[-1]

            if "/history" in store.keys() and self.verbose:
                print("history in file. Appending to file")

            if "/oneline" in store.keys() and self.verbose:
                print("oneline in file. Appending to file")

            if "/formation_channels" in store.keys() and self.verbose:
                print("formation_channels in file. Appending to file")

            if "/history_lengths" in store.keys() and self.verbose:
                print("history_lengths in file. Appending to file")

            # TODO: I need to shift the indices of the binaries or should I reindex them?
            # since I'm storing the information, reindexing them should be fine.

            if last_index_in_file == 0:
                reindex = {
                    i: j
                    for i, j in zip(selection,
                                    np.arange(last_index_in_file, last_index_in_file + len(selection), 1),
                    )
                }
            else:
                reindex = {
                    i: j
                    for i, j in zip(selection,
                                    np.arange(last_index_in_file + 1, last_index_in_file + len(selection) + 1,1,),
                    )
                }

            if "metallicity" not in self.oneline.columns:
                warnings.warn(
                    "No metallicity column in oneline dataframe! Using the metallicity of the population file and adding it to the oneline."
                )
                if len(self.metallicities) > 1:
                    raise ValueError(
                        "The population file contains multiple metallicities. Please add a metallicity column to the oneline dataframe!"
                    )

            if self.verbose:
                print("Writing selected systems to population file...")

            # write oneline of selected systems
            for i in tqdm(
                range(0, len(selection), self.chunksize),
                total=len(selection) // self.chunksize,
                disable=not self.verbose,
            ):
                
                tmp_df = self.oneline[selection[i : i + self.chunksize]]
                if "metallicity" in tmp_df.columns:
                    tmp_df["metallicity"] = tmp_df["metallicity"].astype("float")
                else:
                    tmp_df["metallicity"] = self.metallicities[0]
                tmp_df.rename(index=reindex, inplace=True)
                store.append(
                    "oneline",
                    tmp_df,
                    format="table",
                    data_columns=True,
                    min_itemsize=oneline_min_itemsize,
                    index=False,
                )

            if self.verbose:
                print("Oneline: Done")

            # write history of selected systems
            for i in tqdm(
                range(0, len(selection), history_chunksize),
                total=len(selection) // history_chunksize,
                disable=not self.verbose,
            ):
                tmp_df = self.history[selection[i : i + history_chunksize]]
                tmp_df.rename(index=reindex, inplace=True)
                store.append(
                    "history",
                    tmp_df,
                    format="table",
                    data_columns=True,
                    min_itemsize=history_min_itemsize,
                    index=False,
                )

            if self.verbose:
                print("History: Done")

            # write formation channels of selected systems
            if self.formation_channels is not None:
                for i in tqdm(
                    range(0, len(selection), self.chunksize),
                    total=len(selection) // self.chunksize,
                    disable=not self.verbose,
                ):
                    tmp_df = self.formation_channels.loc[selection[i : i + self.chunksize]]
                    tmp_df.rename(index=reindex, inplace=True)
                    store.append(
                        "formation_channels",
                        tmp_df,
                        format="table",
                        data_columns=True,
                        min_itemsize={"channel_debug": 100, "channel": 100},
                    )

            ## METADATA

            # write the history lengths
            for i in tqdm(
                range(0, len(selection), self.chunksize),
                total=len(selection) // self.chunksize,
                disable=not self.verbose,
            ):
                tmp_df = self.history.lengths.loc[selection[i : i + self.chunksize]]
                tmp_df.rename(index=reindex, inplace=True)
                store.append(
                    "history_lengths", pd.DataFrame(tmp_df), format="table", index=False
                )

            # write mass_per_metallicity
            if "/mass_per_metallicity" in store.keys():
                self_mass = self.mass_per_metallicity
                self_mass["number_of_systems"] = len(selection)
                tmp_df = pd.concat([store["mass_per_metallicity"], self_mass])
                mass_per_metallicity = tmp_df.groupby(tmp_df.index).sum()
                store.put("mass_per_metallicity", mass_per_metallicity)

            else:
                self_mass = self.mass_per_metallicity
                self_mass["number_of_systems"] = len(selection)
                store.put("mass_per_metallicity", self_mass)

        # write ini parameters
        self._save_ini_params(filename)

    @property
    def formation_channels(self):
        """
        Retrieves the formation channels from the population file.

        Returns:
            pandas.DataFrame or None: The formation channels if available, otherwise None.
        """
        with pd.HDFStore(self.filename, mode="r") as store:
            if "/formation_channels" in store.keys():
                self._formation_channels = pd.read_hdf(
                    self.filename, key="formation_channels"
                )
            else:
                if self.verbose:
                    warnings.warn("No formation channels in the population file!")
                self._formation_channels = None

        return self._formation_channels

    def calculate_formation_channels(self, mt_history=False):
        """Calculate the formation channels of the population.

        mt_history is a boolean that determines if the detailed mass-transfer history
        from the HMS-HMS grid is included in the formation channels.


        Parameters
        ----------
        mt_history : bool, optional
            If `True`, include the mass-transfer history in the formation channels. Default is False.

        Raises
        ------
        ValueError
            If the mt_history_HMS_HMS column is not present in the oneline dataframe.
        """

        if self.verbose:
            print("Calculating formation channels...")

        # load the HMS-HMS interp class
        HMS_HMS_event_dict = {
            "initial_MT": "initial_MT",
            "stable_MT": "oRLO1",
            "no_MT": "None",
            "unstable_MT": "oCE1/oDoubleCE1",
        }

        unique_binary_indices = self.indices

        # check if formation channels already exist
        with pd.HDFStore(self.filename, mode="a") as store:
            if "/formation_channels" in store.keys():
                print("Formation channels already exist in the parsed population file!")
                print("Channels will be overwriten")
                del store["formation_channels"]

        def get_events(group):
            # for now, only append information for RLO1; unstable_MT information already exists
            if "oRLO1" in group["interp_class_HMS_HMS"].tolist():
                combined_events = (
                    group["event"].iloc[0] + "_" + group["interp_class_HMS_HMS"].iloc[0]
                )
                tmp = [combined_events]
                tmp.extend(group["event"].iloc[1:])
                combined_events = "_".join(tmp)
            else:
                combined_events = "_".join(group["event"])
            return pd.Series({"channel_debug": combined_events})

        def mt_history(row):
            if (
                pd.notna(row["mt_history_HMS_HMS"])
                and row["mt_history_HMS_HMS"] == "Stable contact phase"
            ):
                return row["channel"].replace("oRLO1", "oRLO1-contact")
            elif (
                pd.notna(row["mt_history_HMS_HMS"])
                and row["mt_history_HMS_HMS"] == "Stable reverse mass-transfer phase"
            ):
                return row["channel"].replace("oRLO1", "oRLO1-reverse")
            else:
                return row["channel"]

        previous = 0

        for i in tqdm(
            range(0, len(unique_binary_indices), self.chunksize),
            disable=not self.verbose,
        ):
            selection = unique_binary_indices[i : i + self.chunksize]

            # create the dataframe for the chunk
            df = pd.DataFrame(index=selection, columns=["channel_debug", "channel"])
            end = (
                previous
                + self.history_lengths.iloc[i : i + self.chunksize].sum().iloc[0]
            )

            # get the history of chunk events and transform the interp_class_HMS_HMS
            interp_class_HMS_HMS = self.oneline.select(
                start=i, stop=i + self.chunksize, columns=["interp_class_HMS_HMS"]
            )
            events = self.history.select(start=previous, stop=end, columns=["event"])

            mask = ~pd.isna(interp_class_HMS_HMS["interp_class_HMS_HMS"].values)
            interp_class_HMS_HMS.loc[mask, "interp_class_HMS_HMS"] = (
                interp_class_HMS_HMS[mask]
                .apply(lambda x: HMS_HMS_event_dict[x["interp_class_HMS_HMS"]], axis=1)
                .values
            )
            del mask

            previous = end
            # combine based on the index, this allows for an easier apply later
            merged = pd.merge(
                events.dropna(), interp_class_HMS_HMS, left_index=True, right_index=True
            )
            del events, interp_class_HMS_HMS

            merged.index.name = "binary_index"
            df["channel_debug"] = merged.groupby("binary_index").apply(get_events)
            del merged
            df["channel"] = (
                df["channel_debug"]
                .str.replace("_redirect_from_ZAMS", "")
                .str.replace("_redirect_from_CO_HMS_RLO", "")
                .str.replace("_redirect_from_CO_HeMS_RLO", "")
                .str.replace("_redirect_from_CO_HeMS", "")
                .str.replace("_CO_contact", "")
            )

            if mt_history:
                columns = self.oneline.columns
                if "mt_history_HMS_HMS" not in columns:
                    raise ValueError(
                        "mt_history_HMS_HMS not saved in the oneline dataframe!"
                    )
                else:
                    tmp_df = pd.DataFrame(
                        index=selection, columns=["channel", "mt_history_HMS_HMS"]
                    )
                    tmp_df["channel"] = df["channel"]
                    x = self.oneline.select(
                        start=i, stop=i + self.chunksize, columns=["mt_history_HMS_HMS"]
                    )
                    tmp_df["mt_history_HMS_HMS"] = x
                    df["channel"] = tmp_df.apply(mt_history, axis=1)
                    del tmp_df
                    del x

            self._write_formation_channels(self.filename, df)
            del df
            
        if self.verbose:
            print("formation_channels written to population file!")

    def _write_formation_channels(self, filename, df):
        """Write the formation channels to the population file
        
        This will append the formation channels to the population file, while restricting the maximum
        length of the channel and channel_debug columns to 100 characters.

        Parameters
        ----------
        filename : str
            The name of the file to write the formation channels to.
        df : pd.DataFrame
            The dataframe containing the formation channels.
        """
        str_length = 100
        
        with pd.HDFStore(filename, mode="a") as store:
            df['channel'] = df['channel'].str.slice(0, str_length)
            df['channel_debug'] = df['channel_debug'].str.slice(0, str_length)
            store.append(
                "formation_channels",
                df,
                format="table",
                data_columns=True,
                min_itemsize={"channel_debug": str_length, "channel": str_length},
            )

    def __len__(self):
        """Get the number of systems in the population.

        Returns
        -------
        int
            The number of systems in the population.

        """
        return self.number_of_systems

    @property
    def columns(self):
        """
        Returns a dictionary containing the column names of the history and oneline dataframes.

        Returns
        -------
        dict
            A dictionary with keys 'history' and 'oneline', where the values are the column names of the respective dataframes.
        """
        return {"history": self.history.columns, "oneline": self.oneline.columns}

    def create_transient_population(
        self, func, transient_name, oneline_cols=None, hist_cols=None
    ):
        """Given a function, create a TransientPopulation

        This method creates a transient population using the provided function.
        `func` is given the history, oneline, and formation channels dataframes as arguments.
        The function should return a dataframe containing the transient population, which needs to contain the columns 'time' and 'metallicity'.

        Processing is done in chunks to avoid memory issues and a pandas DataFrame is stored at `'/transients/transient_name'` in the population file.

        The creation of the transient population can be sped up by limiting the oneline_cols and hist_cols to only the columns needed for the function.
        If you do not provide these, all columns will be used, which can be slow for large populations.

        Parameters
        ----------
        func : function
            Function to apply to the parsed population to create a transient population.
            The function needs to take 3 arguments:
                - history_chunk : pd.DataFrame
                - oneline_chunk : pd.DataFrame
                - formation_channels_chunk : pd.DataFrame
            and return a pd.DataFrame containing the transient population, which needs to contain the columns 'time' and 'metallicity'.

        oneline_cols : list of str, optional
            Columns to extract from the oneline dataframe. Default is all columns.

        hist_cols : list of str, optional
            Columns to extract from the history dataframe. Default is all columns.

        Returns
        -------
        TransientPopulation
            A TransientPopulation object for interfacing with the transient population

        Raises
        ------
        ValueError
            If the transient population requires a time column or if the transient population contains duplicate columns.

        Examples
        --------
        See the tutorials for examples of how to use this method.

        """

        with pd.HDFStore(self.filename, mode="a") as store:
            if f"/transients/{transient_name}" in store.keys():
                print("overwriting transient population")
                del store["transients/" + transient_name]

        min_itemsize = {
            "channel": 100,
        }
        if hist_cols is not None:
            if "time" not in hist_cols:
                raise ValueError("The transient population requires a time column!")

            min_itemsize.update(
                {
                    key: val
                    for key, val in HISTORY_MIN_ITEMSIZE.items()
                    if key in hist_cols
                }
            )
        else:
            hist_cols = self.history.columns
            min_itemsize.update(HISTORY_MIN_ITEMSIZE)

        if oneline_cols is not None:
            min_itemsize.update(
                {
                    key: val
                    for key, val in ONELINE_MIN_ITEMSIZE.items()
                    if key in oneline_cols
                }
            )
        else:
            oneline_cols = self.oneline.columns
            min_itemsize.update(ONELINE_MIN_ITEMSIZE)

        # setup a mapping to the size of each history colummn
        history_lengths = self.history_lengths
        unique_binary_indices = self.indices

        previous = 0
        for i in tqdm(
            range(0, len(unique_binary_indices), self.chunksize),
            disable=not self.verbose,
        ):
            end = previous + history_lengths[i : i + self.chunksize].sum().iloc[0]

            oneline_chunk = self.oneline.select(
                start=i, stop=i + self.chunksize, columns=oneline_cols
            )

            history_chunk = self.history.select(
                start=previous, stop=end, columns=hist_cols
            )

            if self.formation_channels is not None:
                formation_channels_chunk = self.formation_channels[
                    i : i + self.chunksize
                ]
            else:
                formation_channels_chunk = None

            syn_df = func(history_chunk, oneline_chunk, formation_channels_chunk)

            if len(syn_df.columns) != len(syn_df.columns.unique()):
                raise ValueError("Transient population contains duplicate columns!")

            # filter out the columns in min_itemsize that are not in the dataframe
            min_itemsize = {
                key: val for key, val in min_itemsize.items() if key in syn_df.columns
            }

            with pd.HDFStore(self.filename, mode="a") as store:
                store.append(
                    "transients/" + transient_name,
                    syn_df,
                    format="table",
                    data_columns=True,
                    min_itemsize=min_itemsize,
                )

            previous = end

        synth_pop = TransientPopulation(
            self.filename, transient_name, verbose=self.verbose
        )
        return synth_pop

    def plot_binary_evolution(self, index):
        """Plot the binary evolution of a system

        This method is not currently implemented.
        """
        pass


class TransientPopulation(Population):
    """A class representing a population of transient events.

    This class allows you to calculate additional properties of the population,
    such as the efficiency of events per Msun for each solar metallicity, and to
    calculate the cosmic weights of the transient population.

    Attributes
    ----------
    population : pandas.DataFrame
        DataFrame containing the whole transient population.
    transient_name : str
        Name of the transient population.
    efficiency : pandas.DataFrame
        DataFrame containing the efficiency of events per Msun for each solar metallicity.
    columns : list
        List of columns in the transient population.


    Methods
    -------
    select(where=None, start=None, stop=None, columns=None)
        Select a subset of the transient population.
    get_efficiency_over_metallicity()
        Compute the efficiency of events per Msun for each solar metallicity.
    calculate_cosmic_weights(SFH_identifier, MODEL_in=None)
        Calculate the cosmic weights of the transient population.
    plot_efficiency_over_metallicity()
        Plot the efficiency of events per Msun for each solar metallicity.
    plot_delay_time_distribution()
        Plot the delay time distribution of the transient population.
    plot_popsyn_pver_grid_slice()
        Plot the transient population over the parameter space of a grid slice
    """

    def __init__(self, filename, transient_name, verbose=False):
        """Initialise the TransientPopulation object.

        This method initializes the TransientPopulation object by linking it to the population file.
        The transient population linked is located at '/transients/{transient_name}' in the population file.


        Parameters
        ----------
        filename : str
            The name of the file containing the population. The file should be in HDF5 format.
        transient_name : str
            The name of the transient population within the file.
        verbose : bool, optional
            If `True`, additional information will be printed during the initialization process.

        Raises
        ------
        ValueError
            If the specified transient population name is not found in the file.

        Notes
        -----
        The population data is stored in an HDF5 file format. The file should contain a group named '/transients' which
        holds all the transient populations. The specified transient population name should be a valid group name within
        '/transients'. If the transient population has associated efficiencies, they will be loaded as well.

        Examples
        --------
        >>> filename = 'population_data.h5'
        >>> transient_name = 'BBH'
        >>> population = TransientPopulation(filename, transient_name, verbose=True)
        """
        super().__init__(filename, verbose=verbose)

        with pd.HDFStore(self.filename, mode="r") as store:
            if "/transients/" + transient_name not in store.keys():
                raise ValueError(
                    f"{transient_name} is not a valid transient population in {filename}!"
                )

            self.transient_name = transient_name
            if "/transients/" + transient_name + "/efficiencies" in store.keys():
                self._load_efficiency(filename)

    @property
    def population(self):
        """Returns the entire transient population as a pandas DataFrame.

        This method retrieves the transient population data from a file and returns it as a pandas DataFrame.
        Please note that if the transient population is too large, it may consume a significant amount of memory.

        Returns
        -------
        pd.DataFrame
            A DataFrame containing the transient population data.
        """
        return pd.read_hdf(self.filename, key="transients/" + self.transient_name)

    def _load_efficiency(self, filename):
        """Load the efficiency from the file

        Parameters:
            filename (str): The path to the file containing the efficiency data.

        Returns:
            None

        Raises:
            None

        """
        with pd.HDFStore(filename, mode="r") as store:
            self.efficiency = store[
                "transients/" + self.transient_name + "/efficiencies"
            ]
            if self.verbose:
                print("Efficiency table read from population file!")

    def _save_efficiency(self, filename):
        """Save the efficiency to the file.

        Args:
            filename (str): The name of the file to save the efficiency to.

        Returns:
            None
        """
        with pd.HDFStore(filename, mode="a") as store:
            store.put(
                "transients/" + self.transient_name + "/efficiencies", self.efficiency
            )

    @property
    def columns(self):
        """Return the columns of the transient population.

        Returns:
            list: A list of column names in the transient population.
        """
        if not hasattr(self, "_columns"):
            with pd.HDFStore(self.filename, mode="r") as store:
                self._columns = store.select(
                    "transients/" + self.transient_name, start=0, stop=0
                ).columns
        return self._columns

    def select(self, where=None, start=None, stop=None, columns=None):
        """
        Select a subset of the transient population.

        This method allows you to filter and extract a subset of rows from the oneline table stored in an HDF file.
        You can specify conditions to filter the rows, define the range of rows to select, and choose specific columns to include in the subset.

        Parameters
        ----------
        where : str, optional
            A condition to filter the rows of the oneline table. Default is None.
        start : int, optional
            The starting index of the subset. Default is None.
        stop : int, optional
            The ending index of the subset. Default is None.
        columns : list, optional
            The column names to include in the subset. Default is None.

        Returns
        -------
        pd.DataFrame
            The selected subset of the oneline table.

        Examples
        --------
        # Select rows based on a condition
        >>> df = transpop.select(where="S2_mass_i > 30")

        # Select rows from index 10 to 20
        >>> df = transpop.select(start=10, stop=20)

        # Select specific columns
        >>> df = transpop.select(columns=['time', 'metallicity'])
        """
        return pd.read_hdf(
            self.filename,
            key="transients/" + self.transient_name,
            where=where,
            start=start,
            stop=stop,
            columns=columns,
        )

    def get_efficiency_over_metallicity(self):
        """
        Compute the efficiency of events per Msun for each solar_metallicities.

        This method calculates the efficiency of events per solar mass for each solar metallicity value in the transient population.
        It first checks if the efficiencies have already been computed and if so, overwrites them.
        Then, it iterates over each metallicity value and calculates the efficiency by dividing the count of events with the underlying stellar mass.
        The efficiencies are stored in a DataFrame with the metallicity values as the index and the 'total' column representing the efficiency for all channels.
        If the 'channel' column is present, it also computes the merger efficiency per channel and adds the results to the DataFrame.
        """
        if hasattr(self, "efficiency"):
            print("Efficiencies already computed! Overwriting them!")
            with pd.HDFStore(self.filename, mode="a") as store:
                if (
                    "/transients/" + self.transient_name + "/efficiencies"
                    in store.keys()
                ):
                    del store["transients/" + self.transient_name + "/efficiencies"]

        metallicities = self.mass_per_metallicity.index.to_numpy()
        efficiencies = []
        for met in metallicities:
            count = self.select(where="metallicity == {}".format(met)).shape[0]

            # just sums the number of events
            underlying_stellar_mass = self.mass_per_metallicity["underlying_mass"][met]

            eff = count / underlying_stellar_mass
            efficiencies.append(eff)
            print(f"Efficiency at Z={met:1.2E}: {eff:1.2E} Msun^-1")

        self.efficiency = pd.DataFrame(
            index=metallicities, data={"total": np.array(efficiencies)}
        )

        # if the channel column is present compute the merger efficiency per channel
        if "channel" in self.columns:
            channels = np.unique(self.select(columns=["channel"]).values)
            for ch in channels:
                efficiencies = []
                for met in metallicities:
                    count = np.sum(
                        (
                            self.select(
                                where="(metallicity == {})".format(met),
                                columns=["channel"],
                            )
                            == ch
                        ).values
                    )
                    if count > 0:
                        underlying_stellar_mass = self.mass_per_metallicity["underlying_mass"][
                            met
                        ]
                        eff = count / underlying_stellar_mass
                    else:
                        eff = 0
                    efficiencies.append(eff)
                self.efficiency[ch] = np.array(efficiencies)

        # save the efficiency
        self._save_efficiency(self.filename)

    def calculate_cosmic_weights(self, SFH_identifier, MODEL_in=None):
        """
        Calculate the cosmic weights of the transient population.

        This method calculates the cosmic weights of the transient population based on the provided star formation history identifier and model parameters.
        It performs various calculations and stores the results in an HDF5 file at the location '/transients/{transient_name}/rates/{SFH_identifier}'.
        This allows for multiple star formation histories to be used with the same transient population.

        The default MODEL parameters are used if none are provided, which come from the IllustrisTNG model.

        Parameters
        ----------
        SFH_identifier : str
            Identifier for the star formation history.
        MODEL_in : dict, optional
            Dictionary containing the model parameters. If not provided, the default model parameters will be used.

        Returns
        -------
        Rates
            An instance of the Rates class.

        Raises
        ------
        ValueError
            If a parameter name in MODEL_in is not valid.

        Notes
        -----
        This function calculates the cosmic weights of the transient population based on the provided star formation history
        identifier and model parameters. It performs various calculations and stores the results in an HDF5 file.

        The cosmic weights are computed for each event in the population, taking into account the metallicity, redshift,
        and birth time of the events. The weights are calculated using the provided model parameters and the underlying mass
        distribution.

        The calculated weights, along with the corresponding redshifts of the events, are stored in the HDF5 file for further analysis.
        These can be accessed using the Rates class.

        Examples
        --------
        >>> transient_population = TransientPopulation('filename.h5', 'transient_name')
        >>> transient_population.calculate_cosmic_weights('IllustrisTNG', MODEL_in=DEFAULT_MODEL)
        """

        # Set model to DEFAULT or provided MODEL parameters
        # Allows for partial model specification
        if MODEL_in is None:
            MODEL = DEFAULT_MODEL
        else:
            for key in MODEL_in:
                if key not in DEFAULT_MODEL:
                    raise ValueError(key + " is not a valid parameter name!")

            # write the DEFAULT_MODEL with updates parameters to self.MODEL.
            MODEL = DEFAULT_MODEL
            MODEL.update(MODEL_in)

        path_in_file = (
            "/transients/" + self.transient_name + "/rates/" + SFH_identifier + "/"
        )

        with pd.HDFStore(self.filename, mode="a") as store:
            if path_in_file + "MODEL" in store.keys():
                store.remove(path_in_file + "MODEL")
                if self.verbose:
                    print("Cosmic weights already computed! Overwriting them!")
                if path_in_file + "weights" in store.keys():
                    store.remove(path_in_file + "weights")
                if path_in_file + "z_events" in store.keys():
                    store.remove(path_in_file + "z_events")
                if path_in_file + "birth" in store.keys():
                    store.remove(path_in_file + "birth")

        self._write_MODEL_data(self.filename, path_in_file, MODEL)

        rates = Rates(
            self.filename, self.transient_name, SFH_identifier, verbose=self.verbose
        )

        z_birth = rates.centers_redshift_bins
        t_birth = get_cosmic_time_from_redshift(z_birth)
        nr_of_birth_bins = len(z_birth)
        # write birth to the population file
        with pd.HDFStore(self.filename, mode="a") as store:
            store.put(
                path_in_file + "birth", pd.DataFrame(data={"z": z_birth, "t": t_birth})
            )

        get_redshift_from_cosmic_time = redshift_from_cosmic_time_interpolator()
        indices = self.indices

        # sample the SFH for only the events that are within the Hubble time
        # I only need to sample the SFH at each metallicity and z_birth
        # Not for every event!
        SFR_at_z_birth = star_formation_rate(rates.MODEL["SFR"], z_birth)
        # get metallicity bin edges
        met_edges = rates.edges_metallicity_bins

        # get the fractional SFR at each metallicity and z_birth
        fSFR = SFR_Z_fraction_at_given_redshift(
            z_birth,
            rates.MODEL["SFR"],
            rates.MODEL["sigma_SFR"],
            met_edges,
            rates.MODEL["Z_max"],
            rates.MODEL["select_one_met"],
        )

        # simulated mass per given metallicity corrected for the unmodeled
        # single and binary stellar mass
        M_model = rates.mass_per_metallicity.loc[rates.centers_metallicity_bins / Zsun][
            "underlying_mass"
        ].values

        # speed of light
        c = const.c.to("Mpc/yr").value  # Mpc/yr

        # delta cosmic time bin
        deltaT = rates.MODEL["delta_t"] * 10**6  # yr

        for i in tqdm(
            range(0, len(indices), self.chunksize),
            desc="event loop",
            disable=not self.verbose,
        ):
            selected_indices = (
                self.select(start=i, stop=i + self.chunksize, columns=["index"])
                .index.to_numpy()
                .flatten()
            )
            if len(selected_indices) == 0:
                continue

            # selected_indices = indices[i:i+self.chunksize]
            delay_time = (
                self.select(
                    start=i, stop=i + self.chunksize, columns=["time"]
                ).to_numpy()
                * 1e-3
            )  # Gyr

            t_events = t_birth + delay_time
            hubble_time_mask = t_events <= cosmology.age(1e-08).value * 0.9999999

            # get the redshift of the events
            z_events = np.full(t_events.shape, np.nan)
            z_events[hubble_time_mask] = get_redshift_from_cosmic_time(
                t_events[hubble_time_mask]
            )

            D_c = get_comoving_distance_from_redshift(z_events)  # Mpc

            # the events have to be in solar metallicity
            met_events = (
                self.select(start=i, stop=i + self.chunksize, columns=["metallicity"])
                .to_numpy()
                .flatten()
                * Zsun
            )

            weights = np.zeros((len(met_events), nr_of_birth_bins))
            for i, met in enumerate(rates.centers_metallicity_bins):
                mask = met_events == met
                weights[mask, :] = (
                    4.0
                    * np.pi
                    * c
                    * D_c[mask] ** 2
                    * deltaT
                    * (fSFR[:, i] * SFR_at_z_birth)
                    / M_model[i]
                )  # yr^-1

            with pd.HDFStore(self.filename, mode="a") as store:
                store.append(
                    path_in_file + "weights",
                    pd.DataFrame(data=weights, index=selected_indices),
                    format="table",
                )
                store.append(
                    path_in_file + "z_events",
                    pd.DataFrame(data=z_events, index=selected_indices),
                    format="table",
                )
        return rates

    def plot_efficiency_over_metallicity(self, **kwargs):
        """
        Plot the efficiency over metallicity.

        Parameters
        ----------
        channel : bool, optional
            If True, plot the subchannels. Default is False.
        """
        if not hasattr(self, "efficiency"):
            raise ValueError(
                "First you need to compute the efficiency over metallicity!"
            )
        plot_pop.plot_merger_efficiency(
            self.efficiency.index.to_numpy() * Zsun, self.efficiency, **kwargs
        )

    def plot_delay_time_distribution(
        self, metallicity=None, ax=None, bins=100, color="black"
    ):
        """
        Plot the delay time distribution of the transient population.

        This method plots the delay time distribution of the transient population. If a specific metallicity is provided,
        the delay time distribution of the population at that metallicity will be plotted. Otherwise, the delay time distribution
        of the entire population will be plotted.

        Parameters
        ----------
        metallicity : float or None
            The metallicity value to select a specific population. If None, the delay time distribution of the entire population will be plotted.
        ax : matplotlib.axes.Axes or None
            The axes object to plot the distribution on. If None, a new figure and axes will be created.
        bins : int
            The number of bins to use for the histogram.
        color : str
            The color of the histogram.

        Raises
        ------
        ValueError
            If the specified metallicity is not present in the population.

        Notes
        -----
        - The delay time distribution is normalized by the total mass of the population if no metallicity is specified.
        Otherwise, it is normalized by the mass of the population at the specified metallicity.

        """
        if ax is None:
            fig, ax = plt.subplots()

        if metallicity is None:
            time = self.select(columns=["time"]).values
            time = time * 1e6  # yr
            h, bin_edges = np.histogram(time, bins=bins)
            h = h / np.diff(bin_edges) / self.mass_per_metallicity["underlying_mass"].sum()

        else:
            if not any(np.isclose(metallicity, self.solar_metallicities)):
                raise ValueError("The metallicity is not present in the population!")

            time = self.select(
                where="metallicity == {}".format(metallicity), columns=["time"]
            ).values
            time = time * 1e6  # yr
            h, bin_edges = np.histogram(time, bins=bins)
            h = (
                h
                / np.diff(bin_edges)
                / self.mass_per_metallicity["underlying_mass"][metallicity]
            )

        ax.step(bin_edges[:-1], h, where="post", color=color)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Time [yr]")
        ax.set_ylabel("Number of events/Msun/yr")

    def plot_popsyn_over_grid_slice(self, grid_type, met_Zsun, **kwargs):
        """
        Plot the transients over the grid slice.

        Parameters
        ----------
        grid_type : str
            The type of grid to plot.
        met_Zsun : float
            The metallicity of the Sun.
        **kwargs
            Additional keyword arguments to pass to the plot_pop.plot_popsyn_over_grid_slice function.

        """

        plot_pop.plot_popsyn_over_grid_slice(
            pop=self, grid_type=grid_type, met_Zsun=met_Zsun, **kwargs
        )

    def _write_MODEL_data(self, filename, path_in_file, MODEL):
        """
        Write the MODEL data to the HDFStore file.

        Parameters
        ----------
        filename : str
            The path to the HDFStore file.
        path_in_file : str
            The path within the HDFStore file to store the MODEL data.
        MODEL : dict
            The MODEL data to be stored.

        """
        with pd.HDFStore(filename, mode="a") as store:
            if MODEL["dlogZ"] is not None:
                store.put(path_in_file + "MODEL", pd.DataFrame(MODEL))
            else:
                store.put(path_in_file + "MODEL", pd.DataFrame(MODEL, index=[0]))
            if self.verbose:
                print("MODEL written to population file!")


class Rates(TransientPopulation):
    """Class representing rates of a transient population.

    Attributes
    ----------
    SFH_identifier : str
        The identifier for the star formation history.
    base_path : str
        The base path for accessing the rates data.
    MODEL : dict
        The model data for the star formation history.
    weights : pandas.DataFrame
        The weights of the transient population.
    z_birth : pandas.DataFrame
        The redshift of the birth bins.
    z_events : pandas.DataFrame
        The redshift of the events.
    intrinsic_rate_density : pandas.DataFrame
        The intrinsic rate density of the transient population.
    observable_population_names : list
        The names of the observable populations.
    edges_metallicity_bins : np.ndarray
        The edges of the metallicity bins of the star formation history.
    centers_metallicity_bins : np.ndarray
        The centers of the metallicity bins of the star formation history.
    edges_redshift_bins : np.ndarray
        The edges of the redshift bins of the star formation history.
    centers_redshift_bins : np.ndarray
        The centers of the redshift bins of the star formation history.

    Methods
    -------
    select_rate_slice(key, start=None, stop=None)
        Selects a slice of a rates dataframe at key.
    calculate_intrinsic_rate_density(mt_channels=False)
        Compute the intrinsic rate density of the transient population.
    calculate_observable_population(observable_func, observable_name)
        Calculate the observable population based on the provided function.
    observable_population(observable_name)
        Return the observable population based on the provided name.
    plot_hist_properties(prop, intrinsic=True, observable=None, bins=50, channel=None, **kwargs)
        Plot a histogram of a given property available in the transient population.
    """

    def __init__(self, filename, transient_name, SFH_identifier, verbose=False):
        """
        Initialize the Rates object.

        This method initializes a Rates object by linking it to the population file
        with the specified transient name and star formation history identifier.
        The path in the file is '/transients/{transient_name}/rates/{SFH_identifier}'.

        Parameters:
        -----------
        filename : str
            The path to the file containing the transient population data.
        transient_name : str
            The name of the transient.
        SFH_identifier : str
            The identifier for the star formation history.
        verbose : bool, optional
            Whether to print verbose output. Default is False.
        """

        super().__init__(filename, transient_name, verbose=verbose)
        self.SFH_identifier = SFH_identifier

        self.base_path = (
            "/transients/" + self.transient_name + "/rates/" + self.SFH_identifier + "/"
        )

        with pd.HDFStore(self.filename, mode="r") as store:
            if (
                "/transients/"
                + self.transient_name
                + "/rates/"
                + self.SFH_identifier
                + "/MODEL"
                not in store.keys()
            ):
                raise ValueError(
                    f"{self.SFH_identifier} is not a valid SFH_identifier in {filename}!"
                )

        # load in the SFH_model
        self._read_MODEL_data(self.filename)

    def _read_MODEL_data(self, filename):
        """
        Reads the MODEL data from the specified file.

        Parameters
        ----------
        filename : str
            The path to the file containing the MODEL data.
        """
        with pd.HDFStore(filename, mode="r") as store:
            tmp_df = store[self.base_path + "MODEL"]
            if len(tmp_df) > 1:
                self.MODEL = tmp_df.iloc[0].to_dict()
                self.MODEL["dlogZ"] = [tmp_df["dlogZ"].min(), tmp_df["dlogZ"].max()]
            else:
                self.MODEL = tmp_df.iloc[0].to_dict()

            if self.verbose:
                print("MODEL read from population file!")

    @property
    def weights(self):
        """
        Retrieves the weights from the HDFStore.

        The rows are indexed by the binary index of the events, while to columns are indexed by the redshift of the birth bins.

        Returns
        -------
        pandas.DataFrame
            The weights DataFrame.
        """
        with pd.HDFStore(self.filename, mode="r") as store:
            return store[self.base_path + "weights"]

    @property
    def z_birth(self):
        """
        Retrieves the 'birth' data from the HDFStore.

        The 'birth' DataFrame contains the redshift and age of the Universe of the birth bins with columns 'z' and 't'.

        Returns
        -------
        pandas.DataFrame
            The 'birth' DataFrame, which contains the redshift and age of the Universe of the birth bins.
        """
        with pd.HDFStore(self.filename, mode="r") as store:
            return store[self.base_path + "birth"]

    @property
    def z_events(self):
        """
        Returns the 'z_events' data from the HDFStore.

        The 'z_events' data contains the redshifts at which the events occur.
        The rows of the returned DataFrame are indexed by the binary index of the events.
        The columns of the returned DataFrame are indexed by the redshift of the birth bins.

        Returns
        -------
        pandas.DataFrame
            The 'z_events' data from the HDFStore.

        """
        with pd.HDFStore(self.filename, mode="r") as store:
            return store[self.base_path + "z_events"]

    def select_rate_slice(self, key, start=None, stop=None):
        """Selects a slice of a rates dataframe at key.

        This method allows you to select a slice in rows from the diffferent rates dataframes.
        The slice is selected based on the start and stop indices.
        The key specifies which rates dataframe to select, and must be one of ['weights', 'z_events', 'birth'].

        Parameters
        ----------
        key : str
            The key to select the slice from. Must be one of ['weights', 'z_events', 'birth'].
        start : int, optional
            The starting index of the slice. Defaults to None.
        stop : int, optional
            The ending index of the slice. Defaults to None.

        Returns
        -------
        pandas.DataFrame
            The selected slice of rates.

        Raises
        ------
        ValueError
            If the key is not one of ['weights', 'z_events', 'birth'].
        """
        if key not in ["weights", "z_events", "birth"]:
            raise ValueError("key not in [weights, z_events, birth]")

        with pd.HDFStore(self.filename, mode="r") as store:
            return store.select(self.base_path + key, start=start, stop=stop)

    def calculate_intrinsic_rate_density(self, mt_channels=False):
        """
        Compute the intrinsic rate density over redshift of the transient population.

        Besides returning the intrinsic rate density, this method also stores the results in the HDF5 file for further analysis.
        This can be accessed using the intrinsic_rate_density attribute of the Rates class.

        Parameters
        ----------
        mt_channels : bool, optional
            Flag indicating whether to calculate the intrinsic rate density for each channel separately. Default is False.

        Returns
        -------
        pandas.DataFrame
            DataFrame containing the intrinsic rate density values.
        """
        z_events = self.z_events.to_numpy()
        weights = self.weights.to_numpy()
        z_horizon = self.edges_redshift_bins
        n = len(z_horizon)

        if mt_channels:
            channels = self.select(columns=["channel"])
            unique_channels = np.unique(channels)
        else:
            unique_channels = []

        intrinsic_rate_density = pd.DataFrame(index=z_horizon[:-1], columns=["total"])

        normalisation = np.zeros(n - 1)

        for i in tqdm(range(1, n), total=n - 1, disable=not self.verbose):
            normalisation[i - 1] = get_shell_comoving_volume(
                z_horizon[i - 1], z_horizon[i], "infinite"
            )

        for i in tqdm(range(1, n), total=n - 1, disable=not self.verbose):
            mask = (z_events > z_horizon[i - 1]) & (z_events <= z_horizon[i])
            for ch in unique_channels:
                mask_ch = channels.to_numpy() == ch
                intrinsic_rate_density.loc[z_horizon[i - 1], ch] = (
                    np.nansum(weights[mask & mask_ch]) / normalisation[i - 1]
                )

            intrinsic_rate_density.loc[z_horizon[i - 1], "total"] = (
                np.nansum(weights[mask]) / normalisation[i - 1]
            )

        with pd.HDFStore(self.filename, mode="a") as store:
            store.put(self.base_path + "intrinsic_rate_density", intrinsic_rate_density)

        return intrinsic_rate_density

    def calculate_observable_population(self, observable_func, observable_name):
        """
        Calculate an observable population.

        The observable population is calculated based on the provided observable function.
        It should recalculate your weights based on an observability probability given certain transinet parameters.

        This function requires the following input chunks:
        1. transient_pop_chunk
        2. z_events_chunk
        3. weights_chunk
        The observable function should output a DataFrame with the same shape as the weights_chunk.

        The observable population is stored in the HDF5 file at
        the location '/transients/{transient_name}/rates/observable/{observable_name}'.

        Parameters
        ----------
        observable_func : function
            The observability function that takes the TransientPopulation as input.
        observable_name : str
            The name of the observable.

        Note
        ----
        - If the observable population already exists in the file, it will be overwritten.
        """

        with pd.HDFStore(self.filename, mode="a") as store:
            # remove the observable population if it already exists
            if (
                "/transients/"
                + self.transient_name
                + "/rates/observable/"
                + observable_name
                in store.keys()
            ):
                if self.verbose:
                    print("Overwriting observable population!")
                del store[
                    "transients/"
                    + self.transient_name
                    + "/rates/observable/"
                    + observable_name
                ]

        # loop over the transient population and calculate the new weights, while writing to the file
        for i in tqdm(
            range(0, len(self), self.chunksize),
            total=len(self) // self.chunksize,
            disable=not self.verbose,
        ):
            transient_pop_chunk = self.select(start=i, stop=i + self.chunksize)
            weights_chunk = self.select_rate_slice(
                "weights", start=i, stop=i + self.chunksize
            )
            z_events_chunk = self.select_rate_slice(
                "z_events", start=i, stop=i + self.chunksize
            )
            new_weights = observable_func(
                transient_pop_chunk, z_events_chunk, weights_chunk
            )

            with pd.HDFStore(self.filename, mode="a") as store:
                store.append(
                    "transients/"
                    + self.transient_name
                    + "/rates/observable/"
                    + observable_name,
                    new_weights,
                    format="table",
                )

    def observable_population(self, observable_name):
        """Return the observable population based on the provided name.

        This method returns the observable population based on the provided name,
        which can take a while to load if the population is large.
        It loads the observable population from '/transients/{transient_name}/rates/observable/{observable_name}' in the HDF5 file.

        Parameters
        ----------
        observable_name : str
            The name of the observable population to return.

        Returns
        -------
        pandas.DataFrame
            The observable population based on the provided name.
        """

        with pd.HDFStore(self.filename, mode="r") as store:
            if (
                "/transients/"
                + self.transient_name
                + "/rates/observable/"
                + observable_name
                not in store.keys()
            ):
                raise ValueError(
                    f"{observable_name} is not a valid observable population!"
                )
            else:
                return store[
                    "transients/"
                    + self.transient_name
                    + "/rates/observable/"
                    + observable_name
                ]

    @property
    def observable_population_names(self):
        """Return the names of the observable populations in the associated file.

        Returns
        -------
        list
            The names of the observable populations.
        """

        with pd.HDFStore(self.filename, mode="r") as store:
            return [
                key.split("/")[-1]
                for key in store.keys()
                if "/transients/" + self.transient_name + "/rates/observable/" in key
            ]

    @property
    def intrinsic_rate_density(self):
        """Return the intrinsic rate density of the transient population at the specified SFH_identifier and transient_name.

        The data is read from the HDF5 file at '/transients/{transient_name}/rates/{SFH_identifier}/intrinsic_rate_density'.

        Returns
        -------
        pandas.DataFrame
            The intrinsic rate density of the transient population.

        """
        with pd.HDFStore(self.filename, mode="r") as store:
            if self.base_path + "intrinsic_rate_density" not in store.keys():
                raise ValueError(
                    "First you need to compute the intrinsic rate density!"
                )
            else:
                return store[self.base_path + "intrinsic_rate_density"]

    def plot_hist_properties(
        self, prop, intrinsic=True, observable=None, bins=50, channel=None, **kwargs
    ):
        """Plot a histogram of a given property available in the transient population.

        This method plots a histogram of a given property available in the transient population.
        The property can be intrinsic or observable, and the histogram can be plotted for a specific channel if provided.

        Parameters
        ----------
        prop : str
            The property to plot the histogram for.
        intrinsic : bool, optional
            If True, plot the intrinsic property. Default is True.
        observable : str, optional
            The observable population name to plot the histogram for. Default is None.
        bins : int, optional
            The number of bins to use for the histogram. Default is 50.
        channel : str, optional
            The channel to plot the histogram for. Default is None.
            A channel column must be present in the transient population.
        **kwargs
            Additional keyword arguments to pass to the plot

        Raises
        ------
        ValueError
            If the specified property is not a valid property in the transient population.
            If the specified observable is not a valid observable population.
            If the specified channel is not present in the transient population.

        """

        if prop not in self.columns:
            raise ValueError(
                f"{prop} is not a valid property in the transient population!"
            )

        # get the property and its associated weights in the population.

        df = self.select(columns=[prop])
        df["property"] = df[prop]
        del df[prop]
        if intrinsic:
            df["intrinsic"] = np.sum(self.weights, axis=1)

        if observable is not None:
            with pd.HDFStore(self.filename, mode="r") as store:
                if (
                    "/transients/"
                    + self.transient_name
                    + "/rates/observable/"
                    + observable
                    not in store.keys()
                ):
                    raise ValueError(
                        f"{observable} is not a valid observable population!"
                    )
                else:
                    df["observable"] = np.sum(
                        store[
                            "transients/"
                            + self.transient_name
                            + "/rates/observable/"
                            + observable
                        ],
                        axis=1,
                    )
        if channel is not None:
            df["channel"] = self.select(columns=["channel"])
            df = df[df["channel"] == channel]
            if len(df) == 0:
                raise ValueError(
                    f"{channel} is not present in the transient population!"
                )
            plot_pop.plot_hist_properties(df, bins=bins, **kwargs)

        else:
            # plot the histogram using plot_pop.plot_hist_properties
            plot_pop.plot_hist_properties(df, bins=bins, **kwargs)

    @property
    def edges_metallicity_bins(self):
        """Return the edges of the metallicity bins.

        Returns
        -------
        array float
            Returns the edges of all metallicity bins. We assume metallicities
            were binned in log-space.

        """
        met_val = np.log10(self.centers_metallicity_bins)
        bin_met = np.zeros(len(met_val) + 1)
        # if more than one metallicty bin
        if len(met_val) > 1:
            bin_met[0] = met_val[0] - (met_val[1] - met_val[0]) / 2.0
            bin_met[-1] = met_val[-1] + (met_val[-1] - met_val[-2]) / 2.0
            bin_met[1:-1] = met_val[:-1] + (met_val[1:] - met_val[:-1]) / 2.0
        # one metallicty bin
        elif len(met_val) == 1:
            if isinstance(self.MODEL["dlogZ"], float):
                bin_met[0] = met_val[0] - self.MODEL["dlogZ"] / 2.0
                bin_met[-1] = met_val[0] + self.MODEL["dlogZ"] / 2.0
            elif isinstance(self.MODEL["dlogZ"], list) or isinstance(
                self.MODEL["dlogZ"], np.array
            ):
                bin_met[0] = self.MODEL["dlogZ"][0]
                bin_met[-1] = self.MODEL["dlogZ"][1]

        return 10**bin_met

    @property
    def centers_metallicity_bins(self):
        """Return the centers of the metallicity bins.

        Returns
        -------
        array float
            Returns sampled metallicities of the population. This corresponds
            to the center of each metallicity bin.
        """
        return np.sort(self.metallicities)

    @property
    def edges_redshift_bins(self):
        """Compute redshift bin edges.

        Returns
        -------
        array floats
            We devide the cosmic time history of the Universe in equally spaced
            bins of cosmic time of self.MODEL['delta_t'] (100 Myr default) an compute the
            redshift corresponding to edges of these bins.

        """
        return get_redshift_bin_edges(self.MODEL["delta_t"])

    @property
    def centers_redshift_bins(self):
        """Compute redshift bin centers.

        Returns
        -------
        array floats
            We devide the cosmic time history of the Universe in equally spaced
            bins of cosmic time of self.MODEL['delta_t'] (100 Myr default) an compute the
            redshift corresponding to center of these bins.

        """
        return get_redshift_bin_centers(self.MODEL["delta_t"])
