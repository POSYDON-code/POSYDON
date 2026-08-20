"""Lazy wrapper around HDF5 datasets with optional dtype conversion."""

__authors__ = [
    "Seth Gossage <seth.gossage@northwestern.edu>"
]

import numpy as np
import pandas as pd


class LazyHDF5:
    """
    Lazy wrapper around an HDF5 dataset with optional dtype conversion.

    This class provides a lightweight interface for accessing data from an
    HDF5 dataset without immediately loading the entire dataset into memory.
    Data are retrieved lazily when indexed. Optionally, a set of dtype
    conversions can be applied when data are accessed.

    If dtype mappings are provided, retrieved data are cast to the specified
    dtypes either per-field (for structured arrays) or for the selected field
    when accessed by name.

    Assignments (via __setitem__) trigger full materialization of the dataset
    in memory, after which the internal storage is replaced by the in-memory
    array.

    Parameters
    ----------
    dataset : h5py.Dataset or array-like
        The underlying dataset providing the data. Typically an HDF5 dataset
        object supporting NumPy-style indexing.
    dtype_set : dict, optional
        Mapping of field names to NumPy dtypes used to cast the returned data.
        This is typically used for structured arrays where individual fields
        require specific dtype conversions.

    Notes
    -----
    - Data are only read from the dataset when accessed via ``__getitem__`` or
      when converted to a NumPy array.
    - Writing via ``__setitem__`` loads the entire dataset into memory before
      modifying it.
    - The ``dtype`` property reflects the converted dtype if ``dtype_set`` is
      provided.
    """
    def __init__(self, dataset, dtype_set=None):
        """Initialize the lazy wrapper around an HDF5 dataset.

        Parameters
        ----------
        dataset : h5py.Dataset or array-like
            The underlying dataset providing the data.
        dtype_set : dict or None (default: None)
            Mapping of field names to NumPy dtypes used to cast the data.

        """
        self._dataset = dataset
        self._dtype_set = dtype_set
        if self._dtype_set is not None:
            self._dtype_list = list(self._dtype_set.items())

    def __getitem__(self, idx):
        data = self._dataset[idx]
        if self._dtype_set is not None:
            if isinstance(idx, str):
                data = data.astype(self._dtype_set[idx])
            else:
                data = data.astype(self._dtype_list)
        return data

    def __setitem__(self, idx, value):
        # materialize full array in memory
        arr = self.__array__()
        # write new value
        arr[idx] = value

        self._dataset = arr


    def __array__(self):
        data = self._dataset[()]
        if self._dtype_set is not None:
            data = data.astype(self._dtype_list)
        return data

    def astype(self, dtype): # pragma: no cover
        """Return a new LazyHDF5 with the data cast to the given dtype.

        Parameters
        ----------
        dtype : numpy dtype or str
            The dtype to cast the data to.

        Returns
        -------
        LazyHDF5
            A new LazyHDF5 wrapping the cast array.

        """
        return LazyHDF5(np.asarray(self).astype(dtype), self._dtype_set)

    @property
    def dtype(self):
        """Return the dtype of the wrapped data.

        Returns
        -------
        numpy dtype
            The dtype of the dataset, or the converted dtype if a
            `dtype_set` was provided.

        """
        if self._dtype_set is not None:
            return np.dtype(self._dtype_list)
        return self._dataset.dtype

    @property
    def shape(self): # pragma: no cover
        """Return the shape of the wrapped dataset.

        Returns
        -------
        tuple
            The shape of the underlying dataset.

        """
        return self._dataset.shape

    def __len__(self): # pragma: no cover
        return len(self._dataset)

    def to_df(self): # pragma: no cover
        """Convert the wrapped data to a pandas DataFrame.

        Returns
        -------
        pandas.DataFrame
            The data as a DataFrame.

        """
        return pd.DataFrame(self.__array__())
