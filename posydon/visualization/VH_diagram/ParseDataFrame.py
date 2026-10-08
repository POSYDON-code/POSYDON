"""Counting distincts occurence of binary simulation from a file."""

__authors__ = [
    "Wène Kouarfate <Wene.Kouarfate@etu.unige.ch>"
]

import os
from collections import Counter

import numpy as np
import pandas as pd


class ParseDataFrame:
    """Handle the binary parsing."""
    def __init__(self,
                 filename,
                 path="./",
                 key = 'history',
                 column_list = ['state','event','S1_state','S2_state'],
                 index_name = 'binary_index',
                 start=None,
                 stop=None,
                 chunk_size = 500000):
        """Initialize a ParseDataFrame instance.

        Parameters
        ----------
        filename : str
            Name of the file to parse.
        path : str, optional
            Path to the file.
        key : str, optional
            Key under which the dataframe is stored in the HDF file.
        column_list : list of str, optional
            Columns to read from the file.
        index_name : str, optional
            Name of the index column.
        start : int, optional
            Row number to start reading from.
        stop : int, optional
            Row number to stop reading at.
        chunk_size : int, optional
            Number of rows read per chunk.

        """
        self.path = path
        self.filename = filename
        self.key = key
        self.column_list = column_list
        self.index_name = index_name
        self.chunk_size = chunk_size
        self.start = start
        self.stop = stop

        self.index_list = dict()
        self.counts = Counter()

        self._parse_dataf_groupby_col3_chunk()
        self.count_dict = Counter({self.index_list[k]: self.counts[k]
                          for k in self.counts.keys()})


    def _f_lambda(self, df_gb):
        """Function to be given as key argument to DataFrameGroupBy.apply().

        Parameters
        ----------
        df_gb : pandas.DataFrame
            Group of the dataframe to process.

        """
        h = hash(tuple(df_gb.to_numpy().ravel()))#.to_numpy() recommanded by pandas doc instead of .values

        self.counts[h] += 1
        self.index_list.setdefault(h, df_gb.index[0])

        return None


    def _parse_dataf_groupby_col3_chunk(self):
        """Parse the dataframe by grouping on the third column, chunk by chunk."""
        file_path = os.path.join(self.path, self.filename)
        rdf = pd.DataFrame(columns = self.column_list, index=pd.Index([], name=self.index_name))

        for dataf in pd.read_hdf(file_path, self.key,
                                 columns=self.column_list,
                                 start = self.start,
                                 stop = self.stop,
                                 chunksize=self.chunk_size):

            dataf = pd.concat([rdf, dataf])
            rdf = dataf.loc[[dataf.index[-1]]]
            dataf = dataf.drop(dataf.index[-1])

            gb_df_col = dataf.groupby(by=dataf.index.name)
            gb_df_col.apply(self._f_lambda)

        self._f_lambda(rdf)

    def get_frequencies(self):
        """Return the frequencies of each binary as a percentage of the total.

        Returns
        -------
        collections.Counter
            Counter mapping each binary index to its frequency in percent.

        """
        total = sum(self.counts.values())
        return Counter({self.index_list[k]: 100 * self.counts[k] / total
                          for k in self.counts.keys()})

    def get_most_numpy(self, k):
        """Return the k most common binaries as a numpy array.

        Parameters
        ----------
        k : int
            Number of most common binaries to return.

        Returns
        -------
        numpy.ndarray
            Array of the k most common binaries.

        """
        #one can then acess columns for VHDiagramm_m
        return np.array(self.count_dict.most_common(k))

    def parse_dataf_gb_iter_chunk(dataf, index_list, cnt):
        """Parse the dataframe by iterating over the groups chunk by chunk.

        A more relevant parser imo but turns out to take more time than
        groupby/apply.

        Parameters
        ----------
        dataf : pandas.DataFrame
            Dataframe to parse.
        index_list : dict
            Mapping from hash to binary index.
        cnt : collections.Counter
            Counter of occurences of each binary.

        """
        file_path = os.path.join(self.path, self.filename)
        rdf = pd.DataFrame(columns = self.column_list, index=pd.Index([], name=self.index_name))

        for dataf in pd.read_hdf(file_path, self.key,
                                     columns=self.column_list,
                                     start = self.start,
                                     stop = self.stop,
                                     chunksize=self.chunk_size):

            dataf = pd.concat([rdf, dataf])
            rdf = dataf.loc[[dataf.index[-1]]]
            dataf = dataf.drop(dataf.index[-1])

            gb_df_col = dataf.groupby(dataf.index.name)
            for i,s in gb_df_col.__iter__():
                h = hash(tuple(df_gb.to_numpy().ravel()))
                self.counts[h] += 1
                self.index_list.setdefault(h, df_gb.index[0])
