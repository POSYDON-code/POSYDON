"""Counting distincts occurence of binary simulation from a file"""

__authors__ = [
    "Wène Kouarfate <Wene.Kouarfate@etu.unige.ch>"
]

import numpy as np
import pandas as pd
import os
from collections import Counter

class ParseDataFrame:
    """Handle the binary parsing"""
    def __init__(self,
                 filename,
                 path="./",
                 key = 'history',
                 column_list = ['state','event','S1_state','S2_state'], 
                 index_name = 'binary_index',
                 start=None,
                 stop=None,
                 chunk_size = 500000):
        
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
        """function to be given as key argument to DataFrameGroupBy.apply()"""
        h = hash(tuple(df_gb.to_numpy().ravel()))#.to_numpy() recommanded by pandas doc instead of .values
        
        self.counts[h] += 1
        self.index_list.setdefault(h, df_gb.index[0])
        
        return None


    def _parse_dataf_groupby_col3_chunk(self):
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
        total = sum(self.counts.values())
        return Counter({self.index_list[k]: 100 * self.counts[k] / total
                          for k in self.counts.keys()})
    
    def get_most_numpy(self, k):
        #one can then acess columns for VHDiagramm_m
        return np.array(self.count_dict.most_common(k))
    
    def parse_dataf_gb_iter_chunk(dataf, index_list, cnt):
        """a more relevant parser imo but turns out to take more time than than groupby/apply"""
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