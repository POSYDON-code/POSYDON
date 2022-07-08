"""The definition of the SimulationModel class."""


__authors__ = [
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Maxime Rambosson <Maxime.Rambosson@etu.unige.ch>",
]


import pandas as pd
import os


class SimulationModel:
    """Model containing the dataset of simulations."""

    def __init__(self, filename, path="./"):
        """Initialize a SimulationModel instance."""
        self._df = None
        self.path = path
        self.filename = filename

    def load_csv(self):
        """Load dataframe as CSV with .gz compression."""
        self._df = pd.read_hdf(
            os.path.join(self.path, self.filename), compression="gzip",
            key='history', low_memory=False
        )

    def get_by_binary_index(self, index):
        """Return a copy of dataframe's slice with binary_index == index.

        Parameters
        ----------
        index : int
            Binary index wanted.

        Returns
        -------
        pandas.DataFrame
            Copy of a sub-dataframe with binary_index == index

        """
        return self._df.loc[index].copy()
