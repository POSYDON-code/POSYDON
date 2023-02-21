"""Evolve multiple BinaryPopulations together.

e.g. with multiple metallicities
"""



__authors__ = [
    "Kyle Akira Rocha <kylerocha2024@u.northwestern.edu>",
    "Simone Bavera <Simone.Bavera@unige.ch>",
]


from posydon.popsyn.io import parse_inifile, binarypop_kwargs_from_ini
from posydon.popsyn.binarypopulation import BinaryPopulation


class SyntheticPopulation:

    def __init__(self, path):
        """
        Parameters
        ----------

        path : str
            Path to the inifile to parse. You can supply a list in the
            metallicity parameter to evolve more than one population.

        """

        synthetic_pop_params = binarypop_kwargs_from_ini(path)

        self.metallicity = synthetic_pop_params['metallicity']

        if not isinstance( self.metallicity, list):
            self.metallicity = [self.metallicity]

        self.binary_populations = []
        for met in self.metallicity:
            ini_kw = binarypop_kwargs_from_ini(path)
            ini_kw['metallicity'] = met
            ini_kw['temp_directory'] = self.create_met_prefix(met) + ini_kw['temp_directory']
            self.binary_populations.append(  BinaryPopulation(**ini_kw)  )


    def evolve(self):
        """Evolve population(s) at given Z(s)."""
        for ind, pop in enumerate( self.binary_populations ):
            print( f'Z={pop.kwargs["metallicity"]:.2e} Z_sun' )
            pop.evolve()
            met_prefix = f'{pop.kwargs["metallicity"]:.2e}_Zsun_'
            pop.save( met_prefix + 'population.h5' )

    @staticmethod
    def create_met_prefix(met):
        """Append a prefix to the name of directories for batch saving."""
        return f'{met:.2e}_Zsun_'
