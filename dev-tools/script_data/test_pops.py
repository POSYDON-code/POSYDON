import os
import shutil
import tracemalloc

import pandas as pd
from binaries_suite import write_binary_to_screen
from tabulate import tabulate

from posydon.config import PATH_TO_POSYDON
from posydon.popsyn.binarypopulation import BinaryPopulation
from posydon.popsyn.synthetic_population import PopulationRunner


def print_pop_settings(population):

    print("Population settings:")

    ignore_kwargs = ["extra_columns", "only_select_columns", "scalar_names",
                     "include_S1", "S1_kwargs", "include_S2", "S2_kwargs",
                     "population_properties", "warnings_verbose", "history_verbose",
                      "error_checking_verbose", "use_MPI", "read_samples_from_file",
                      "RANK", "size", "optimize_ram", "ram_per_cpu",
                      "dump_rate", "tqdm", "temp_directory", "breakdown_to_df"]

    for key, val in population.kwargs.items():
        if key in ignore_kwargs:
            continue
        else:
            print(f"\t {key} : {val}")

def test_popruns():

    path_to_default_params = os.path.join(PATH_TO_POSYDON, "dev-tools/script_data/test_population_params.ini")

    pop = BinaryPopulation.from_ini(path_to_default_params, verbose=False)
    print("Evolving BinaryPopulation...\n")

    print_pop_settings(pop)

    pop.evolve(breakdown_to_df=False, optimize_ram=False, tqdm=False)

    #print(pop.manager.binaries)
    write_binary_to_screen(pop.manager.binaries[0])

    print("Done!")

    #del pop

    #poprun = PopulationRunner(my_ini_filename, verbose=False)
    #print('Number of binary populations:',len(poprun.binary_populations))
    #print('Metallicities:', poprun.solar_metallicities)
    #print('Number of binaries:', poprun.binary_populations[0].number_of_binaries)
    #print("Evolving PopulationRunner...")
    #poprun.evolve(overwrite=True)
    #print("Done!")


if __name__ == "__main__":

    #tracemalloc.start()

    test_popruns()

    #current, peak = tracemalloc.get_traced_memory()
    #print(f"Current: {current/1e6:.1f} MB, Peak: {peak/1e6:.1f} MB")
