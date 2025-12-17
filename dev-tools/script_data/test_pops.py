import os
import shutil
import tracemalloc

from posydon.config import PATH_TO_POSYDON
from posydon.popsyn.binarypopulation import BinaryPopulation
from posydon.popsyn.synthetic_population import PopulationRunner


def test_popruns():

    #path_to_default_params = os.path.join(PATH_TO_POSYDON, "posydon/popsyn/population_params_default.ini")
    my_ini_filename = "population_params.ini"

    # Copy the defult .ini file into the current directory with new .ini filename defined above
    #my_ini_path = shutil.copyfile(path_to_default_params, my_ini_filename)

    pop = BinaryPopulation.from_ini(my_ini_filename, verbose=True)
    print("Evolving BinaryPopulation...")

    pop.evolve(optimize_ram=False, tqdm=False)
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
